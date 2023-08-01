#!/usr/bin/env python3

import sys
import csv
import gzip
import dxpy
import tarfile
import sqlite3
import pandas as pd

from pathlib import Path
from typing import TypedDict
from general_utilities.association_resources import run_cmd
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger


class ConsequenceSeverity(TypedDict):
    score: int
    type: str


def define_score(csqs: str) -> ConsequenceSeverity:
    csqs_split = csqs.split('&')
    # This is a Eugene-decided level of severity partially decided based on VEP severity score.
    # This should contain all possible CSQ annotations for a variant other than those reserved
    # for large structural variants
    vep_consequences = {'stop_gained': {'score': 1, 'type': 'PTV'},
                        'frameshift_variant': {'score': 2, 'type': 'PTV'},
                        'splice_acceptor_variant': {'score': 3, 'type': 'PTV'},
                        'splice_donor_variant': {'score': 4, 'type': 'PTV'},
                        'stop_lost': {'score': 5, 'type': 'STOP_LOST'},
                        'start_lost': {'score': 6, 'type': 'START_LOST'},
                        'inframe_insertion': {'score': 7, 'type': 'INFRAME'},
                        'inframe_deletion': {'score': 8, 'type': 'INFRAME'},
                        'missense_variant': {'score': 9, 'type': 'MISSENSE'},
                        'protein_altering_variant': {'score': 10, 'type': 'INFRAME'},
                        'splice_region_variant': {'score': 11, 'type': 'NONCODING'},
                        'incomplete_terminal_codon_variant': {'score': 12, 'type': 'INFRAME'},
                        'start_retained_variant': {'score': 13, 'type': 'SYN'},
                        'stop_retained_variant': {'score': 14, 'type': 'SYN'},
                        'synonymous_variant': {'score': 15, 'type': 'SYN'},
                        '5_prime_UTR_variant': {'score': 16, 'type': 'UTR'},
                        '3_prime_UTR_variant': {'score': 17, 'type': 'UTR'},
                        'intron_variant': {'score': 18, 'type': 'INTRONIC'},
                        'upstream_gene_variant': {'score': 19, 'type': 'UPSTREAM'},
                        'downstream_gene_variant': {'score': 20, 'type': 'DOWNSTREAM'},
                        'intergenic_variant': {'score': 21, 'type': 'INTERGENIC'},
                        'no_score': {'score': 22, 'type': 'ERROR'}}
    ret_csq = vep_consequences['no_score']
    for c in csqs_split:
        if c in vep_consequences:
            if vep_consequences[c]['score'] < ret_csq['score']:
                ret_csq = vep_consequences[c]
    return ret_csq


def final_process_record(rec: dict, severity: ConsequenceSeverity) -> dict:
    # Has to have a "#" to be compatible with VCF I/O
    rec['#CHROM'] = rec['CHROM']
    # Rename some columns for printing purposes
    rec['parsed_csq'] = severity['type']
    # Records who do not have equal REF/ALT length are assigned as InDels
    if len(rec['REF']) != len(rec['ALT']):
        rec['is_indel'] = True
    else:
        rec['is_indel'] = False
    # This corrects an issue with sites with 100% missingness that BCFtools doesn't handle correctly
    if rec['AF'] == '.':
        rec['AF'] = 0
        rec['FILTER'] = 'FAIL'
    else:
        rec['AF'] = float(rec['AF'])
    # Set sensible defaults for a variety of fields:
    # "gnomad_maf", "REVEL" "SIFT", "PolyPhen", "LoF"
    rec['REVEL'] = rec['REVEL'] if rec['REVEL'] != '.' else 'NaN'  # NaN is default VCF spec for missing floats
    rec['SIFT'] = rec['SIFT'] if rec['SIFT'] != '.' else 'NA'  # NA is default VCF spec for missing strings
    rec['PolyPhen'] = rec['PolyPhen'] if rec['PolyPhen'] != '.' else 'NA'  # NA is default VCF spec for missing strings
    rec['LoF'] = rec['LoF'] if rec['LoF'] != '.' else 'NA'  # NA is default VCF spec for missing strings
    rec['AA'] = rec['AA'] if rec['AA'] != '.' else 'NA'  # NA is default VCF spec for missing strings
    rec['AApos'] = rec['AApos'] if rec['AApos'] != '.' else 'NA'  # NA is default VCF spec for missing strings
    return rec


logger = MRCLogger(__name__).get_logger()

# Handle inputs
chromosome = sys.argv[1]
data_dir = f'/{sys.argv[2]}'
info_path = Path(sys.argv[3])
is_hg19 = sys.argv[4].lower()
if is_hg19 == 'false':
    is_hg19 = False
elif is_hg19 == 'true':
    is_hg19 = True
else:
    raise ValueError(f'is_hg19 value of "{is_hg19}" does not parse to True / False')

# The below checks for the various annotation(s) that we use for each imputation version (because they couldn't make it
# easy and do one version for all of them...)
# GEL
if info_path.match('*.tar.gz'):
    tar = tarfile.open(info_path, "r:gz")
    tar.extractall()
    db_file = [Path(x) for x in Path('./').glob('*.imputed.BOLT.stats.tsv.gz')]  # Lazily find the results TSV and load it
    maf_info_db = pd.read_csv(db_file[0], sep="\t", dtype={'CHR': str})
    maf_info_db = maf_info_db[maf_info_db['CHR'] == chromosome]
    maf_info_db.rename(columns={'BOLT_MAF': 'MAF'})

    # Find the bgen.bgi file that we want to annotate:
    logger.info('Locating .bgi file...')
    search_regex = f'ukb\\d+_c{chromosome}_b\\d+_v\\d+.bgen.bgi'
    bgi_search = dxpy.find_one_data_object(classname="file", name=search_regex, name_mode="regexp",
                                           folder=data_dir, project=dxpy.PROJECT_CONTEXT_ID)

    bgi_dxfile = dxpy.DXFile(dxid=bgi_search['id'], project=bgi_search['project'])
    bgi_file = Path(bgi_dxfile.describe(fields={'name': True})['name'])
    dxpy.download_dxfile(bgi_dxfile.get_id(), bgi_file.name)
    file_stem_name = bgi_file.name.replace('.bgen.bgi', '')

    # bgen to vcf – iterate through the bgen index using sqlite3
    bgi_connection = sqlite3.connect(bgi_file)
    # Variant table (as a sqlite3 db) has columns:
    # ['chromosome', 'position', 'rsid', 'number_of_alleles', 'allele1', 'allele2',
    # 'file_start_position', 'size_in_bytes']

    # Allows me to programmatically access them after running a sqlite3 query
    logger.info('Merging .bgi and INFO scores...')
    bgi_connection.row_factory = sqlite3.Row
    bgi_table = pd.read_sql('SELECT * FROM Variant', con=bgi_connection)

    # And merge to the maf/info db:
    # BOLT ALLELE1/ALLELE0 is flipped...
    var_table = maf_info_db.merge(bgi_table,
                                  left_on=['varID', 'ALLELE1', 'ALLELE0'],
                                  right_on=['rsid', 'allele1', 'allele2'],
                                  how='left')

    # Keep columns we need
    var_table = var_table[['rsid', 'position', 'allele1', 'allele2', 'MAF', 'INFO']]

# Original
elif info_path.match('*.mfi.txt'):
    mfi_search = dxpy.find_one_data_object(classname="file", name=f'{info_path}', name_mode="exact",
                                           folder=data_dir, project=dxpy.PROJECT_CONTEXT_ID)
    mfi_dxfile = dxpy.DXFile(dxid=mfi_search['id'], project=mfi_search['project'])
    mfi_file = Path(mfi_dxfile.describe(fields={'name': True})['name'])
    dxpy.download_dxfile(mfi_dxfile.get_id(), mfi_file.name)
    file_stem_name = mfi_file.name.replace('.mfi.txt', '')
    var_table = pd.read_csv(mfi_file, sep="\t", names=['varID', 'rsid', 'position', 'allele1',
                                                       'allele2', 'MAF', 'minor_allele', 'INFO'])

    var_table = var_table[['rsid', 'position', 'allele1', 'allele2', 'MAF', 'INFO']]

# TOPMed
elif info_path.match('*.vcf.gz'):
    vcf_search = dxpy.find_one_data_object(classname="file", name=f'{info_path}', name_mode="exact",
                                           folder=f'{data_dir}/helper_files/', project=dxpy.PROJECT_CONTEXT_ID)
    vcf_dxfile = dxpy.DXFile(dxid=vcf_search['id'], project=vcf_search['project'])
    vcf_file = Path(vcf_dxfile.describe(fields={'name': True})['name'])
    dxpy.download_dxfile(vcf_dxfile.get_id(), vcf_file.name)
    file_stem_name = vcf_file.name.replace('.sites.vcf.gz', '')
    query_cmd = f'bcftools query -f "%ID\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\t%INFO/R2\\n" -o /test/info.tsv ' \
                f'/test/{vcf_file}'
    run_cmd(query_cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
    var_table = pd.read_csv('info.tsv', sep="\t", names=['rsid', 'position', 'allele1',
                                                         'allele2', 'MAF', 'INFO'])

else:
    raise ValueError(f'File {info_path} does not appear to be a tar of results, an .mfi.txt, or .vcf.gz file')

# # Process downloaded data:
logger.info('Unpacking resource files...')
cmd = 'gunzip hs38DH.fa.gz'
run_cmd(cmd)

# Loftee
cmd = "tar -zxf loftee_hg38.tar.gz"
run_cmd(cmd)
Path('loftee_hg38.tar.gz').unlink()

# VEP
Path("vep_caches/").mkdir()  # This is for legacy reasons to make sure all tests work...
cmd = "tar -zxf homo_sapiens_vep_108_GRCh38.tar.gz -C vep_caches/"
run_cmd(cmd)
Path('homo_sapiens_vep_108_GRCh38.tar.gz').unlink()

# Setup conversion to Hg38 is hg19, if required
if is_hg19:
    run_cmd('pip3 install pyBigWig==0.3.18', is_docker=False)  # Current pyBigWig install is borked, install 0.3.18
    run_cmd('pip3 install CrossMap', is_docker=False)
    run_cmd('wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz', is_docker=False)

logger.info('Generating split VCF files for VEP annotation...')
variant_count = 0
file_num = 0
total_variant_count = len(var_table)

# Filter to variants we actually want
var_table = var_table.query('MAF >= 9E-4 & INFO >= 0.5')

# Open a blank '0' dummy file to make iteration easier
vcf_writer = gzip.open(Path(f'{file_stem_name}_chunk{file_num}.vcf.gz'), 'wt')
for var in var_table.itertuples():
    if variant_count % 10000 == 0:
        vcf_writer.close()
        # if hg19, we need to CrossMap
        if is_hg19 and variant_count != 0:
            run_cmd(f'CrossMap.py gvcf --compress hg19ToHg38.over.chain.gz {file_stem_name}_chunk{file_num}.vcf.gz hs38DH.fa '
                    f'{file_stem_name}_chunk{file_num}.crossMap.vcf', is_docker=False)
            run_cmd(f'bcftools sort -o /test/{file_stem_name}_chunk{file_num}.crossMap.sorted.vcf.gz '
                    f'/test/{file_stem_name}_chunk{file_num}.crossMap.vcf.gz',
                    is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
            Path(f'{file_stem_name}_chunk{file_num}.crossMap.vcf.gz').unlink()
            Path(f'{file_stem_name}_chunk{file_num}.vcf.gz').replace(f'{file_stem_name}_chunk{file_num}.orig.vcf.gz')
            Path(f'{file_stem_name}_chunk{file_num}.crossMap.sorted.vcf.gz').replace(f'{file_stem_name}_chunk{file_num}.vcf.gz')

        file_num += 1
        vcf_writer = gzip.open(Path(f'{file_stem_name}_chunk{file_num}.vcf.gz'), 'wt')
        _ = vcf_writer.write(f'##fileformat=VCFv4.2\n')
        _ = vcf_writer.write(f'##FILTER=<ID=PASS,Description="All filters passed">\n')
        _ = vcf_writer.write(f'##contig=<ID={chromosome}>\n')
        _ = vcf_writer.write(f'##INFO=<ID=MAF,Number=1,Type=Float,Description="Allele Frequency estimate for each alternate allele">\n')
        _ = vcf_writer.write(f'##INFO=<ID=INFO,Number=1,Type=Float,Description="Imputation INFO score">\n')
        _ = vcf_writer.write(f'##INFO=<ID=RSID,Number=1,Type=String,Description="Variant RSID">\n')
        _ = vcf_writer.write(f'##INFO=<ID=A1,Number=1,Type=String,Description="Original A1">\n')
        _ = vcf_writer.write(f'##INFO=<ID=A2,Number=1,Type=String,Description="Original A2">\n')
        _ = vcf_writer.write(f'##INFO=<ID=OPOS,Number=1,Type=Integer,Description="Original position">\n')
        _ = vcf_writer.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    variant_count += 1
    new_id = f'{chromosome}_{var.position}_{var.allele1}_{var.allele2}'
    _ = vcf_writer.write(f'{chromosome}\t{var.position}\t{new_id}\t{var.allele1}\t'
                         f'{var.allele2}\t.\tPASS\tRSID={var.rsid};'
                         f'MAF={var.MAF};INFO={var.INFO};'
                         f'A1={var.allele1};A2={var.allele2};OPOS={var.position}\n')  # to help with hg19 if required

vcf_writer.close()
del var_table

logger.info(f'Printed {variant_count} variants out of {total_variant_count} possible variants ({(variant_count / total_variant_count) * 100:0.2f}%)')


# Run VEP
def run_vep(job_num: int):

    # Run VEP
    vep_cmd = f'perl -Iensembl-vep/cache/Plugins/loftee/ -Iensembl-vep/cache/Plugins/loftee/maxEntScan/ ' \
              f'ensembl-vep/vep --offline --cache --assembly GRCh38 ' \
              f'--dir_cache /test/vep_caches/ --everything --allele_num --distance 250000 ' \
              f'-i /test/{file_stem_name}_chunk{job_num}.vcf.gz --format vcf ' \
              f'--fasta /test/hs38DH.fa ' \
              f'-o /test/{file_stem_name}_chunk{job_num}.vep.vcf.gz ' \
              f'--compress_output bgzip --vcf ' \
              f'--dir_plugins ensembl-vep/cache/Plugins/ ' \
              f'--plugin LoF,loftee_path:ensembl-vep/cache/Plugins/loftee,' \
              f'human_ancestor_fa:/test/loftee_hg38/human_ancestor.fa.gz,' \
              f'conservation_file:/test/loftee_hg38/loftee.sql,' \
              f'gerp_bigwig:/test/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw ' \
              f'--plugin REVEL,/test/new_tabbed_revel_grch38.tsv.gz'
    run_cmd(vep_cmd, True, docker_image='egardner413/mrcepid-burdentesting')

    # Generate VEP TSV
    split_cmd = f'bcftools +split-vep -df "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%FILTER\\t%INFO/MAF\\t%MANE_SELECT\\t%Feature\\t%Gene\\t%BIOTYPE\\t%CANONICAL\\t%SYMBOL\\t%Consequence\\t%REVEL\\t%SIFT\\t%PolyPhen\\t%LoF\\t%Amino_acids\\t%Protein_position\\t%RSID\\t%INFO/INFO\\t%DISTANCE\\t%INFO/A1\\t%INFO/A2\\t%INFO/OPOS\\n" ' \
                f'-o /test/{file_stem_name}_chunk{job_num}.vep_table.tsv ' \
                f'/test/{file_stem_name}_chunk{job_num}.vep.vcf.gz'
    run_cmd(split_cmd, True, docker_image='egardner413/mrcepid-burdentesting')

    # And parse the TSV file:
    with Path(f'{file_stem_name}_chunk{job_num}.vep_table.tsv').open('r') as vep_reader, \
            Path(f'{file_stem_name}_chunk{job_num}.vep_table.annote.tsv').open('w') as annote_file:

        # These are all possible fields from the vep table that we generated in run_vep()
        # And then read it in as a csv.DictReader()
        csv_reader_header = (
            "CHROM", "POS", "REF", "ALT", "ID", "FILTER", "AF",
            "mane_transcript", "ENST_ID", "ENSG_ID", "biotype", "is_canonical", "symbol", "csq",
            "REVEL", "SIFT", "PolyPhen", "LoF", "AA", "AApos", "RSID", "INFO", "DISTANCE", "A1", "A2", "OPOS")
        vep_csv = csv.DictReader(vep_reader, delimiter="\t", fieldnames=csv_reader_header, quoting=csv.QUOTE_NONE)

        # Next, open a file that will contain (in tabix format tsv) the info we want to add back to the vcf
        # And these are all possible output fields that we want
        annote_writer_header = (
            "#CHROM", "POS", "REF", "ALT", "ID", "AF", "mane_transcript", "ENST_ID", "ENSG_ID", "biotype",  # OG Fields
            "symbol", "csq", "REVEL", "SIFT", "PolyPhen", "LoF", "AA", "AApos", "RSID", "INFO", "DISTANCE", # OG Fields
            "parsed_csq", "is_multiallelic", "is_indel")  # New Fields

        if is_hg19:  # Add extra fields for HG38 coordinates if hg19 imputation
            annote_writer_header = annote_writer_header + ('HG38_CHROM', 'HG38_POS', 'HG38_REF', 'HG38_ALT', 'CrossMapStatus')

        annote_writer = csv.DictWriter(annote_file, delimiter="\t", fieldnames=annote_writer_header,
                                       extrasaction='ignore', restval='NA')

        # Now we need to iterate through the .tsv file that we made
        # A single variant can be spread across multiple rows as it can be annotated for multiple transcripts in the same gene
        # The annotations ARE always one after another (all entries for one variant are sequential), so we don't have to
        # worry about the order of the file.
        # But we do have to collect multiple entries for one variant, and then decide which one is the most "important". So we:
        # 1. Iterate through some records until we find a record that is NOT the same ref/alt
        # 2. Decide if the severity of the current record is "worse" than the currently held_rec
        # 3. Write the record (function: final_process_record())
        # 4. Repeat steps 1 - 3 for the next record until we reach the end of the file
        held_rec_name = None
        held_rec = None
        held_severity_score = None
        # Iterate through records (step 1)
        total_recs = 0
        for rec in vep_csv:
            total_recs += 1
            # Set a unique record ID
            current_rec_name = '%s_%s_%s_%s' % (rec['CHROM'], rec['POS'], rec['REF'], rec['ALT'])
            # Make sure distance is actually an integer:
            try:
                rec['DISTANCE'] = int(rec['DISTANCE'])
            except ValueError:  # Sometimes distance is '.' because no gene closer than 250kbp
                rec['DISTANCE'] = 999999
            if current_rec_name != held_rec_name:  # If ID is not the same, write currently held record and reset (steps 3 - 4)
                if (held_rec_name != None):  # Obviously, don't print if going through the first rec since there is no stored INFO yet
                    # Write the record with the most severe consequence (step 3)
                    held_rec = final_process_record(held_rec, held_severity_score)
                    # if hg19, switch fields back to where they're supposed to be
                    if is_hg19:
                        held_rec['HG38_POS'] = held_rec['POS']
                        held_rec['HG38_REF'] = held_rec['REF']
                        held_rec['HG38_ALT'] = held_rec['ALT']
                        held_rec['HG38_CHROM'] = held_rec['#CHROM']
                        held_rec['#CHROM'] = chromosome  # Sometimes crossmap goes to a different chromosome
                        held_rec['POS'] = held_rec['OPOS']
                        held_rec['REF'] = held_rec['A1']
                        held_rec['ALT'] = held_rec['A2']
                        held_rec['CrossMapStatus'] = 'PASS'
                    _ = annote_writer.writerow(held_rec)
                # Reset to a new record (step 4)
                held_rec_name = current_rec_name
                held_rec = rec  # Make sure at least one record makes it through
                # This function decides how "severe" a given CSQ annotation for a record is. See the function for more details
                held_severity_score = define_score(held_rec['csq'])
            else:
                # Calculate severity of this record
                current_severity_score = define_score(rec['csq'])
                # check to see if we should prioritise the new record based on the following ordered criteria (step 2):
                # Below are named in DECREASING selection importance
                # 1. protein_coding transcript
                # 2. MANE Transcript
                # 3. VEP Canonical Transcript
                # 4. CSQ Severity
                # All "else" statements are when two records have identical annotations for the given category above
                if rec['biotype'] == 'protein_coding' and held_rec['biotype'] != 'protein_coding':
                    held_rec = rec
                    held_severity_score = define_score(held_rec['csq'])
                elif rec['biotype'] != 'protein_coding' and held_rec['biotype'] == 'protein_coding':
                    held_rec = held_rec
                else:
                    if rec['mane_transcript'] != '.' and held_rec['mane_transcript'] == '.':
                        held_rec = rec
                        held_severity_score = define_score(held_rec['csq'])
                    elif rec['mane_transcript'] == '.' and held_rec['mane_transcript'] != '.':
                        # This doesn't actually do anything, just to keep things obvious / Python happy
                        held_rec = held_rec
                    else:
                        # Prioritise up/downstream variants by closest (w/in 250kbp)
                        if (current_severity_score['type'] == 'UPSTREAM' or current_severity_score['type'] == 'DOWNSTREAM') \
                                and (held_severity_score['type'] == 'UPSTREAM' or held_severity_score['type'] == 'DOWNSTREAM'):
                            if rec['DISTANCE'] < held_rec['DISTANCE']:
                                held_rec = rec
                                held_severity_score = define_score(held_rec['csq'])
                            elif rec['DISTANCE'] > held_rec['DISTANCE']:
                                # This doesn't actually do anything, just to keep things obvious / Python happy
                                held_rec = held_rec
                        else:
                            if current_severity_score['score'] < held_severity_score['score']:
                                held_rec = rec
                                held_severity_score = define_score(held_rec['csq'])

        if total_recs == 0:
            raise Exception(f'No records were iterated for job {job_num}')

        # And print the last record since it cannot be compared to an old record above:
        held_rec = final_process_record(held_rec, held_severity_score)
        if is_hg19:
            held_rec['HG38_POS'] = held_rec['POS']
            held_rec['HG38_REF'] = held_rec['REF']
            held_rec['HG38_ALT'] = held_rec['ALT']
            held_rec['HG38_CHROM'] = held_rec['#CHROM']
            held_rec['#CHROM'] = chromosome  # Sometimes crossmap goes to a different chromosome
            held_rec['POS'] = held_rec['OPOS']
            held_rec['REF'] = held_rec['A1']
            held_rec['ALT'] = held_rec['A2']
            held_rec['CrossMapStatus'] = 'PASS'
        _ = annote_writer.writerow(held_rec)

        # if hg19, we need to cat variants that didn't liftover on the bottom with NA fields
        if is_hg19 and Path(f'{file_stem_name}_chunk{job_num}.crossMap.vcf.unmap').exists():
            query_cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%FILTER\\t%INFO/MAF\\t%RSID\\t%INFO/INFO\\n" ' \
                        f'-o /test/{file_stem_name}_chunk{job_num}.unmap.tsv ' \
                        f'/test/{file_stem_name}_chunk{job_num}.crossMap.vcf.unmap'
            run_cmd(query_cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
            with Path(f'{file_stem_name}_chunk{job_num}.unmap.tsv').open('r') as unmap_tsv:
                unmap_reader_header = (
                    "CHROM", "POS", "REF", "ALT", "ID", "FILTER", "AF",
                    "RSID", "INFO")
                unmap_csv = csv.DictReader(unmap_tsv, delimiter="\t", fieldnames=unmap_reader_header,
                                           quoting=csv.QUOTE_NONE)
                for var in unmap_csv:
                    var['HG38_CHROM'] = var['CHROM']
                    var['HG38_POS'] = 'NA'
                    var['HG38_REF'] = 'NA'
                    var['HG38_ALT'] = 'NA'
                    var['CrossMapStatus'] = 'FAIL'
                    var['#CHROM'] = chromosome
                    var['AF'] = float(var['AF'])
                    var['REVEL'] = 'NaN'
                    if len(var['REF']) != len(var['ALT']):
                        var['is_indel'] = True
                    else:
                        var['is_indel'] = False
                    _ = annote_writer.writerow(var)
    sort_cmd = f'sort -k1,1 -k2,2n {file_stem_name}_chunk{job_num}.vep_table.annote.tsv > {file_stem_name}_chunk{job_num}.vep_table.annote.sorted.tsv'
    run_cmd(sort_cmd, is_docker=False)


logger.info('Performing VEP annotation...')
thread_utility = ThreadUtility(incrementor=5, thread_factor=1)
for file in range(1, file_num+1):
    thread_utility.launch_job(run_vep,
                              job_num=file)
thread_utility.collect_futures()

# Cat everything together:
logger.info('Concatenating individual VEP runs...')
cmd = f'cat {file_stem_name}_chunk*.vep_table.annote.sorted.tsv > {file_stem_name}.vep_table.annote.sorted.tsv'
run_cmd(cmd, is_docker=False)


def get_manh_pos(chrom, pos, tot_len):
    cum_len = ref_table[ref_table['CHROM'] == chrom]['CUMLEN']
    manh_pos = (pos + cum_len) / tot_len
    return manh_pos


if is_hg19:
    ref_file = Path('hs37d5.fa.fai')
else:
    ref_file = Path('hs38DH.fa.fai')

logger.info('Performing final processing and .tsv generation...')
with Path(f'{file_stem_name}.vep_table.annote.sorted.tsv').open('r') as annote_file, \
        Path(f'{file_stem_name}.annote.tsv').open('w') as final_annotations_file,\
        ref_file.open('r') as ref_index:
    pd_names = ["#CHROM", "POS", "REF", "ALT", "ID", "AF", "mane_transcript", "ENST_ID", "ENSG_ID", "biotype",  # OG Fields
                "symbol", "csq", "REVEL", "SIFT", "PolyPhen", "LoF", "AA", "AApos", "RSID", "INFO", "DISTANCE",  # OG Fields
                "parsed_csq", "is_multiallelic", "is_indel"]
    dtype_dict = {'#CHROM': str}
    if is_hg19:
        pd_names.extend(['HG38_CHROM', 'HG38_POS', 'HG38_REF', 'HG38_ALT', 'CrossMapStatus'])
        dtype_dict['HG38_CHROM'] = str
    annote_table = pd.read_csv(annote_file, sep="\t", dtype=dtype_dict,
                               names=pd_names)
    annote_table = annote_table.rename(columns={'#CHROM': 'CHROM'})
    annote_table = annote_table.sort_values(by='POS')
    # Set manh.pos for plotting
    ref_table = pd.read_csv(ref_index, sep="\t", names=['CHROM', 'LEN', 'CUMLEN', 'WIDTH', 'WIDTH2'], dtype={'CHROM':'str'})
    ref_table['CHROM'] = ref_table['CHROM'].apply(lambda x: x.replace('chr', ''))  # Hg38 'chr' prefix
    tot_len = int(ref_table[ref_table['CHROM'] == 'Y']['LEN'] + ref_table[ref_table['CHROM'] == 'Y']['CUMLEN'])
    annote_table['MANH.POS'] = annote_table.apply(lambda x: get_manh_pos(x["CHROM"], x['POS'], tot_len), axis=1)
    # Some columns get turned into floats for some reason... convert back
    annote_table['DISTANCE'] = annote_table['DISTANCE'].apply(lambda x: f'{x:0.0f}')
    if is_hg19:
        annote_table['HG38_POS'] = annote_table['HG38_POS'].apply(lambda x: f'{x:0.0f}')
        annote_table['HG38_COORD'] = annote_table.apply(lambda x: 'NA' if pd.isna(x['HG38_POS']) else f'{x["HG38_CHROM"]}_{x["HG38_POS"]}', axis=1)
    annote_table.to_csv(final_annotations_file, sep="\t", na_rep='NA', index=False)
    bgzip_cmd = f'bgzip /test/{file_stem_name}.annote.tsv'
    run_cmd(bgzip_cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')
    index_cmd = f'tabix -c C -s 1 -b 2 -e 2 /test/{file_stem_name}.annote.tsv.gz'
    run_cmd(index_cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting')

with Path('output_files.txt').open('w') as output_writer:
    output_writer.write(f'{file_stem_name}.annote.tsv.gz\n')
    output_writer.write(f'{file_stem_name}.annote.tsv.gz.tbi\n')
