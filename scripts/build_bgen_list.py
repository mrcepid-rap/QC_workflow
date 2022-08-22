import re
import csv
import dxpy

found_files = dxpy.find_data_objects(classname='file', folder='/filtered_bgen/',)

bgen_index = {}

for file in found_files:
    dxfile = dxpy.DXFile(dxid=file['id'], project=file['project'])

    file_match = re.match('chr([\\dXY]{1,2})\\.filtered\\.(\\S+)', dxfile.describe()['name'])
    if file_match is not None:
        chrom = file_match.group(1)
        suffix = file_match.group(2)

        if suffix == 'bgen':
            file_type = 'bgen_dxid'
        elif suffix == 'bgen.bgi':
            file_type = 'bgen_index_dxid'
        elif suffix == 'vep.tsv.gz':
            file_type = 'vep_dxid'
        elif suffix == 'sample':
            file_type = 'sample_dxid'
        elif suffix == 'vep.tsv.gz.tbi':
            file_type = 'vep_index_dxid'

        if chrom in bgen_index:
            bgen_index[chrom][file_type] = file['id']
        else:
            bgen_index[chrom] = {'chrom': chrom, file_type: file['id']}

bgen_index_file = open('bgen_locs.tsv', 'w')
bgen_writer = csv.DictWriter(bgen_index_file, fieldnames=['chrom', 'bgen_dxid', 'bgen_index_dxid', 'sample_dxid', 'vep_dxid', 'vep_index_dxid'], delimiter="\t")
bgen_writer.writeheader()

for key in bgen_index.keys():
    bgen_writer.writerow(bgen_index[key])

bgen_index_file.close()