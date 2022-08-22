import csv
import gzip
from typing import List


def convert_to_list(column: str) -> List[int]:

    split_list = column.split(',')
    split_list = list(map(int, split_list[0:len(split_list) - 1]))
    return split_list


hgReader = csv.DictReader(gzip.open('data_files/hgTables.txt.gz', 'rt'), delimiter="\t")
cds_writer = open('data_files/cds_lengths.txt', 'w')

for row in hgReader:

    if row['transcriptType'] == 'protein_coding':
        gene_ID = row['name'].split('.')[0]

        num_exons = int(row['blockCount'])
        block_sizes = convert_to_list(row['blockSizes'])
        block_starts = convert_to_list(row['chromStarts'])
        gene_start = int(row['#chromStart'])
        cds_start = int(row['thickStart'])
        cds_end = int(row['thickEnd'])

        exon_starts = [(gene_start + start) for start in block_starts]
        exon_ends = []
        for i in range(0, num_exons):
            exon_ends.append(exon_starts[i] + block_sizes[i])

        cds_length = 0
        for i in range(0, num_exons):
            # Do nothing, is non-coding
            if exon_starts[i] < cds_start and exon_ends[i] < cds_start:
                cds_length += 0
            # Do nothing, is non-coding
            elif exon_starts[i] > cds_end and exon_ends[i] > cds_end:
                cds_length += 0
            # Single exon gene
            elif exon_starts[i] < cds_start and exon_ends[i] > cds_end:
                cds_length += cds_end - cds_start
            # CDS + UTR
            elif exon_starts[i] < cds_start < exon_ends[i]:
                cds_length += exon_ends[i] - cds_start
            # CDS + UTR
            elif exon_starts[i] < cds_end < exon_ends[i]:
                cds_length += cds_end - exon_starts[i]
            # Every other exon
            else:
                cds_length += exon_ends[i] - exon_starts[i]

        cds_writer.write(f'{gene_ID}\t{cds_length}\n')

cds_writer.close()




