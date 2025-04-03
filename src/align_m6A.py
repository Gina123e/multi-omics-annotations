'''
This script aims to align m6A sits onto transcripts.
@author: ljf
@affliation: zhejiang University
@date:2025.0328
'''
import pandas as pd
import numpy as np
import re
from tqdm import tqdm
def is_m6A_motif(sequence):
    # 定义 m6A motif 模式 (DRACH)
    motif_pattern = r'^[GAT][AG]A[C][ACT]$'
    # 使用正则表达式匹配序列
    return bool(re.match(motif_pattern, sequence))

fasta = pd.read_csv('../../../../SNP/gencode.v38_tx.fa', sep=',',header=None)
fasta.columns = ['id', 'fa']

#read m6A_intersect_genome.gtf
GTF = pd.read_csv('../../../../SNP/gencode.v38.annotation.gtf', header=None, sep='\t', comment='#')
GTF.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
GTF_tx = GTF[GTF['feature'] == 'transcript']
GTF_tx['transcript_id'] = GTF_tx['attribute'].str.extract(r'transcript_id "([^"]+)"')

GTF_exon = GTF[GTF['feature'] == 'exon']
GTF_exon['transcript_id'] = GTF_exon['attribute'].str.extract(r'transcript_id "([^"]+)"')
GTF_exon['exon_number'] = GTF_exon['attribute'].str.extract(r'exon_number (\d+)')
GTF_exon['exon_id'] = GTF_exon['attribute'].str.extract(r'exon_id "([^"]+)"')

m6A_intersect_genome = pd.read_csv('m6A_1_intersect.gtf', usecols=[0,1,5,14],header=None, sep='\t')

m6A_intersect_genome.columns = ['chr', 'start',  'strand','info']
m6A_intersect_genome['transcript_id'] = m6A_intersect_genome['info'].str.extract(r'transcript_id "([^"]+)"')
m6A_intersect_genome['exon_number'] = m6A_intersect_genome['info'].str.extract(r'exon_number (\d+)')

f1 = open('./m6A_on_tx_classic.txt', 'w')
line = ['id', 'pos', 'm6A_chr', 'm6A_strand', 'm6A_pos','motif']
f2 = open('./m6A_on_tx_nonclassic.txt', 'w')

f1.write('\t'.join(map(str,line)) +'\n')
f2.write('\t'.join(map(str,line)) +'\n')
for i in tqdm(range(m6A_intersect_genome.shape[0])):
    tx = m6A_intersect_genome.iloc[i]['transcript_id']
    m6A_pos = m6A_intersect_genome.iloc[i]['start']
    m6A_chr =  m6A_intersect_genome.iloc[i]['chr']
    exon_number = int(m6A_intersect_genome.iloc[i]['exon_number'])
    m6A_strand = m6A_intersect_genome.iloc[i]['strand']
    rna_seq = fasta[fasta['id'] == tx]

    exons = GTF_exon[GTF_exon['transcript_id'] == tx]
    exons_start, exons_end = list(exons['start']), list(exons['end'])
    pos = 0
    if rna_seq.shape[0]:
        rna_fa = rna_seq.iloc[0]['fa']
        for e in range(exon_number - 1):
            pos += int(exons_end[e]) - int(exons_start[e]) + 1
        if m6A_strand =='+':
            pos += m6A_pos - exons_start[exon_number-1] + 1
        else:
            pos += exons_end[exon_number-1] - m6A_pos + 1
        
        motif = rna_fa[pos-3:pos+2]
        line = [tx, pos, m6A_chr, m6A_strand, m6A_pos,motif]
        if is_m6A_motif(motif):
            #line = [tx, m6A_strand, pos]
            f1.write('\t'.join(map(str,line)) +'\n')
        else:#不是m6A motif打印出来
            #line = [tx, pos, m6A_chr, m6A_strand, m6A_pos,motif]
            f2.write('\t'.join(map(str,line)) + '\n')
f1.close()
f2.close()





# ## fasta
# hg38fa = pd.read_csv('../../../../SNP/gencode.v38_tx_strand.csv', header=0, sep='\t')

# ## 读取transcript coords
# transcript_coords = pd.read_csv('../../../../ljf_projects/1.multi_omics_annotation/genecode/transcript_coord.txt', header=None, sep= '\t')
# transcript_coords.columns = ['tx', 'chr',  'strand', 'coord']
# coordsfile = 'HEK293T_on_tx.txt'
# f = open(coordsfile, 'w')
# for i in tqdm(range(m6A_intersect_tx.shape[0])):
# #for i in range(0):
#     aligned_tx = m6A_intersect_tx.iloc[i]['tx_id']
#     RNA =  hg38fa[hg38fa['tx'] == aligned_tx].iloc[0]['fa'].strip()
#     m6A = m6A_intersect_tx.iloc[i]['start']  #0-based
#     tx_strand = m6A_intersect_tx.iloc[i]['strand']
#     if tx_strand == '-':
#         RNA = RNA[::-1]
    
#     aligned_tx_coords = [int(each) for each in transcript_coords[transcript_coords['tx'] == aligned_tx]['coord'].iloc[0].split(',')]
#     tx_length = sum(1 + np.array(aligned_tx_coords[1::2]) - np.array(aligned_tx_coords[::2]))
#     if tx_length != len(RNA):
#         print('coords error')
#     m6A_relative_coord = 0
    
#     for e in range( int(len(aligned_tx_coords)/2)):
#         if m6A >= aligned_tx_coords[2*e] and m6A <= aligned_tx_coords[2*e+1]:
#             for h in range(e):
#                 m6A_relative_coord += aligned_tx_coords[2*h+1] -  aligned_tx_coords[2*h] + 1
#             m6A_relative_coord += m6A  - aligned_tx_coords[2*e] + 1
#             #print('start2: '+ str(e))
#             break
#     if m6A_relative_coord > 0:
#         if tx_strand == '+':
#             if RNA[m6A_relative_coord-1] == 'A':
#                 f.write( aligned_tx + '\t' + str(tx_strand) +'\t' + str(m6A_relative_coord) + '\n')
#         else:

#             true_m6A = tx_length - m6A_relative_coord + 1
#             if RNA[true_m6A-1] == 'A':
#                 f.write( aligned_tx + '\t'+ str(tx_strand) +'\t' + str(true_m6A) + '\n')
# f.close()