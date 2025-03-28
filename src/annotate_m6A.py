'''
This script aims to align m6A sits onto transcripts.
@author: ljf
@affliation: zhejiang University
@date:2025.0328
'''
import pandas as pd
import numpy as np
import re

#read m6A_intersect_genome.gtf
m6A_intersect_genome = pd.read_csv('../data/m6A_intersect_genome.gtf', usecols=[0,1,2,5,8,14],header=None, sep='\t')
m6A_intersect_genome.columns = ['chr', 'start', 'end', 'strand','feature','info']
m6A_intersect_tx = m6A_intersect_genome[m6A_intersect_genome['feature'] == 'transcript']

m6A_intersect_tx['tx_id'] = m6A_intersect_tx['info'].str.extract(r'transcript_id "([^"]+)"')
## 读取transcript coords
transcript_coords = pd.read_csv('../genecode/transcript_coord.txt', header=None, sep= '\t')
transcript_coords.columns = ['tx', 'chr',  'strand', 'coord']
coordsfile = '../results/m6A_on_tx.txt'
f = open(coordsfile, 'w')
for i in range(m6A_intersect_tx.shape[0]):
#for i in range(0):
    aligned_tx = m6A_intersect_tx.iloc[i]['tx_id']
    m6A = m6A_intersect_tx.iloc[i]['end']  #0-based
    tx_strand = m6A_intersect_tx.iloc[i]['strand']
    
    aligned_tx_coords = [int(each) for each in transcript_coords[transcript_coords['tx'] == aligned_tx]['coord'].iloc[0].split(',')]
    tx_length = sum(1 + np.array(aligned_tx_coords[1::2]) - np.array(aligned_tx_coords[::2]))
    m6A_relative_coord = 0
    
    for e in range( int(len(aligned_tx_coords)/2)):
        if m6A >= aligned_tx_coords[2*e] and m6A <= aligned_tx_coords[2*e+1]:
            for h in range(e):
                m6A_relative_coord += aligned_tx_coords[2*h+1] -  aligned_tx_coords[2*h] + 1
            m6A_relative_coord += m6A  - aligned_tx_coords[2*e] + 1
            #print('start2: '+ str(e))
            break
    if m6A_relative_coord > 0:
        if tx_strand == '+':
            f.write( aligned_tx + '\t' + str(tx_strand) +'\t' + str(m6A_relative_coord) + '\n')
        else:
            true_m6A = tx_length - m6A_relative_coord + 1
            f.write( aligned_tx + '\t'+ str(tx_strand) +'\t' + str(true_m6A) + '\n')
f.close()


