from tqdm import tqdm 
import numpy as np
import pandas as pd
GTF = pd.read_csv('../genecode/gencode.v38.annotation.gtf', header=None, sep='\t', comment='#')
GTF.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
GTF_tx = GTF[GTF['feature'] == 'transcript']
GTF_tx['transcript_id'] = GTF_tx['attribute'].str.extract(r'transcript_id "([^"]+)"')

GTF_exon = GTF[GTF['feature'] == 'exon']
GTF_exon['transcript_id'] = GTF_exon['attribute'].str.extract(r'transcript_id "([^"]+)"')
GTF_exon['exon_number'] = GTF_exon['attribute'].str.extract(r'exon_number (\d+)')
GTF_exon['exon_id'] = GTF_exon['attribute'].str.extract(r'exon_id "([^"]+)"')

GTF_UTR = GTF[GTF['feature'] == 'UTR']
GTF_UTR['transcript_id'] = GTF_UTR['attribute'].str.extract(r'transcript_id "([^"]+)"')
GTF_UTR['exon_number'] = GTF_UTR['attribute'].str.extract(r'exon_number (\d+)')
GTF_UTR['exon_id'] = GTF_UTR['attribute'].str.extract(r'exon_id "([^"]+)"')


GTF_stop = GTF[GTF['feature'] == 'stop_codon']
GTF_stop['transcript_id'] = GTF_stop['attribute'].str.extract(r'transcript_id "([^"]+)"')
GTF_stop['exon_number'] = GTF_stop['attribute'].str.extract(r'exon_number (\d+)')
GTF_stop['exon_id'] = GTF_stop['attribute'].str.extract(r'exon_id "([^"]+)"')

GTF_start = GTF[GTF['feature'] == 'start_codon']
GTF_start['transcript_id'] = GTF_start['attribute'].str.extract(r'transcript_id "([^"]+)"')
GTF_start['exon_number'] = GTF_start['attribute'].str.extract(r'exon_number (\d+)')
GTF_start['exon_id'] = GTF_start['attribute'].str.extract(r'exon_id "([^"]+)"')

# 用CDS鉴定 UTR在哪儿
GTF_CDS = GTF[GTF['feature'] == 'CDS']
GTF_CDS['transcript_id'] = GTF_CDS['attribute'].str.extract(r'transcript_id "([^"]+)"')
GTF_CDS['exon_number'] = GTF_CDS['attribute'].str.extract(r'exon_number (\d+)')
GTF_CDS['exon_id'] = GTF_CDS['attribute'].str.extract(r'exon_id "([^"]+)"')

from tqdm import tqdm 
UTR5 =[]
UTR3 = []
stop_codon =[]
start_codon = []
f = open('tx_annotation.txt', 'w')
line = ['tx',  'len', 'strand','5UTR', '3UTR', 'start_codon', 'stop_codon','exon_junction']
f.write('\t'.join(map(str, line))  + '\n')
for i in tqdm(range(GTF_tx.shape[0])):
#for i in range(120):
    tx_id = GTF_tx.iloc[i]['transcript_id']
    
    tx_strand = GTF_tx.iloc[i]['strand']
    tx_exons = GTF_exon[GTF_exon['transcript_id'] ==tx_id]
    tx_exons_start, tx_exons_end = list(tx_exons['start']), list(tx_exons['end'])
    #print(tx_exons_start)
    tx_len = 0
    exon_junction = []
    for l in range(len(tx_exons_start)):
        tx_len += tx_exons_end[l] - tx_exons_start[l] + 1
        exon_junction.append(tx_len)
    tx_CDS = GTF_CDS[GTF_CDS['transcript_id'] == tx_id]
    tx_CDS_exons = list(tx_CDS['exon_number'])
    tx_CDS_start, tx_CDS_end = list(tx_CDS['start']), list(tx_CDS['start'])
    if len(tx_CDS_exons):
        five_exon = int(tx_CDS_exons[0])
        three_exon = int(tx_CDS_exons[-1])
            # 转录本的所有外显子
        tx_UTR = GTF_UTR[GTF_UTR['transcript_id'] == tx_id]
        tx_UTR_start, tx_UTR_end = list(tx_UTR['start']), list(tx_UTR['end'])
        tx_UTR_exon_number = [int(each) for each in list(tx_UTR['exon_number'])]

        # stop
        start_codon_start = 0
        stop_codon_start = 0
        five = 0
        three =0
        
        tx_stop = GTF_stop[GTF_stop['transcript_id'] == tx_id]
        tx_start = GTF_start[GTF_start['transcript_id'] == tx_id]
        #tx_start = GTF_start[GTF_start['transcript_id'] == tx_id]
        if tx_start.shape[0]:
            tx_start_exon = int(tx_start.iloc[0]['exon_number'])
            tx_start_start, tx_start_end = tx_start.iloc[0]['start'], tx_start.iloc[0]['end']
        else:
            tx_start_exon = -1
            start_codon_start, start_codon_end = -1, -1 
        if tx_stop.shape[0]:
            tx_stop_exon = int(tx_stop.iloc[0]['exon_number'])
            tx_stop_start, tx_stop_end = tx_stop.iloc[0]['start'], tx_stop.iloc[0]['end']
        else:
            tx_stop_exon = -1
            stop_codon_start, stop_codon_end =  -1, -1

        #stop/start
        if tx_stop_exon > 0:
            for k in range(tx_stop_exon-1):
                stop_codon_start += tx_exons_end[k] - tx_exons_start[k] + 1
                #print(k, stop_codon_start)
            if tx_strand == '+':
                stop_codon_start += tx_stop_start - tx_exons_start[tx_stop_exon-1] + 1
            else:
                stop_codon_start +=  tx_exons_end[tx_stop_exon-1] - tx_stop_end + 1
            stop_codon_end = stop_codon_start + 2
        
        # start codon
        if tx_start_exon > 0:
            for k in range(tx_start_exon-1):
                start_codon_start += tx_exons_end[k] - tx_exons_start[k] + 1
                #print(tx_id, start_codon_start)
            if tx_strand == '+':
                start_codon_start += tx_start_start - tx_exons_start[tx_start_exon-1] + 1
                #print(tx_id, start_codon_start)
            else:
                start_codon_start +=  tx_exons_end[tx_start_exon-1] - tx_start_end + 1
            start_codon_end = start_codon_start + 2
        line = [tx_id, tx_len, tx_strand, five, three, (start_codon_start,start_codon_end), (stop_codon_start,stop_codon_end), exon_junction]
        #print(line)
        f.write('\t'.join(map(str,line)) + '\n')
        
