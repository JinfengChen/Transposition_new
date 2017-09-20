#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python RunTEMP_sumtable.py --input RIL275_TEMP

Summary mPing insertions and sequence depth for RIL275
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam/GN1.bam       141.39  25.74   76.38   303.56
def lib_size(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line)
                ril  = re.sub(r'\D+', '', os.path.split(unit[0])[1])
                data[ril] = [unit[1], unit[2], unit[3], unit[4]]
    return data

#Sample  #Read   Average Total   Depth   Mapped_Depth    Mapped_rate     #Library        FileName
#GN1_?   23383054        101     2361688454      6.34862487634   6.19321574194   0.975520819479  1       ../../input/fastq/Bam/RIL1_0_CGTACG_FC153L5.recal.bam
def seq_depth(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'):
                unit = re.split(r'\t', line)
                ril  = re.sub(r'\D+', '', unit[0])
                depth   = '%.2f' % float(unit[4])
                mapped  = '%.2f' % float(unit[5])
                maprate = '%.2f' % float(unit[6])
                data[ril] = [depth, mapped, maprate]
    return data

#Chr1	HEG4	TEMP	6244628	6245128	.	-	.	
#ID=TE_Insertion_1;Note=Non-reference;Name=mPing;Class=singleton;
#VariantSupport=1;Frequency=0.1000;Junction1Support=0;Junction2Support=0;Five_Support=0;Three_Support=1;
def parse_gff(infile):
    both_end = 0
    sing_end = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line)
                anno = re.split(r';', unit[8])
                lj   = re.sub(r'Junction1Support=', '', anno[6])
                rj   = re.sub(r'Junction2Support=', '', anno[7])
                ls   = re.sub(r'Five_Support=', '', anno[8])
                rs   = re.sub(r'Three_Support=', '', anno[9])
                if int(lj) + int(ls) > 0 and int(rj) + int(rs) > 0:
                    both_end += 1
                else:
                    sing_end += 1
    return [both_end, sing_end]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   

    depth  = seq_depth('RIL.bam.unique.stat')
    size   = lib_size('RIL275.InsertSize.list') 
    #RIL275_RelocaTEi/RelocaTEi_GN1/repeat/results/ALL.all_nonref_insert.gff 
    #RIL275_TEMP/TEMP_GN1/GN1.insertion.refined.bp.summary.gff
    data   = defaultdict(lambda : list())
    output = '%s.summary.table' %(args.input)
    dirs   = glob.glob('%s/TEMP_*' %(args.input))
    for d in dirs:
        ril    = re.split(r'\_', os.path.split(d)[1])[1]
        ril_id = re.sub(r'\D+', '', ril)
        gff  = '%s/%s.insertion.refined.bp.summary.gff' %(d, ril)
        confident, single_end = parse_gff(gff)
        data[ril_id] = [confident, single_end]

    count = 0
    ofile = open(output, 'w')
    print >> ofile, 'Sample\tInsertSize\tSize_STD\tDepth\tMapped_Depth\tMapped_Rate\tNon_Ref_mPing\tConfident\tEvidence_From_One_End'
    for strain in sorted(data.keys(), key=int):
        non_ref = data[strain][0] + data[strain][1]
        print >> ofile, 'RIL%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(strain, size[strain][0], size[strain][1], depth[strain][0], depth[strain][1], depth[strain][2], non_ref, data[strain][0], data[strain][1])
    ofile.close()


if __name__ == '__main__':
    main()

