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
python RunRelocaTEi_CombinedGFF.py --input RIL275_RelocaTEi

Combine RelocaTEi results of 275 RILs into one gff, which will used to do analysis on unique mPing and excision analysis.
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


#Chr1    not.give        transposable_element_attribute  2129220 2129222 +       .       .       
#ID=Chr1.2129222.spanners;avg_flankers=1;spanners=0;type=homozygous;TE=mping;TSD=TAA
def parse_class_gff(infile):
   total = 0
   hom   = 0
   het   = 0
   som   = 0
   r_hom = re.compile(r'homozygous')
   r_het = re.compile(r'heterozygous')
   r_som = re.compile(r'somatic')
   with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                anno = re.split(r';', unit[8])
                flag = re.sub(r'type=', '', anno[3])
                total += 1
                if r_hom.search(flag):
                    hom += 1
                elif r_het.search(flag):
                    het += 1
                elif r_som.search(flag):
                    som += 1
   return [total, hom, het, som]

#Chr1	not.give	RelocaTE_i	2129220	2129222	.	.	.	
#ID=repeat_Chr1_2129220_2129222;TSD=TAA;Right_junction_reads:1;Left_junction_reads:1;Right_support_reads:0;Left_support_reads:0;
def parse_gff(infile, ril_id, ofile):
    ril_name = 'RIL%s_0' %(ril_id)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                #print '%s\t%s\t%s' %(unit[0], ril_name, unit[8])
                anno = re.split(r';', unit[8])
                unit[1] = ril_name
                anno.insert(1, 'Strain=%s' %(ril_name))
                unit[8] = ';'.join(anno)
                print >> ofile, '\t'.join(unit)
                #rj   = re.sub(r'Right_junction_reads:', '', anno[2])
                #lj   = re.sub(r'Left_junction_reads:', '', anno[3])
                #rs   = re.sub(r'Right_support_reads:', '', anno[4])
                #ls   = re.sub(r'Left_support_reads:', '', anno[5])
                

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
   

    #RIL275_RelocaTEi/RelocaTEi_GN1/repeat/results/ALL.all_nonref_insert.gff 
    data    = defaultdict(lambda : list())
    #output0 = '%s.CombinedGFF.characterized.gff' %(args.input)
    #output1 = '%s.CombinedGFF.ALL.gff' %(args.input)
    output2 = '%s.CombinedGFF.Shared.gff' %(args.input)
    #ofile0  = open(output0, 'w')
    #ofile1  = open(output1, 'w')
    ofile2  = open(output2, 'w')
    dirs    = glob.glob('%s/*_RelocaTEi' %(args.input))
    for d in dirs:
        ril    = re.split(r'\_', os.path.split(d)[1])[0]
        ril_id = re.sub(r'\D+', '', ril)
        gff       = '%s/repeat/results/ALL.all_nonref_insert.gff' %(d)
        gff_class = '%s/repeat/results/ALL.all_nonref_insert.characTErized.gff' %(d)
        gff_shared  = '%s/repeat/results/ALL.all_ref_insert.no_ping_pong.gff' %(d)
        #parse_gff(gff, ril_id, ofile1)
        #parse_gff(gff_class, ril_id, ofile0)
        parse_gff(gff_shared, ril_id, ofile2)
    #ofile0.close()
    #ofile1.close()
    ofile2.close()

if __name__ == '__main__':
    main()

