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
python RunRelocaTEi_sumtable.py --input RIL275_RelocaTEi

Summary mPing insertions and sequence depth for 275 RILs.
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
def parse_class_gff(infile, ping):
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
                tsd  = re.sub(r'TSD=', '', anno[-1])
                mping = '%s.%s' %(unit[0], unit[3])
                if not len(tsd) == 3:
                    continue
                if ping.has_key(mping):
                    continue
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
def parse_gff(infile, ping):
    both_end = 0
    sing_end = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line)
                anno = re.split(r';', unit[8])
                tsd  = re.sub(r'TSD=', '', anno[1])
                mping = '%s.%s' %(unit[0], unit[3])
                rj   = re.sub(r'Right_junction_reads:', '', anno[2])
                lj   = re.sub(r'Left_junction_reads:', '', anno[3])
                rs   = re.sub(r'Right_support_reads:', '', anno[4])
                ls   = re.sub(r'Left_support_reads:', '', anno[5])
                if ping.has_key(mping):
                    continue
                if int(lj) + int(ls) > 0 and int(rj) + int(rs) > 0 and len(tsd) == 3:
                    both_end += 1
                else:
                    sing_end += 1
    return [both_end, sing_end]

#Chr1    not.give        transposable_element_attribute  4220010 4220012
def read_gff(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ping = '%s.%s' %(unit[0], unit[3])
                data[ping] = 1
    return data


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
   
    ping = read_gff('Compare_unique/HEG4.ALL.ping.gff')
    depth  = seq_depth('RIL.bam.unique.stat')
    size   = lib_size('RIL275.InsertSize.list') 
    #RIL275_RelocaTEi/RelocaTEi_GN1/repeat/results/ALL.all_nonref_insert.gff 
    #RIL275_TEMP/TEMP_GN1/GN1.insertion.refined.bp.summary.gff
    data   = defaultdict(lambda : list())
    output = '%s.summary_clean.table' %(args.input)
    dirs   = glob.glob('%s/RelocaTEi_*' %(args.input))
    for d in dirs:
        ril    = re.split(r'\_', os.path.split(d)[1])[1]
        ril_id = re.sub(r'\D+', '', ril)
        #gff  = '%s/%s.insertion.refined.bp.summary.gff' %(d, ril)
        gff       = '%s/repeat/results/ALL.all_nonref_insert.gff' %(d)
        gff_class = '%s/repeat/results/ALL.all_nonref_insert.characTErized.gff' %(d)
        confident, single_end = parse_gff(gff, ping)
        n_tsd, hom, het, som  = parse_class_gff(gff_class, ping)
        data[ril_id] = [confident, single_end, n_tsd, hom, het, som]

    count = 0
    ofile = open(output, 'w')
    print >> ofile, 'Sample\tInsertSize\tSize_STD\tDepth\tMapped_Depth\tMapped_Rate\tNon_Ref_mPing\tConfident\tEvidence_From_One_End\tCharacterized\tHomozygous\tHeterozygous\tSomatic'
    for strain in sorted(data.keys(), key=int):
        non_ref = data[strain][0] + data[strain][1]
        print >> ofile, 'RIL%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(strain, size[strain][0], size[strain][1], depth[strain][0], depth[strain][1], depth[strain][2], non_ref, data[strain][0], data[strain][1], data[strain][2], data[strain][3], data[strain][4], data[strain][5])
    ofile.close()


if __name__ == '__main__':
    main()

