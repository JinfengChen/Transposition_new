#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import time
from Bio import SeqIO

def usage():
    test="name"
    message='''
python RunTEMP_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam

Run TEMP on population level data.
--input: directory of bam files, GN160.bam
--genome: reference genome sequence, need to index using bwa first
--repeat: fasta repeat sequence, do not have too similar element like mPing and ping, which are >99% identity
--project: output directory

    '''
    print message

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource walltime=100:00:00,mem=5G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#readgroup:RIL1_0_CGTACG_FC153L5	platform:Illumina	map:../input/bam/GN1.bam	readlen:100.18	lib:RIL1_0_CGTACG_FC153L5	num:10001	lower:76.38	upper:303.56	mean:141.39	std:25.74	SWnormality:-71.20	exe:samtools view
def getsize(infile):
    size = [] 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                size = [unit[8], unit[9], unit[6], unit[7]]
                for i in range(len(size)):
                    size[i] = re.sub(r'.*:', '', size[i]) 
    return size

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-e', '--existingRun')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa'
  
    if not args.repeat:
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping_pogo.fa'
    
    if not args.project:
        args.project = 'Rice_pop'
    
    TEMP        = '/rhome/cjinfeng/software/tools/TEMP'
    Reference   = os.path.abspath(args.genome)
    Repeat      = os.path.abspath(args.repeat)
    existingTE  = '%s.RepeatMasker.out.bed' %(Reference)       
 
    project = os.path.abspath(args.project)
    if not os.path.exists(project):
        os.mkdir(project)
    #ofile = open ('%s.TEMP.run.sh' %(project), 'w')
    #print project
    bams  = glob.glob('%s/*.bam' %(args.input))
    count = 0
    for bam in sorted(bams):
        #print bam
        cmd    = []
        count += 1
        prefix = os.path.split(os.path.splitext(bam)[0])[1]
        tempdir= '%s/%s' %(project, 'TEMP_%s' %(prefix))
        if not os.path.exists(tempdir):
            os.mkdir(tempdir)
        tempout = os.path.abspath('%s.insertion.refined.bp.summary' %(prefix))
        if not os.path.isfile(tempout) and not os.path.isfile('%s/%s' %(tempdir, tempout)):
            cfg      = '%s/%s.config' %(args.input, prefix)
            size     = []
            if os.path.isfile(cfg):
                size = getsize(cfg)
            lib_size = size[0] if len(size) == 4 else 500
            #print '%s\t%s\t%s\t%s' %(bam, lib_size, lib_size, lib_sd)
            #bash $TEMP/scripts/TEMP_Insertion.sh -i $bam -s $TEMP/scripts -r $TE_Fasta -t $TE_BED -m 3 -f 400 -c 8
            tempout = os.path.abspath('%s.insertion.refined.bp.summary' %(prefix))
            cmd.append('bash %s/scripts/TEMP_Insertion.sh -i %s -s %s/scripts -r %s -t %s -m 3 -f %s -c 1 > %s.log 2> %s.std' %(TEMP, bam, TEMP, Repeat, existingTE, lib_size, prefix, prefix))
            cmd.append('python /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/TEMP/TEMP2GFF.py --insertion %s' %(tempout))
            cmd.append('rm %s.bam %s.bam.bai' %(os.path.abspath(prefix), os.path.abspath(prefix)))
            cmd.append('mv %s.* %s' %(os.path.abspath(prefix), os.path.abspath(tempdir)))
            cmd.append('mv %s %s' %(os.path.abspath(tempdir), project))
            #print cmd 
            for c in cmd:
                #print >> ofile, c
                os.system(c)
    #ofile.close()
    #runjob('%s.TEMP.run.sh' %(project), 5)
 
if __name__ == '__main__':
    main()

