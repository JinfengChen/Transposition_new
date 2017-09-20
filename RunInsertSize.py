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
python RunInsertSize.py --input /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/

Get insert size for bam in directory of input.
--input: directory of bam files, GN160.bam
--project: output directory

    '''
    print message


def runjob(script, lines):
    #cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource walltime=100:00:00,mem=5G --convert no %s' %(lines, script)
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 5 --lines %s --interval 120 --task 1 --mem 20G --time 100:00:00 --convert no %s' %(lines, script)
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
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'Rice_pop'
    
    project = os.path.abspath(args.project)
    shell = open ('%s.InsertSize.sh' %(project), 'w')
    #print project
    bams  = glob.glob('%s/*.bam' %(args.input))
    count = 0
    for bam in sorted(bams):
        #print bam
        prefix = os.path.split(os.path.splitext(bam)[0])[1]
        cfg    = '%s/%s.config' %(args.input, prefix)
        if not os.path.exists(cfg) or os.path.getsize(cfg) == 0:
            #print bam
            bam2cfg  ='perl /rhome/cjinfeng/BigData/software/SVcaller/breakdancer-1.1.2/bin/bam2cfg.pl %s > %s' %(bam, cfg)
            print >> shell, bam2cfg
            count += 1
            #print 'done'
    shell.close()

    if count > 0:
        #print 'Run Job'
        runjob('%s.InsertSize.sh' %(project), 1)
    else:
        print 'No Job'
    #os.system('rm -R *.InsertSize.sh*')
 
    ofile = open ('%s.InsertSize.list' %(project), 'w')
    for bam in sorted(bams):
        prefix = os.path.split(os.path.splitext(bam)[0])[1]
        cfg    = '%s/%s.config' %(args.input, prefix)
        size = getsize(cfg)
        #size, std, lower, upper
        print >> ofile, '%s\t%s\t%s\t%s\t%s' %(bam, size[0], size[1], size[2], size[3])
    ofile.close()

if __name__ == '__main__':
    main()

