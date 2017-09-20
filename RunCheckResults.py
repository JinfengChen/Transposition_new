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
python RunCheckResults.py --input RIL275_RelocaTEi --tools RelocaTEi
python RunCheckResults.py --input RIL275_TEMP --tools TEMP
python RunCheckResults.py --input RIL275_RelocaTE --tools RelocaTE

Check if the final result file is empty or have too few insertions
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def test_file(infile):
    data = defaultdict(str)
    if os.path.getsize(infile) == 0:
        return 0
    else:
        line_n = 0
        with open (infile, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2:
                    line_n += 1 
        return line_n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-t', '--tools')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    #RIL275_RelocaTEi/RelocaTEi_GN1/repeat/results/ALL.all_nonref_insert.gff 
    #RIL275_TEMP/TEMP_GN1/GN1.insertion.refined.bp.summary.gff
    dirs = glob.glob('%s/%s_*' %(args.input, args.tools))
    for d in dirs:
        ril = re.split(r'\_', os.path.split(d)[1])[1]
        if args.tools == 'TEMP':
            gff = '%s/%s.insertion.refined.bp.summary.gff' %(d, ril)
            n   = test_file(gff)
            if int(n) <= 10:
                print ril, gff, n
        elif args.tools == 'RelocaTEi':
            gff = '%s/repeat/results/ALL.all_nonref_insert.gff' %(d)
            n   = test_file(gff)
            if int(n) <= 10:
                print ril, gff, n
        elif args.tools == 'RelocaTE':
            gff = '%s/mping/results/%s.mping.all_inserts.gff' %(d, ril)
            n   = test_file(gff)
            if int(n) <= 10:
                print ril, gff, n

if __name__ == '__main__':
    main()

