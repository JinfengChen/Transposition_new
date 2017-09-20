#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Sample  Coverage        non-ref Homo    Het     Som
#RIL39   93.62213026     183     183     0       0
#RIL106  8.246933301     228     227     1       0
def readtable(infile):
    data   = defaultdict(str)
    output = '%s.sorted.table' %(os.path.splitext(infile)[0])
    #print output1
    #output = 'sorted.table'
    ofile = open(output, 'w')
    with open (infile, 'r') as filehd:
        header = filehd.readline().rstrip()
        print >> ofile, header
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit   = re.split(r'\t',line)
                ril_id = int(re.sub(r'RIL', '', unit[0]))
                data[ril_id] = line
    
    for i in sorted(data.keys(), key=int):
        print >> ofile, data[i]
    ofile.close()
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    readtable(args.input)

if __name__ == '__main__':
    main()

