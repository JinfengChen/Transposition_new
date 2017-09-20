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
#RIL1    6.111434656     200     184     16      0
#RIL2    5.109878359     128     126     2       0
def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if unit[0].startswith(r'RIL'):
                    ril = re.sub(r'RIL', r'', unit[0])
                    data[ril] = '\t'.join(unit[1:])
                elif unit[0].startswith(r'Sample'):
                    data[unit[0]] = '\t'.join(unit[1:])
    return data




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--table1')
    parser.add_argument('--table2')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.table1) > 0
    except:
        usage()
        sys.exit(2)

    t1 = readtable(args.table1)
    t2 = readtable(args.table2)
    if t1.has_key('Sample'):
        print 'Sample\t%s' %(t1['Sample'])
        del t1['Sample']
    if t2.has_key('Sample'):
        print 'Sample\t%s' %(t2['Sample'])
        del t2['Sample']
    for ril in sorted(t1.keys(), key=int):
        len2 = len(t2[t2.keys()[0]])
        new_inf = t2[ril] if t2.has_key(ril) else '\t'.join(['NA']*int(len2))
        print 'RIL%s\t%s\t%s' %(ril, t1[ril], new_inf)

if __name__ == '__main__':
    main()

