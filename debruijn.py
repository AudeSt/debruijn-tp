# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
import random
from random import randint
import statistics
import argparse
import os
import sys
from operator import itemgetter
import networkx as nx
import matplotlib
random.seed(9001)

__author__ = "Aude STERLIN"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Aude STERLIN"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aude STERLIN"
__email__ = "sterlin.aude@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
        :Parameters:
            path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
       Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

#==============================================================
# Etape 1
fastq=get_arguments().fastq_file

def read_fastq(fastq_file):
    '''
    Sequences generator.
    Testing :

    result = read_fastq(fastq_file)
    print('first line :',next(result))
    print('second line :',next(result))

    '''
    with open(fastq_file) as file :
        for line in file.readlines():
            if line[0] in ['A','T','C','G']:
                yield line

Sequences = read_fastq(fastq)

seq0=next(Sequences)

def cut_kmer(seq, kmer_size):
    '''
    Cutting all the kmers from a sequence
    '''
    last_index=-1
    first_base_index=last_index+1
    last_base_index=first_base_index+kmer_size
#     print(kmer_size, first_base_index,last_base_index)
    while last_base_index<=len(seq):
        yield seq[first_base_index:last_base_index]
        last_index=last_base_index-1
        first_base_index=last_index+1
        last_base_index=first_base_index+kmer_size

    yield seq[first_base_index:].replace('\n','')

# cut=cut_kmer('AAuAAt',3)
# print('1',next(cut))
# print('2',next(cut))
# print('3',next(cut))

def build_kmer_dict(fastq_file, kmer_size):
    '''
    building a kmer dict
    '''
    kmer_dict_ref={}
    all_seq=read_fastq(fastq_file)
    for seq in read_fastq(fastq_file):
        kmer_dict={}
        for kmer in cut_kmer(seq, kmer_size):
            if kmer not in list(kmer_dict.keys())+[''] :
                kmer_dict[kmer]=seq.count(kmer)

        for kmer in list(kmer_dict.keys()):
            if kmer not in list(kmer_dict_ref.keys()):
                kmer_dict_ref[kmer]=kmer_dict[kmer]
            else :
                kmer_dict_ref[kmer]=kmer_dict_ref[kmer]+kmer_dict[kmer]

    return kmer_dict_ref

        
print(build_kmer_dict(fastq,5))
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    main()

