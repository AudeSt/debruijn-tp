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
import matplotlib.pyplot as plt
import pickle
matplotlib.use('Agg')
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
args=get_arguments()
fastq=args.fastq_file
kmer_length=args.kmer_size

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
                yield line[:-1]

Sequences = read_fastq(fastq)

def cut_kmer(seq, kmer_size):
    '''
    Cutting all the kmers from a sequence
    '''
    index=0
#     print(kmer_size, first_base_index,last_base_index)
    while index+kmer_size<=len(seq):
        yield seq[index:index+kmer_size]
        index+=1
        
def count_in_seq(kmer,seq):
    count=0
    for i in range(len(seq)-len(kmer)+1):
        seq_slice=seq[i:i+len(kmer)]
        if seq_slice==kmer:
            count+=1
#     print('c',count)
    return count

def build_seq_dict(seq, kmer_size):
    '''
    kmer_dict but for 1 read only
    '''
    kmer_dict={}
    for kmer in cut_kmer(seq, kmer_size):
#         print(seq.count(kmer))
        if kmer not in list(kmer_dict.keys()):
            kmer_dict[kmer]=count_in_seq(kmer,seq)
            
    return kmer_dict

def build_kmer_dict(fastq_file, kmer_size):
    '''
    building a kmer dict
    '''
    kmer_dict={}
    for seq in read_fastq(fastq_file):
#         print('seq',seq)
        seq_dict=build_seq_dict(seq, kmer_size)
#         print('list of seq keys :',list(seq_dict.keys()) )
        for kmer in list(seq_dict.keys()):
#             print('kmer in keys :', kmer)
            
            if kmer not in list(kmer_dict.keys()):
                kmer_dict[kmer]=seq_dict[kmer]
            else :
                kmer_dict[kmer]=kmer_dict[kmer]+seq_dict[kmer]

    return kmer_dict

kmer_dict=build_kmer_dict(fastq,kmer_length)
# print(kmer_dict)
# print(len(kmer_dict.keys()) == 4)
# print ("TCA" in kmer_dict)
# print ("CAG" in kmer_dict)
# print ("AGA" in kmer_dict)
# print ("GAG" in kmer_dict)
# print (kmer_dict["AGA"] == 2)

# def build_predecesseurs(kmer, node_list):
#     precedents = [i+kmer[:-1] for i in ['A','T','C','G']]
#     predecesseurs=[]
#     for from_kmer in node_list :
#         if from_kmer in precedents :
#             predecesseurs.append(from_kmer)
#     return predecesseurs

def build_graph(kmer_dict):
    G=nx.DiGraph()
    kmers=list(kmer_dict.keys())
    print('kmers',kmers)
    for kmer in kmers:
        G.add_node(kmer)
        kmers_rm=kmers.remove(kmer)
        for to_kmer in kmers:
            if to_kmer[:-1]==kmer[1:]:
                G.add_edge(kmer, to_kmer, weight=kmer_dict[kmer])

    l = nx.get_edge_attributes(G,'weight')
#     print(l)
    p=nx.spring_layout(G, scale=8)
    nx.draw(G, pos=p, node_size=2)
    nx.draw_networkx_edge_labels(G, pos=p, edge_labels=l, font_size=6)


    plt.savefig('Graph.png')
    return G

# G1=build_graph(kmer_dict)
file = open(os.path.abspath(os.path.join(os.path.dirname(__file__), "kmer.pck")),'rb')
kmer_dict = pickle.load(file)
print(kmer_dict)
graph = build_graph(kmer_dict)
print(graph.number_of_nodes() == 4)
print(graph.number_of_edges() == 4)
print(graph.nodes())
print("AGA" in graph)
print("GAG" in graph)
print(graph.edges["AGA", "GAG"]['weight'] == 2)
file.close()
#==============================================================
#Etape 2



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
