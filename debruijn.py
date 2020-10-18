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
import matplotlib
import matplotlib.pyplot as plt
import pickle
import networkx as nx
import hashlib
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
# args=get_arguments()
# fastq=args.fastq_file
# kmer_length=args.kmer_size

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

# Sequences = read_fastq(fastq)

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
    '''
    Counts the number of kmers in the sequence
    '''
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

# kmer_dict=build_kmer_dict(fastq,kmer_length)

def is_linked(prefixe,suffixe):
    if prefixe[1:]==suffixe[:-1]:
        return True

def build_graph(kmer_dict):
    G=nx.DiGraph()
    kmers=list(kmer_dict.keys())
#     print(kmers)
    for kmer in kmers:
        G.add_node(kmer[:-1])
        G.add_node(kmer[1:])

    for from_node in G.nodes():
        reminding_nodes=list(G.nodes()).copy()
        reminding_nodes.remove(from_node)

        for to_node in reminding_nodes :
#             print(from_node, to_node)
#             print(is_linked(from_node,to_node))
            if is_linked(from_node,to_node):
                if from_node+to_node[-1] in kmers:
                    G.add_edge(from_node,to_node,weight=kmer_dict[from_node+to_node[-1]])


    l = nx.get_edge_attributes(G,'weight')
#     print(l)
#     p=nx.spring_layout(G, scale=8)
#     nx.draw(G, pos=p, node_size=2)
#     nx.draw_networkx_edge_labels(G, pos=p, edge_labels=l, font_size=6)
    p=nx.spring_layout(G, scale=8)
    nx.draw(G, pos=p, node_size=10)

    plt.savefig('Graph.png')
    return G

# G1=build_graph(kmer_dict)

# file = open(os.path.abspath(os.path.join(os.path.dirname(__file__), "kmer.pck")),'rb')
# kmer_dict = pickle.load(file)
# graph = build_graph(kmer_dict)
# #TCAGAGA
# #TCA  TC CA
# #CAG CA AG
# #AGA AG GA
# #GAG GA AG
# #AGA AG GA
# print( graph.number_of_nodes() == 4)
# print( graph.number_of_edges() == 4)
# print( "AG" in graph)
# print( "GA" in graph)
# print( graph.edges["AG", "GA"]['weight'] == 2)
# file.close()
#==============================================================
#Etape 2

def get_starting_nodes(G):
    nodes=list(G.nodes())
    edges=list(G.edges())
    edge_targets=[edges[i][1] for i in range(len(edges))]
#     print(edges, edge_targets)
    
    starting_nodes=[]
    for node in nodes :
#         print(node)
        if node not in edge_targets :
            starting_nodes.append(node)
    return(starting_nodes)

# starting_nodes=get_starting_nodes(G1)
# print(starting_nodes)

# graph = nx.DiGraph()
# graph.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
# nodes = get_starting_nodes(graph)    
# print( len(nodes) == 2)
# print( 1 in nodes)
# print( 3 in nodes)

def get_sink_nodes(G):
    nodes=list(G.nodes())
    edges=list(G.edges())
    edge_starts=[edges[i][0] for i in range(len(edges))]
#     print(edges, edge_starts)
    
    sink_nodes=[]
    for node in nodes :
#         print(node)
        if node not in edge_starts :
            sink_nodes.append(node)
    return(sink_nodes)
# sink_nodes=get_sink_nodes(G1)
# print(sink_nodes)

# graph = nx.DiGraph()
# graph.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
# nodes = get_sink_nodes(graph)
# print( len(nodes) == 2)
# print( 6 in nodes)
# print( 7 in nodes)


def get_contigs(G, starting_nodes, sink_nodes):
    nodes=list(G.nodes())
    edges=list(G.edges())
    contig_list=[]
    for start in starting_nodes:
        for sink in sink_nodes:
            for path in nx.all_simple_paths(G, source=start, target=sink):
                contig=path[0]
                for node in path[1:]:
                    contig=contig+node[-1]
                contig_size=len(contig)
                contig_list.append((contig, contig_size))
#     contigs=[(contig, taille du contig)*plein]
    return contig_list

# contig_list=get_contigs(G1, starting_nodes, sink_nodes)

# print(contig_list)
# graph = nx.DiGraph()
# graph.add_edges_from([("TC", "CA"), ("AC", "CA"), ("CA", "AG"), ("AG", "GC"), ("GC", "CG"), ("CG", "GA"), ("GA", "AT"), ("GA", "AA")])
# contig_list = get_contigs(graph, ["TC", "AC"], ["AT" , "AA"])
# results = ["TCAGCGAT", "TCAGCGAA", "ACAGCGAT", "ACAGCGAA"]
# print( len(contig_list) == 4)
# for contig in contig_list:
#     print( contig[0] in results)
#     print( contig[1] == 8)

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, output_file_name):
    output_file=open(output_file_name, "w")
#     print(output_file)
    for i in range(len(contig_list)):
        contig_tuple=contig_list[i]
        contig=contig_tuple[0]
        contig_size=contig_tuple[1]
        output_file.write('>contig_{0} len={1}'.format(i,contig_size)+'\n')
        output_file.write(fill(contig)+'\n')

    output_file.close()
    return 0


# save_contigs(contig_list,' eva71.fasta')

#==============================================================
# Main program
#==============================================================

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    fastq=args.fastq_file
    kmer_length=args.kmer_size
    output=os.curdir + os.sep+args.output_file
    print(output)
    #Part 1
    Sequences = read_fastq(fastq)
    kmer_dict=build_kmer_dict(fastq,kmer_length)
    G=build_graph(kmer_dict)
    
    #Part 2
    starting_nodes=get_starting_nodes(G)
    sink_nodes=get_sink_nodes(G)
    contig_list=get_contigs(G, starting_nodes, sink_nodes)
    
    save_contigs(contig_list,output)

    
if __name__ == '__main__':
    main()

# main()