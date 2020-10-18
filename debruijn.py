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
import itertools
import scipy
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
    nx.draw_networkx_edge_labels(G, pos=p, edge_labels=l, font_size=6)
    plt.savefig('Graph.png')
    return G

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
    return contig_list

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

#==============================================================
#Etape 3

def std(values):
    return statistics.stdev(values)

def path_average_weight(G, path):
#     print(path)
    weights=[]
    for i in range(len(path)-1):
        w=G[path[i]][path[i+1]]["weight"]
        weights.append(w)
    average_weight=sum(weights)/len(weights)
    return average_weight

def remove_paths(G, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        edges_list=[(path[i],path[i+1]) for i in range(len(path)-1)]
        G.remove_edges_from(edges_list)
        nodes_list=[path[i] for i in range(1,len(path)-1)]
        G.remove_nodes_from(nodes_list)
    
    if delete_entry_node:
        G.remove_node(path[0])
    
    if delete_sink_node:
        G.remove_node(path[-1])

    return G


def select_best_path(G, path_list, path_length_list, path_average_weight_list, delete_entry_node=False, delete_sink_node=False):
    
    # comparaison of the average weights :
    weight_std=std(path_average_weight_list)

    if weight_std > 0 :

        max_weight=max(path_average_weight_list)
        min_weight=min(path_average_weight_list)
        removed_path_index=[]

        for i in range(len(path_average_weight_list)):

            if path_average_weight_list[i]<max_weight:

                G=remove_paths(G, [path_list[i]], delete_entry_node, delete_sink_node)
                removed_path_index.append(i)

        updated_path_list=path_list.copy()
        updated_length_list=path_length_list.copy()
        updated_average_weight_list=path_average_weight_list.copy()

        for index in removed_path_index:

            updated_path_list.remove(path_list[index])
            updated_length_list.remove(path_length_list[index])
            updated_average_weight_list.remove(path_average_weight_list[index])

        path_list=updated_path_list.copy()
        path_length_list=updated_length_list.copy()
        path_average_weight_list=updated_average_weight_list.copy()

    if len(path_list)>1:

        # comparaison of the lengths :
        length_std=std(path_length_list)

        if length_std > 0:

            max_length=max(path_length_list)
#             min_length=min(path_length_list)
            removed_path_index=[]

            for i in range(len(path_list)):

                if path_length_list[i]<max_length:
#                 if path_length_list[i]>min_length:

                    G=remove_paths(G, [path_list[i]], delete_entry_node, delete_sink_node)
                    removed_path_index.append(i)

            updated_path_list=path_list.copy()
            updated_length_list=path_length_list.copy()
            updated_average_weight_list=path_average_weight_list.copy()

            for index in removed_path_index:

                updated_path_list.remove(path_list[index])
                updated_length_list.remove(path_length_list[index])
                updated_average_weight_list.remove(path_average_weight_list[index])

            path_list=updated_path_list.copy()
            path_length_list=updated_length_list.copy()
            path_average_weight_list=updated_average_weight_list.copy()
            
    if len(path_list)>1: #random
        
        random_index=random.randint(0,len(path_list)-1)

        for i in range(len(path_list)):
                
                if path_list[i]!=path_list[random_index]:
                    
                    G=remove_paths(G, [path_list[i]], delete_entry_node, delete_sink_node)
        
        path_list=[path_list[random_index]]
        path_length_list=[path_length_list[random_index]]
        path_average_weight_list=[path_average_weight_list[random_index]]

    return G

def solve_bubble(G, ancestor_node, descendant_node):
    p_list=nx.all_simple_paths(G, source=ancestor_node, target=descendant_node) #finding all the paths between the ancestor and the descendant
    all_paths=[path for path in p_list]
    len_list=[len(path) for path in all_paths]
    avg_weight_list=[path_average_weight(G, path) for path in all_paths]
    G=select_best_path(G, all_paths, len_list, avg_weight_list)
    
    return G

def have_common_ancestor(G, node_1, node_2):
    common_ancestor=False
    nodes=list(G.nodes())
    nodes.remove(node_1)
    nodes.remove(node_2)
    i=0
    common_ancestors=[]
    
    for i in range(len(nodes)):
        node=nodes[i]
        node_1_path_list=[path for path in nx.all_simple_paths(G, source=node, target=node_1)]
        node_2_path_list=[path for path in nx.all_simple_paths(G, source=node, target=node_2)]
        if (node_1_path_list !=[])&(node_2_path_list!=[]):
            common_ancestors.append(node)
    if common_ancestors!=[]:
        total_distance_to_nodes=[]
        for ancestor in common_ancestors:
            distance_to_1=nx.shortest_path_length(G, ancestor, node_1)
            distance_to_2=nx.shortest_path_length(G, ancestor, node_2)
            total_distance_to_nodes.append(distance_to_1+distance_to_2)
        
        min_distance=min(total_distance_to_nodes)
        i=0
        while total_distance_to_nodes[i] != min_distance:
            i+=1
        closest_ancestor=common_ancestors[i]
        return closest_ancestor
    else :
        return False

def get_edges_to(G, to_node):
    edges_to=[]
    for from_node in G.nodes():
        if (from_node, to_node) in G.edges():
            edges_to.append((from_node, to_node))
    return edges_to

def simplify_bubbles(G):
    nodes=list(G.nodes())
    bubble=False
    i=0
    for node in nodes:
#         print('node',node)
        in_edges=get_edges_to(G,node)
        
        if len(in_edges)>=2:
            combinations=list(itertools.combinations(in_edges, 2))[0]
            i=0
#             print('combi',combinations)
            c=combinations[i]
#             print('c',c)
            node_1=c[0]
            node_2=c[1]
#             print(nodes)
#             print(node_1)
#             print(node_2)
            while (have_common_ancestor(G, node_1, node_2)==False)&(i<len(combinations)-1):
                i+=1
                c=combinations[i]
                node_1=c[0]
                node_2=c[1]
#             print(node_1,node_2)
            common_ancestor=have_common_ancestor(G, node_1, node_2)
            if common_ancestor!=False:
                G=solve_bubble(G, common_ancestor, node)
#                 print(list(G.nodes()))
                G=simplify_bubbles(G)

    return G


def have_common_descendant(G, entry_1, entry_2):
    common_descendants=False
    nodes=list(G.nodes())
    nodes.remove(entry_1)
    nodes.remove(entry_2)
    i=0
    common_descendants=[]
    
    for i in range(len(nodes)):
        node=nodes[i]
        if (nx.has_path(G,entry_1,node))&(nx.has_path(G,entry_2,node)):
            common_descendants.append(node)
    if common_descendants!=[]:
        total_distance_to_entries=[]
        for descendant in common_descendants:
            distance_to_1=nx.shortest_path_length(G, entry_1, descendant)
            distance_to_2=nx.shortest_path_length(G, entry_2, descendant)
            total_distance_to_entries.append(distance_to_1+distance_to_2)
        
        min_distance=min(total_distance_to_entries)
        i=0
        while total_distance_to_entries[i] != min_distance:
            i+=1
        closest_descendant=common_descendants[i]
        return closest_descendant
    else :
        return False

def solve_entry_tips(G, entry_node_list):

    combinations=list(itertools.combinations(entry_node_list, 2))
    while len(combinations)>0:
        c=combinations[0]
#         print('c',c)
        entry_1=c[0]
        entry_2=c[1]
        closest_descendant=have_common_descendant(G, entry_1, entry_2)

        paths_entry_1_descendant=list(nx.all_simple_paths(G, source=entry_1, target=closest_descendant))
        paths_entry_2_descendant=list(nx.all_simple_paths(G, source=entry_2, target=closest_descendant))

        p_list=paths_entry_1_descendant+paths_entry_2_descendant

        all_paths=[path for path in p_list]
        len_list=[len(path) for path in all_paths]
        avg_weight_list=[path_average_weight(G, path) for path in all_paths]

        G=select_best_path(G, all_paths, len_list, avg_weight_list,delete_entry_node=True)
        combinations.remove(c)
        
    return G


def solve_out_tips(G, out_node_list):
    combinations=list(itertools.combinations(out_node_list, 2))
    while len(combinations)>0:
        c=combinations[0]
#         print('c',c)
        out_1=c[0]
        out_2=c[1]
        closest_ancestor=have_common_ancestor(G, out_1, out_2)

        paths_out_1_ancestor=list(nx.all_simple_paths(G, source=closest_ancestor, target=out_1))
        paths_out_2_ancestor=list(nx.all_simple_paths(G, source=closest_ancestor, target=out_2))

        p_list=paths_out_1_ancestor+paths_out_2_ancestor

        all_paths=[path for path in p_list]
        len_list=[len(path) for path in all_paths]
        avg_weight_list=[path_average_weight(G, path) for path in all_paths]

        G=select_best_path(G, all_paths, len_list, avg_weight_list,delete_sink_node=True)
        combinations.remove(c)
    print(G.nodes())
    return G


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

    #Part 1
    print('Part 1')
    Sequences = read_fastq(fastq)
    kmer_dict=build_kmer_dict(fastq,kmer_length)
    G=build_graph(kmer_dict)
    
    #Part 2
    print('Part 2')
    starting_nodes=get_starting_nodes(G)
    sink_nodes=get_sink_nodes(G)
    contig_list=get_contigs(G, starting_nodes, sink_nodes)
    
    save_contigs(contig_list,output)
    
    #Part 3
    print('Part 3')
    

    
if __name__ == '__main__':
    main()

# main()