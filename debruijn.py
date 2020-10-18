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
import itertools
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
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

    return count

def build_seq_dict(seq, kmer_size):
    '''
    kmer_dict but for 1 read only
    '''
    kmer_dict={}
    for kmer in cut_kmer(seq, kmer_size):

        if kmer not in list(kmer_dict.keys()):
            kmer_dict[kmer]=count_in_seq(kmer,seq)

    return kmer_dict

def build_kmer_dict(fastq_file, kmer_size):
    '''
    building a kmer dict
    '''
    kmer_dict={}
    for seq in read_fastq(fastq_file):

        seq_dict=build_seq_dict(seq, kmer_size)

        for kmer in list(seq_dict.keys()):

            if kmer not in list(kmer_dict.keys()):
                kmer_dict[kmer]=seq_dict[kmer]

            else :
                kmer_dict[kmer]=kmer_dict[kmer]+seq_dict[kmer]

    return kmer_dict

def is_linked(prefixe,suffixe):
    '''
    Checks if there should be an edge
    between prefixe and suffixe.
    '''

    if prefixe[1:]==suffixe[:-1]:

        return True

    return False

def build_graph(kmer_dict):
    '''
    Builds a graph from a dictionnary of kmers.
    '''

    my_graph=nx.DiGraph()
    kmers=list(kmer_dict.keys())

    for kmer in kmers:
        my_graph.add_node(kmer[:-1])
        my_graph.add_node(kmer[1:])

    for from_node in my_graph.nodes():
        reminding_nodes=list(my_graph.nodes()).copy()
        reminding_nodes.remove(from_node)

        for to_node in reminding_nodes :

            if is_linked(from_node,to_node):

                if from_node+to_node[-1] in kmers:

                    my_graph.add_edge(from_node,to_node,weight=kmer_dict[from_node+to_node[-1]])

    return my_graph

#==============================================================
#Etape 2

def get_starting_nodes(my_graph):

    '''
    Gets all the nodes that are at the start of paths from a graph.
    '''
    nodes=list(my_graph.nodes())
    edges=list(my_graph.edges())
    edge_targets=[edges[i][1] for i in range(len(edges))]

    starting_nodes=[]

    for node in nodes :

        if node not in edge_targets :
            starting_nodes.append(node)

    return starting_nodes

def get_sink_nodes(my_graph):
    '''
    Gets all the nodes that are at the end of paths from a graph.
    '''
    nodes=list(my_graph.nodes())
    edges=list(my_graph.edges())
    edge_starts=[edges[i][0] for i in range(len(edges))]

    sink_nodes=[]

    for node in nodes :

        if node not in edge_starts :
            sink_nodes.append(node)

    return sink_nodes

def get_contigs(my_graph, starting_nodes, sink_nodes):
    '''
    Gets contigs made of the paths from a graph.
    The contigs start with the starting_nodes and end with the sink_nodes.
    '''
    contig_list=[]

    for start in starting_nodes:

        for sink in sink_nodes:

            for path in nx.all_simple_paths(my_graph, source=start, target=sink):

                contig=path[0]

                for node in path[1:]:
                    contig=contig+node[-1]

                contig_size=len(contig)
                contig_list.append((contig, contig_size))

    return contig_list

def fill(text, width=80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, output_file_name):
    '''
    Save the contigs in a fasta file.
    '''
    output_file=open(output_file_name, "w")

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
    '''
    Get the standard deviation of a list of values
    '''
    return statistics.stdev(values)

def path_average_weight(my_graph, path):
    '''
    Computes the average weight of a given path from a given graph.
    '''
    weights=[]

    for i in range(len(path)-1):

        weight=my_graph[path[i]][path[i+1]]["weight"]
        weights.append(weight)

    average_weight=sum(weights)/len(weights)

    return average_weight

def remove_paths(my_graph, path_list, delete_entry_node, delete_sink_node):
    '''
    Deletes all the paths from path_list in my_graph.
    If delete_entry_node is set to True, the first node of each path is deleted.
    If delete_sink_node is set to True, the last node of each path is deleted.
    '''

    for path in path_list:

        edges_list=[(path[i],path[i+1]) for i in range(len(path)-1)]
        my_graph.remove_edges_from(edges_list)

        nodes_list=[path[i] for i in range(1,len(path)-1)]
        my_graph.remove_nodes_from(nodes_list)

        if delete_entry_node:
            my_graph.remove_node(path[0])

        if delete_sink_node:
            my_graph.remove_node(path[-1])

    return my_graph


def select_best_path(my_graph, path_list, path_length_list, path_average_weight_list,
                     delete_entry_node=False, delete_sink_node=False):

    '''
    Selects the best path from a graph list.
    Priority is given to the highest average weight, and then to the longest path length.
    If several paths are equal, a random path is selected.
    '''
    # comparaison of the average weights :
    weight_std=std(path_average_weight_list)

    if weight_std > 0 :

        max_weight=max(path_average_weight_list)
        removed_path_index=[]

        for i in range(len(path_average_weight_list)):

            if path_average_weight_list[i]<max_weight:

                my_graph=remove_paths(my_graph, [path_list[i]], delete_entry_node, delete_sink_node)
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
            removed_path_index=[]

            for i in range(len(path_list)):

                if path_length_list[i]<max_length:

                    my_graph=remove_paths(my_graph, [path_list[i]],
                                          delete_entry_node,
                                          delete_sink_node)
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

    if len(path_list)>1: # random selection

        random_index=randint(0,len(path_list)-1)

        for i in range(len(path_list)):

            if path_list[i]!=path_list[random_index]:

                my_graph=remove_paths(my_graph, [path_list[i]], delete_entry_node, delete_sink_node)

        path_list=[path_list[random_index]]
        path_length_list=[path_length_list[random_index]]
        path_average_weight_list=[path_average_weight_list[random_index]]

    return my_graph

def solve_bubble(my_graph, ancestor_node, descendant_node):
    '''
    For a given 'bubble' (several paths between two nodes), selects the best path.
    '''
    #finding all the paths between the ancestor and the descendant
    p_list=nx.all_simple_paths(my_graph, source=ancestor_node, target=descendant_node)
    all_paths=[path for path in p_list]

    if len(all_paths)==1:
        return my_graph

    len_list=[len(path) for path in all_paths]
    avg_weight_list=[path_average_weight(my_graph, path) for path in all_paths]

    my_graph=select_best_path(my_graph, all_paths, len_list, avg_weight_list)

    return my_graph


def have_common_ancestor(my_graph, node_1, node_2):
    '''
    Checks if two nodes have a common ancestor.
    '''
    nodes=list(my_graph.nodes())
    nodes.remove(node_1)
    nodes.remove(node_2)

    common_ancestors=[]

    for node in nodes:

        if (nx.has_path(my_graph,
                        source=node,
                        target=node_1)
           )&(
            nx.has_path(my_graph,source=node,target=node_2)):

            common_ancestors.append(node)

    if common_ancestors!=[]:

        total_distance_to_nodes=[]

        for ancestor in common_ancestors:

            distance_to_1=nx.shortest_path_length(my_graph, ancestor, node_1)
            distance_to_2=nx.shortest_path_length(my_graph, ancestor, node_2)
            total_distance_to_nodes.append(distance_to_1+distance_to_2)

        min_distance=min(total_distance_to_nodes)

        i=0
        while total_distance_to_nodes[i] != min_distance:
            i+=1

        closest_ancestor=common_ancestors[i]
        return closest_ancestor

    else :
        return False

def get_edges_to(my_graph, to_node):
    '''
    Gets all the edges that lead to a node.
    '''
    edges_to=[]

    for from_node in my_graph.nodes():

        if (from_node, to_node) in my_graph.edges():
            edges_to.append((from_node, to_node))

    return edges_to

def simplify_bubbles(my_graph):
    '''
    Goes recursively through the entire graph to remove all the bubbles.
    '''
    nodes=list(my_graph.nodes())

    number_of_unchecked=0

    for node in nodes:
        if len(get_edges_to(my_graph, node))>1:
            number_of_unchecked+=1

    if number_of_unchecked>0:

        node_index=0
        node=nodes[node_index]

        while (len(get_edges_to(my_graph, node))<2) & (node_index<len(nodes)-1):
            node_index+=1
            node=nodes[node_index]

        edges=get_edges_to(my_graph, node)
        combinations=list(itertools.combinations(edges, 2))

        comb=combinations[0]
        node_1=comb[0][0]
        node_2=comb[1][0]

        common_ancestor=have_common_ancestor(my_graph, node_1, node_2)
        if common_ancestor!=False:

            my_graph=solve_bubble(my_graph, common_ancestor, node)

            number_of_unchecked-=1

            if number_of_unchecked > 0:

                return simplify_bubbles(my_graph)

            else :

                return my_graph



def have_common_descendant(my_graph, entry_1, entry_2):
    '''
    Checks if two nodes share a descendant.
    '''
    nodes=list(my_graph.nodes())
    nodes.remove(entry_1)
    nodes.remove(entry_2)
    common_descendants=False
    i=0
    common_descendants=[]

    for node in nodes:

        if (nx.has_path(my_graph,entry_1,node))&(nx.has_path(my_graph,entry_2,node)):
            common_descendants.append(node)

    if common_descendants!=[]:

        total_distance_to_entries=[]

        for descendant in common_descendants:

            distance_to_1=nx.shortest_path_length(my_graph, entry_1, descendant)
            distance_to_2=nx.shortest_path_length(my_graph, entry_2, descendant)
            total_distance_to_entries.append(distance_to_1+distance_to_2)

        min_distance=min(total_distance_to_entries)
        i=0

        while total_distance_to_entries[i] != min_distance:
            i+=1

        closest_descendant=common_descendants[i]

        return closest_descendant

    else :

        return False

def solve_entry_tips(my_graph, entry_node_list):
    '''
    Selects the best entry for the graph, based on all the possible entries.
    '''
    combinations=list(itertools.combinations(entry_node_list, 2))


    while len(combinations)>0:

        comb=combinations[0]

        entry_1=comb[0]
        entry_2=comb[1]
        closest_descendant=have_common_descendant(my_graph,
                                                  entry_1,
                                                  entry_2)

        if closest_descendant!=False:

            paths_entry_1_descendant=list(nx.all_simple_paths(my_graph,
                                                              source=entry_1,
                                                              target=closest_descendant))
            paths_entry_2_descendant=list(nx.all_simple_paths(my_graph,
                                                              source=entry_2,
                                                              target=closest_descendant))

            p_list=paths_entry_1_descendant+paths_entry_2_descendant

            all_paths=[path for path in p_list]
            len_list=[len(path) for path in all_paths]
            avg_weight_list=[path_average_weight(my_graph, path)
                             for path in all_paths]

            my_graph=select_best_path(my_graph, all_paths,
                                      len_list,avg_weight_list,delete_entry_node=True)

        combinations=combinations[1:]

    return my_graph

def solve_out_tips(my_graph, out_node_list):

    '''
    Selects the best sink for the graph, based on all the possible sinks.
    '''
    combinations=list(itertools.combinations(out_node_list, 2))

    while len(combinations)>0:

        comb=combinations[0]
        out_1=comb[0]
        out_2=comb[1]

        closest_ancestor=have_common_ancestor(my_graph, out_1, out_2)

        if closest_ancestor!=False:

            paths_out_1_ancestor=list(nx.all_simple_paths(my_graph,
                                                          source=closest_ancestor,
                                                          target=out_1))
            paths_out_2_ancestor=list(nx.all_simple_paths(my_graph,
                                                          source=closest_ancestor,
                                                          target=out_2))

            p_list=paths_out_1_ancestor+paths_out_2_ancestor

            all_paths=[path for path in p_list]
            len_list=[len(path) for path in all_paths]
            avg_weight_list=[path_average_weight(my_graph, path)
                             for path in all_paths]

            my_graph=select_best_path(my_graph,
                                      all_paths,
                                      len_list,
                                      avg_weight_list,
                                      delete_sink_node=True)

        combinations.remove(comb)


    return my_graph

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

    # Read the file and build a graph
    print('building a graph...')
    kmer_dict=build_kmer_dict(fastq,kmer_length)
    my_graph=build_graph(kmer_dict)

    # Resolving the bubbles
    print('resolving the bubbles...')
    my_graph=simplify_bubbles(my_graph)

    #Resolving the starting and sink nodes
    print('resolving the starting and sink nodes...')
    print('starting nodes...')
    starting_nodes=get_starting_nodes(my_graph)
    print('sink nodes...')
    sink_nodes=get_sink_nodes(my_graph)
    print('solving entry...')
    my_graph=solve_entry_tips(my_graph,starting_nodes)
    print('solving out...')
    my_graph=solve_out_tips(my_graph,sink_nodes)

    #Writing the contig(s)
    print('saving the graph to png...')

    labels = nx.get_edge_attributes(my_graph,'weight')
    positions=nx.spring_layout(my_graph, scale=8)
    nx.draw(my_graph, pos=positions, node_size=10)
    nx.draw_networkx_edge_labels(my_graph, pos=positions, edge_labels=labels, font_size=6)
    plt.savefig('Graph.png')

    print('saving the contigs to fasta...')
    starting_nodes=get_starting_nodes(my_graph)
    sink_nodes=get_sink_nodes(my_graph)
    contig_list=get_contigs(my_graph, starting_nodes, sink_nodes)
    save_contigs(contig_list,output)

if __name__ == '__main__':
    main()
