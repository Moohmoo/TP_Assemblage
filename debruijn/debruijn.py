#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
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

import argparse
from asyncore import write
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "OUSSAREN Mohamed"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["OUSSAREN Mohamed"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "OUSSAREN Mohamed"
__email__ = "oussarenmohamed@gmail.com"
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
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """ Read a fastq file and extract a generator sequence.

        Parameters
        ----------
        fastq : string
            contains file name of a fasq file.

        Return
        ------
        yield object
            return a generator object.
    """
    bases = ["A", "T", "C", "G"]
    with open(fastq_file, "r") as fileas:
        for line in fileas:
            if line.startswith(tuple(bases)) :
                yield line.strip()


def cut_kmer(read, kmer_size):
    """ Cut a sequence in k mers and extract a generator of k-mers.

        Paramaters
        ----------
        read : string 
            the read sequence.
        kmer_size : int
            The length of the kmers.

        Return
        ------
        yield object
            return a generator object.
    """
    i = 0
    while True:
        if (i + kmer_size > len(read) or i == len(read) - 1):
            break
        yield read[i:i + kmer_size]
        i = i + 1


def build_kmer_dict(fastq_file, kmer_size):
    """ Generate a dictionnary contains the number of occurrences per kmers.

        Paramaters
        ----------
        fastq_file : string 
            the name of fastq file.
        kmer_size : int
            The length of the kmers.

        Return
        ------
        dictionnary
            contains the number of occurrences per kmers.
    """
    kmers_dict = {}
    fastq = read_fastq(fastq_file)

    for sequence in fastq :
        set_kmers = cut_kmer(sequence, kmer_size)
        for kmers in set_kmers :
            if kmers not in kmers_dict:
                kmers_dict[kmers] = 0
            kmers_dict[kmers] += 1

    return kmers_dict


def build_graph(kmer_dict):
    """ Build a graph by using networkx package.

        Paramaters
        ----------
        kmer_dict : dictionnary
            contains the number of occurrences per kmers.

        Return
        ------
        object networkx
            contains the graph.
    """
    graph = nx.DiGraph()
    kmers_keys = list(kmer_dict.keys())
    i = 0
    while i + 1 <= len(kmers_keys) :
        l = kmers_keys[i][:-1]
        r = kmers_keys[i][1:]
        graph.add_edge(l, r, weight=kmer_dict[kmers_keys[i]])
        i = i + 1

    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ Remove a path from a given graph.

        Paramaters
        ----------
        graph : object networkx
            contains the graph.
        path_list : list of path
            contains different path to delete.
        delete_entry_node : boolean
            flag to know whether to delete entry node.
        delete_sink_node : boolean
            flag to know whether to delete sink node.

        Return
        ------
        graph object
            contains a graph object cleaned
    """

    index_entry = 0 if delete_entry_node else 1
    index_sink = None if delete_sink_node else -1

    for path in path_list:
        graph.remove_nodes_from(path[index_entry:index_sink])

    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """ Remove a path from a given graph.

        Paramaters
        ----------
        graph : object networkx
            contains the graph.
        path_list : list 
            contains different path to delete.
        path_length : list
            contains length of each path
        weight_avg_list : list
            contains weigth average of each path
        delete_entry_node : boolean
            flag to know whether to delete entry node.
        delete_sink_node : boolean
            flag to know whether to delete sink node.

        Return
        ------
        graph object
            contains a graph object cleaned
    """
    max_weigth = max(weight_avg_list)
    index_weigth = []

    for i in range(len(weight_avg_list)):
        if weight_avg_list[i] == max_weigth:
            index_weigth.append(i)
    
    if len(index_weigth) >= 2:
        max_path_length = max(path_length)
        shift = []

        for i in index_weigth:
            if path_length[i] == max_path_length:
                shift.append(i)

        best_i = random.choice(shift)
    else:
        best_i = index_weigth[0]

    temp = path_list
    temp.pop(best_i)

    graph_cleared = remove_paths(graph, temp, delete_entry_node, delete_sink_node)
    return graph_cleared


def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """ Solve a bubble in networkx graph.

    Paramaters
        ----------
        graph : object networkx
            contains the graph.
        ancestor_node : object
            contains informations about ancestor node
        descandant_node : object
            contains informations about descendant node

        Return
        ------
        graph object
            contains a graph object cleaned
    """
    weights = []
    lengths = []

    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))

    for path in paths:
        print(path)
        lengths.append(len(path))
        weights.append(path_average_weight(graph, path))

    graph = select_best_path(graph, paths, lengths, weights)
    return graph

def simplify_bubbles(graph):
    """ Clean all bubble in a networkx network.
    
        Paramaters
        ----------
        graph : object networkx
            contains the graph.

        Return
        ------
        graph object
            contains a graph object cleaned
    """
    buble_to_clear = False
    nodes = list(graph)

    for node in nodes:
        node_origin = node

        if graph.in_degree(node) > 1:
            predecessors = list(graph.predecessors(node))

            if len(predecessors) <= 2:
                buble_to_clear = True
                break

        if buble_to_clear:
            break

    if not buble_to_clear:
        ancestor = nx.lowest_common_ancestor(graph, predecessors[0], predecessors[1])
        graph = solve_bubble(graph, ancestor, node_origin)

    return graph

def solve_entry_tips(graph, starting_nodes):
    """ Clean all bubble in a networkx graph depending starting nodes.
    
        Paramaters
        ----------
        graph : object networkx
            contains the graph.
        starting_nodes : list 
            contains all starting nodes.

        Return
        ------
        graph object
            contains a graph object cleaned.
    """
    nodes = []
    ancestors = []

    if len(starting_nodes) >= 3:
        for i in starting_nodes:
            for j in starting_nodes:
                nodes.append((i, j))
    else:
        nodes = [tuple(starting_nodes)]

    for node in nodes:
        ancestors.append(nx.lowest_common_ancestor(graph.reverse(),
                                             node[0], node[1]))


    weights = []
    paths = []
    lengths = []

    for i in range(len(nodes)):
        path1 = list(nx.all_simple_paths(graph, nodes[i][0], ancestors[i]))[0]
        path2 = list(nx.all_simple_paths(graph, nodes[i][1], ancestors[i]))[0]

        weights.append(path_average_weight(graph, path1))
        weights.append(path_average_weight(graph, path2))

        paths.append(path1)
        paths.append(path2)

        lengths.append(len(path1))
        lengths.append(len(path2))

    graph = select_best_path(graph, paths, lengths, weights, True)

    return graph

def solve_out_tips(graph, ending_nodes):
    """ Clean all bubble in a networkx graph depending ending nodes.
    
        Paramaters
        ----------
        graph : object networkx
            contains the graph.
        ending_nodes : list 
            contains all ending nodes.

        Return
        ------
        graph object
            contains a graph object cleaned.
    """

    
    nodes = []
    ancestors = []

    if len(ending_nodes) >= 3:
        for i in ending_nodes:
            for j in ending_nodes:
                nodes.append((i, j))
    else:
        nodes = [tuple(ending_nodes)]

    for node in nodes:
        ancestors.append(nx.lowest_common_ancestor(graph,
                                             node[0], node[1]))


    weights = []
    paths = []
    lengths = []

    for i in range(len(ancestors)):
        path1 = list(nx.all_simple_paths(graph, ancestors[i], nodes[i][0]))[0]
        path2 = list(nx.all_simple_paths(graph, ancestors[i], nodes[i][1]))[0]

        weights.append(path_average_weight(graph, path1))
        weights.append(path_average_weight(graph, path2))

        paths.append(path1)
        paths.append(path2)

        lengths.append(len(path1))
        lengths.append(len(path2))

    graph = select_best_path(graph, paths, lengths, weights, delete_sink_node=True)

    return graph
    

def get_starting_nodes(graph):
    """ Return a list of input nodes.

        Paramaters
        ----------
        graph : object networkx
            contains the graph.

        Return
        ------
        list
            contains all the input nodes.
    """
    node_list = [node for node in graph.nodes() if len(tuple(graph.predecessors(node))) == 0]
    return node_list

def get_sink_nodes(graph):
    """ Return a list of output nodes.

        Paramaters
        ----------
        graph : object networkx
            contains the graph.

        Return
        ------
        list
            contains all the output nodes.
    """
    node_list = [node for node in graph.nodes() if len(tuple(graph.successors(node))) == 0]
    return node_list


def get_contigs(graph, starting_nodes, ending_nodes):
    """ Return all contigs depending starting nodes and ending nodes.

        Paramaters
        ----------
        graph : object networkx
            contains the graph.
        starting_nodes : list
            contains a list of nodes corresponding to starting nodes
        ending_nodes : list
            contains a list of nodes corresponding to ending nodes

        Return
        ------
        list
            contains all possible contigs with their length.
    """
    contigs = []
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:

            if nx.has_path(graph, starting_node, ending_node) :
                path_contig = nx.all_simple_paths(graph, starting_node, ending_node)
                path_contig = list(path_contig)[0]

                for path_node in path_contig :
                    contig = ''.join(path_contig[::len(path_node)])
                contigs.append((contig, len(contig)))

    return contigs

def save_contigs(contigs_list, output_file):
    """ Writes an output file containing the contigs according to the fasta format.

        Paramaters
        ----------
        contigs_list : list
            contains all possible contigs with their length.
        output_file : the name of the output file

        Return
        ------
        None
    """
    with open(output_file, "w") as fileout:
        i = 0
        for contig in contigs_list:
            seq_contig = contig[0]
            length_contig = contig[1]
            fileout.write(f">contig_{i} len={length_contig}\n")
            sequence = textwrap.fill(seq_contig, width=80) + "\n"
            fileout.write(sequence)
            i = i + 1


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
