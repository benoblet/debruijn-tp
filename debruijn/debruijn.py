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

# First, specific functions
from operator import itemgetter
from random import randint
from itertools import product
# Second, global packages (ordered following pylint remarks)
import argparse
import os
import sys
import random
import statistics
import matplotlib
# Third, global packages with nicknames
import networkx as nx
# Forth, peculiar points
random.seed(9001)


__author__ = "Benedicte Noblet"
__copyright__ = "Universite Paris Diderot - Universite de Paris"
__credits__ = ["Benedicte Noblet"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Benedicte Noblet"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
    Parameters
    ----------
    path: string
        Path to the file
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
      Returns:
      	An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__,
    				      usage="{0} -h".format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Read fasta file and get sequences one at a time.

    Parameters
    ----------
    fastq_file: string
        Path to fastq file, third line with identifier or not

    Returns
    -------
    iterator
        An iterator operating on nucleotide sequences
    """
    with open(fastq_file, 'r') as myfile:
        for line in myfile:
            # pylint error: see github.com/PyCQA/pylint/issues/1999
            line = next(myfile, None)  # to prevent pylint message
            yield line.strip()
            line = next(myfile, None)
            line = next(myfile, None)


def cut_kmer(read, kmer_size):
    """Cut string in substrings of fixed length.

    Parameters
    ----------
    read: string
        Sequence read on Illumina or other platform
    kmer_size: integer
        Length k for k-mer (substring)

    Returns
    -------
    iterator
        An iterator operating over all substrings
        of length k in read sequence.
    """
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Create dictionnary of all read kmers.

    Parameters
    ----------
    fastq_file: string
        Path to fastq file
    kmer_size: integer
        Length k for k-mer (substring)

    Returns
    -------
    dictionnary
        A dictionnary containing all found kmers
        as keys and number of corresponding kmer
        as dictionnary value.
    """
    kmers_me = {}
    sequences = read_fastq(fastq_file)
    for read in sequences:
    	#print(f"current read : {read}")
        kmers = cut_kmer(read, kmer_size)
        for kmer in kmers:
            #print(f"kmer : {kmer}")
            if kmer in kmers_me:
                kmers_me[kmer] += 1
                #print("added")
            else:
                kmers_me[kmer] = 1
                #print("created")
    return kmers_me


def build_graph(kmer_dict):
    """Build oriented and weighted DeBruijn graph.

    Parameters
    ----------
    kmer_dict: dictionnary
        Dictionnary return by build_kmer_dict()

    Returns
    -------
    networkx DiGraph
        Graph with suffixes and prefixes of length
        k-1 as nodes linked by edges weighted with
        number of found sequences.
    """
    #print(kmer_dict)
    my_graph = nx.DiGraph()
    for kmer in kmer_dict:
    	#print(f"current kmer : {kmer}")
        prefix = kmer[:-1]
    	#print(f"prefixe : {prefix}")
        suffix = kmer[1:]
    	#print(f"suffixe : {suffix}")

    	#print(f"edge value : {kmer_dict[kmer]}")
        my_graph.add_edge(prefix, suffix, weight = kmer_dict[kmer])
    return my_graph


def remove_paths(graph, path_list,
                 delete_entry_node, delete_sink_node):
    """Remove given paths from a graph.

    Parameters
    ----------
    graph: networkx weighted DiGraph
        A network built by a function such as upper
        build_graph() function.
    path_list: list
        A list of tuple with ordered nodes of each path
        we want to delete from graph.
    delete_entry_node, delete_sink_node: booleans
        Value to set if starting (or ending respectively)
        node is to be kept in graph.

    Returns
    -------
    graph: networkx weighted DiGraph
        Input graph lacking paths (nodes and/or edges)
        from path_list.
    """
    # Store all nodes to prevent nodes used in several path
    # from being removed and following "path" using them in
    # path_list not being processed [deprecated comment]
    toberemoved = []
    for nodelist in path_list:
        print(f"current nodelist: {nodelist}")

        # test user wish for starting
        if delete_entry_node:
            print(f"first node to del: {nodelist[0]}")
            toberemoved.append(nodelist[0])
        # test user wish for ending
        if delete_sink_node:
            print(f"last node to del: {nodelist[-1]}")
            toberemoved.append(nodelist[-1])

        # identify central nodes if any
        if len(nodelist) > 2:
            print("Existing central nodes")
            for centralnode in nodelist[1:len(nodelist)-1]: # [1:-1] works also
                toberemoved.append(centralnode)

    # some cleaning of list
    toberemoved = list(set(toberemoved))
    # removal of remaining nodes
    graph.remove_nodes_from(toberemoved)

    return graph


def std(data):
    """Compute standard deviation of data.
    A call for stdev( ) function from statistic library.

    Parameters
    ----------
    data: list
        All float values to study for sample spread calculation.
    Returns
    -------
    float
        Computed standard deviation for values given in data.
    """
    if len(data) == 1:
        # I met an issue in simplify_bubbles
        # statistics.variance needs at least two values
        return 0

    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Identify longest or more frequent path among concurrent ones.

    Paths are parallel ones (i.e. starting and ending at same nodes)
    or competing ones (i.e. starting or ending at a same node).
    Here, we suppose that longest path or more present ones are more
    likely to occur in biological systems when performing reads
    assembly.

    Parameters
    ----------
    graph: networkx weighted DiGraph
        A network built by a function such as upper
        build_graph() function.
    path_list: list
        A list of lists or tuples with ordered nodes of each path.
        to study to compare.
    path_length: list
        A list of ordered path length values. They must fit path_list
        order.
    weight_avg_list: list
        A list of ordered path average weight values for each path
        in path_list. They must follow same order as path_list.
    delete_entry_node, delete_sink_node: booleans
        Value to set if starting (or ending respectively)
        node is to be kept in graph. See above: remove_path( ).

    Returns
    -------
    graph: networkx weighted DiGraph
        Input graph lacking redondant paths (nodes and/or edges) that
        where not the best found solution among existing ones.
    """
    # After simplify_bubbles's issue (10/10/2021 11:10 pm)
    # Do we have several paths in input path_list ?
    if len(path_list) <= 1:
        return graph

    # Do we have a more frequent/probable one ?
    weight_std = std(weight_avg_list)
    if weight_std > 0:
        if len(path_list) != len(weight_avg_list):
            sys.exit("Lists of paths and average weights do not seem" \
                     " to match. See select_best_path( ) function")
        maximum = max(weight_avg_list)
        for rank, path in enumerate(path_list):
            if weight_avg_list[rank] == maximum:
                best_one = path
                break
    else:
        # Do we have a longer one ?
        length_std = std(path_length)
        if length_std > 0:
            if len(path_list) != len(path_length):
                sys.exit("Lists of paths and associated lengths do" \
                         " not seem to match. See select_best_path" \
                         "( ) function")
            maximum = max(path_length)
            for rank, path in enumerate(path_list):
                if path_length[rank] == maximum:
                    best_one = path
                    break
        else:
            # Let's pick up one randomly, start and stop includ
            best_one = path_list[randint(0, len(path_list))]

    # remove chosen_one from list and delete all others
    path_list.remove(best_one)
    graph = remove_paths(graph, path_list,
                         delete_entry_node, delete_sink_node)

    return graph


def path_average_weight(graph, path):
    """Compute weight average for a path.
    The mean of weights allows to compare path of different
    lengths.

    Parameters
    ----------
    graph: networkx weighted DiGraph
        A network built by a function such as upper
        build_graph() function.
    path: tuple(graph, starting, ending)
        An item of the generated list by nx.all_simple_path
        containing all available paths between starting and
        ending (sink) nodes.

    Returns
    -------
    float
        Computed average weight value for edges of
        defined path in function parameter.
    """
    # from TP subject, release available on 10/09/2021 3:00 pm
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    """Remove given paths from a graph.

    Parameters
    ----------
    graph: networkx weighted DiGraph
        A network built by a function such as upper
        build_graph() function.
    ancestor_node, descendant_node: node
        Nodes label

    Returns
    -------
    graph: networkx weighted DiGraph
        Input graph with only one path between ancestor and
        descendant node.
    """
    all_paths = []
    all_length = []
    all_avgweit = []
    for path in nx.all_simple_paths(graph,
                                    ancestor_node, descendant_node):
        all_paths.append(path)
        all_length.append(len(path))
        all_avgweit.append(path_average_weight(graph, path))
    # modify graph based on above lists
    graph = select_best_path(graph, all_paths, all_length, all_avgweit,
                             delete_entry_node=False,
                             delete_sink_node=False)
    return graph


def simplify_bubbles(graph):
    """Remove multiple paths from a graph.
    See function select_best_path for more details on made choices.

    Parameters
    ----------
    graph: networkx weighted DiGraph
        A network built by a function such as upper
        build_graph() function.

    Returns
    -------
    graph: networkx weighted DiGraph
        Input graph with only no multiple paths.
    """
    bubble = False

    # Translation of algorithm given by teacher (written in French)
    for onenode in graph.nodes:
        if onenode not in graph.nodes:
            # we deleted it previously, no need to study it anymore
            continue
        # let's find all previous nodes linked by one edge
        predecessors = list(graph.predecessors(onenode))
        if len(predecessors) > 1:
            # either they are linked by a same node (i.e. bubble)
            # or they have different start node (competing paths).
            for i in range(len(predecessors)-1):
                one = predecessors[i]
                for j in range(i, len(predecessors)):
                    two = predecessors[j]
                    common = nx.lowest_common_ancestor(graph, one, two)
                    if common:
                        # damn, we've got bubble(s)
                        bubble = True
                        # we need to study it before going on
                        graph = simplify_bubbles(solve_bubble(graph, common, onenode))
                        #break
                # recursion error encountered when below lines outside loop:
                # https://stackoverflow.com/questions/52873067/recursionerror \
                # -maximum-recursion-depth-exceeded-in-comparison
                #if bubble:
                    # we call for solve_bubble to remove found existing bubble(s),
                    # and gives what the result to current function to find if
                    # there is(are) remaining bubble(s) to simplify.
                    

    return graph


def solve_entry_tips(graph, starting_nodes):
    # to remove pytest big error
    return graph


def solve_out_tips(graph, ending_nodes):
    # to remove pytest big error
    return graph


def get_starting_nodes(graph):
    """Find starting nodes.

    Parameters
    ----------
    graph: networkx DiGraph
        graph returned by build_graph() function.

    Returns
    -------
    List
        A list of all nodes without any predecessor.
    """
    starters = []
    for node in graph.nodes:
        predecessors = list(graph.predecessors(node))
        if len(predecessors) == 0:
            starters.append(node)
    return starters


def get_sink_nodes(graph):
    """Find terminal nodes.

    Parameters
    ----------
    graph: networkx DiGraph
        graph returned by build_graph() function.

    Returns
    -------
    List
        A list of all nodes without any successor.
    """
    stopsones = []
    for node in graph.nodes:
        successors = list(graph.successors(node))
        if len(successors) == 0:
            stopsones.append(node)
    return stopsones


def get_contigs(graph, starting_nodes, ending_nodes):
    """Find all contigs sequences.
    A contig is defined as a know

    Parameters
    ----------
    graph: networkx DiGraph
        Graph returned by build_graph() function.
    starting_nodes: list
        List of all nodes that have no predecessor.
    ending_nodes: list
        List of all nodes that have no successor.

    Returns
    -------
    List
        List of tuples (contig, length of contig).
    """
    my_contigs = []
    for starting, ending in product(starting_nodes, ending_nodes):
        #print(f"outer loop : {starting}, {ending}")
        #print(f"has path : {nx.has_path(graph, starting, ending)}")
        if not nx.has_path(graph, starting, ending):
            continue
        else:
            for singlepath in nx.all_simple_paths(graph, starting, ending):
                currentig = singlepath[0]
                for follower in singlepath[1:]:
                    currentig += follower[-1]
                my_contigs.append((currentig, len(currentig)))
    return my_contigs


def save_contigs(contigs_list, output_file):
    """Saved contigs in fasta format file.

    Exact format:
    >contig_1 len=longueur du contig
    Sequence..
    >contig_2 len=longueur du contig
    Sequence.

    Parameters
    ----------
    contigs_list: list
        List of tuple (contig sequence, length) as
        returned by get_contigs() function.
    output_file: string
        Filepath for output file.

    Returns
    -------
    None
        List of tuples (contig, length of contig).
    """
    with open(output_file, 'w') as file_out:
        for i, contig_infos in enumerate(contigs_list):
            sequence = contig_infos[0]
            length = contig_infos[1]
            file_out.write(f">contig_{i} len={length}\n")  # starts at 0 to fit predefined pytest
            file_out.write(f"{fill(sequence)}\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


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


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


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
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


#==============================================================
# Bene tests (Chabname's idea)
#==============================================================

def main_test(current=5):

    if current == 1:
        kmer_reader = cut_kmer("TCAGA", 3)
        for element in kmer_reader:
            print(element)
        sys.exit()

    args = get_arguments()

    ledico = build_kmer_dict(args.fastq_file, args.kmer_size)

    graph = build_graph(ledico)

    if current == 2:
        print("fetching start and ending nodes")
    begins = get_starting_nodes(graph)
    endings = get_sink_nodes(graph)

    if current == 2:
        print(f"{len(begins)} noeuds de debut")
        print(f"{len(endings)} noeuds de fin")
        print("calling get_contigs")
    my_contigs = get_contigs(graph, begins, endings)

    if current == 3:
        print("...saving...")
        save_contigs(my_contigs, args.output_file)

    if current == 4:
        print("computed average weight")
        onepath = my_contigs[1]
        print(f"my first contig : {onepath}")
        onemean = path_average_weight(graph, onepath)
        print(f"computed edges average : {onemean:.4f}")

    if current == 5:
        print("graph status before")
        print(graph)
        somepaths = list(product(begins[:5], endings[:5]))
        print(f"{len(somepaths)} nodes to remove")
        graph = remove_paths(graph, somepaths, True, True)
        print(f"new graph: {graph}")


if __name__ == '__main__':
#    main()
    main_test()
