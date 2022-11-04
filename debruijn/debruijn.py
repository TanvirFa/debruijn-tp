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

__author__ = "Tanvir Faruque"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Tanvir Faruque"
__email__ = "tanvir.faruque@etu.u-paris.fr"
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
	with open (fastq_file) as file:
		for i in file:
			yield next(file).strip("\n")
			next(file)
			next(file)


def cut_kmer(read,kmer_size):
	count = 0
	while((count + kmer_size) <= len(read)):
		yield read[count:(count + kmer_size)]
		count += 1


def build_kmer_dict(fastq_file, kmer_size):
	dic_kmer = {}
	for i in read_fastq(fastq_file):
		for j in cut_kmer(i,kmer_size):
			if j in dic_kmer:
				dic_kmer[j] += 1
			else:
				dic_kmer.update({j:1})
	return dic_kmer


def build_graph(kmer_dict):
	graph_kmer = nx.DiGraph()
	for i in kmer_dict:
		graph_kmer.add_edge(i[:-1],i[1:],weight = kmer_dict[i])
	return graph_kmer


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
	for p in path_list:
		if(delete_entry_node and delete_sink_node):
			graph.remove_nodes_from(p) 
		elif(delete_entry_node):
			graph.remove_nodes_from(p[:-1]) 
		elif(delete_sink_node):
			graph.remove_nodes_from(p[1:]) 
		else:
			graph.remove_nodes_from(p[1:-1]) 
	return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node = False, delete_sink_node = False):
    if(statistics.stdev(weight_avg_list) > 0):
        ind = weight_avg_list.index(max(weight_avg_list))
        path_list.pop(ind) 
        graph = remove_paths(graph,path_list,delete_entry_node, delete_sink_node)
    elif(statistics.stdev(path_length) > 0):
        ind = path_length.index(max(path_length))
        path_list.pop(ind)
        graph = remove_paths(graph,path_list,delete_entry_node, delete_sink_node)
    else:
        random.seed(9001)
        r = random.randint(0, len(path_list))
        path_list.pop(r)
        graph = remove_paths(graph,path_list,delete_entry_node, delete_sink_node)
    return graph

def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
	path_list = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
	path_lenght = []
	weight_avg_list = []
	for i in path_list:
		lw = []
		path_lenght.append(len(i))
		l = list(graph.subgraph(i).edges(data = True))
		for j in l:
			lw.append(j[2]['weight'])
		weight_avg_list.append(statistics.mean(lw))
	graph = select_best_path(graph, path_list, path_lenght, weight_avg_list,delete_entry_node = False, delete_sink_node = False)
	return graph

def simplify_bubbles(graph):
	bubble = False 
	for i in graph.nodes:
		noeud_n = i
		liste_predecesseurs = list(graph.predecessors(i))
		if len(liste_predecesseurs) > 1:
			for j in liste_predecesseurs :
				for pre in liste_predecesseurs:
					noeud_ancetre = nx.lowest_common_ancestor(graph,pre,j)
					if noeud_ancetre != None:
							bubble = True
							break
		if bubble:
			break
	# La simplification ayant pour conséquence de supprimer des noeuds du hash
	# Une approche récursive est nécessaire avec networkx
	if bubble:				
		graph = simplify_bubbles(solve_bubble(graph,noeud_ancetre, noeud_n))
	return graph

def drop_entry_tips(graph,list_prec,noeud_n):
	path_lenght = []
	weight_avg_list = []
	path_list = []
	#len(path_list)=1 car il n'y a plus de bulle
	for ancestor_node in list_prec:
		path_list.append(list(nx.all_simple_paths(graph, source=ancestor_node, target = noeud_n))[0])
		for i in path_list:
			lw = []
			path_lenght.append(len(i))
			l = list(graph.subgraph(i).edges(data=True))
			for j in l:
				lw.append(j[2]['weight'])
			weight_avg_list.append(statistics.mean(lw))
	graph = select_best_path(graph, path_list, path_lenght, weight_avg_list,delete_entry_node = True, delete_sink_node = False)
	return graph

def drop_out_tips(graph,list_successor,noeud_n):
	path_lenght = []
	weight_avg_list = []
	path_list = []
	#len(path_list)=1 car il n'y a plus de bulle
	for ancestor_node in list_successor:
		path_list.append(list(nx.all_simple_paths(graph, source=noeud_n, target = list_successor))[0])
		for i in path_list:
			lw = []
			path_lenght.append(len(i))
			l = list(graph.subgraph(i).edges(data=True))
			for j in l:
				lw.append(j[2]['weight'])
			weight_avg_list.append(statistics.mean(lw))
	graph = select_best_path(graph, path_list, path_lenght, weight_avg_list,delete_entry_node = False, delete_sink_node = True)
	return graph

def solve_entry_tips(graph, starting_nodes):
	tips = False 
	for i in graph.nodes:
		noeud_n = i
		count = 0
		prec_list = []
		liste_predecesseurs = list(graph.predecessors(i))
		if len(liste_predecesseurs) > 1:
			for j in liste_predecesseurs :
				if j in starting_nodes:
					count += 1
					prec_list.append(j)
					print(prec_list)
			if count > 1:
				tips = True
				break
		if tips:
			break
	if tips:
	  graph = drop_entry_tips(graph,prec_list,noeud_n)		
	  graph = solve_entry_tips(graph,starting_nodes)
	return graph

def solve_out_tips(graph, ending_nodes):
	tips = False 
	for i in graph.nodes:
		noeud_n = i
		count = 0
		successor_list = []
		liste_successeurs = list(graph.successors(i))
		if len(liste_successeurs) > 1:
			for j in liste_successeurs :
				if j in ending_nodes:
					count+=1
					successor_list.append(j)
			if count > 1:
				tips = True
				break
		if tips:
			 break
	if tips:
		graph = drop_out_tips(graph,successor_list,noeud_n)		
		graph = solve_out_tips(graph,ending_nodes)
	return graph


def get_starting_nodes(graph):
	l_first_nodes = []
	for i in graph.nodes:
		it = list(graph.predecessors(i))
		if (len(it) == 0):
			l_first_nodes.append(i)
	return l_first_nodes

def get_sink_nodes(graph):
	l_last_nodes = []
	for i in graph.nodes:
		it = list(graph.successors(i))
		if (len(it) == 0):
			l_last_nodes.append(i)
	return l_last_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
	list_contig = []
	seq = ""
	for i in starting_nodes:
		for j in ending_nodes:
			if(nx.has_path(graph, i, j)):
				for k in nx.all_simple_paths(graph, source = i, target = j):
					seq += i
					for kmer in k[1:]:
						seq += kmer[-1]
					list_contig.append((seq,len(seq)))
					seq = ""
	return list_contig
		

def save_contigs(contigs_list, output_file):
	file=open(output_file+".fasta","w")
	count = 0
	for contig in contigs_list:
		file.write(textwrap.fill(">contig_"+str(count)+" len="+str(contig[1]),width=80))
		file.write("\n")
		file.write(textwrap.fill(contig[0], width=80))
		file.write("\n")
	file.close()


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
    #en arguments le fichier input et nom fichier output
    kmer_size=10
    graph=build_graph(build_kmer_dict(args[0],kmer_size))
    graph=simplify_bubbles(graph)
    start_nodes=get_starting_nodes(graph)
    end_nodes=get_sink_nodes(graph)
    graph=solve_entry_tips(graph,start_nodes)
    graph=solve_out_tips(graph,end_nodes)
    save_contigs(get_contigs(graph,start_nodes,end_nodes),args[1])

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
