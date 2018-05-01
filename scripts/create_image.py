#!/usr/bin/python2.7
# coding: utf-8
#=======================================================================
#
#	Create a random graph
#	by Ikki Fujiwara, National Institute of Informatics
#	2015-06-23
#
#=======================================================================
# "create-random.py" is licensed under a Creative Commons Attribution 4.0 International License.
# http://creativecommons.org/licenses/by/4.0/

import networkx as nx
import argparse
import sys
import os

argumentparser = argparse.ArgumentParser()
argumentparser.add_argument('infile', type=argparse.FileType('r'), default=sys.stdin)

def main(args):
        infile = args.infile
#        g = nx.read_edgelist(infile)
	g = nx.read_edgelist(infile, nodetype=int)
	basename = os.path.basename(infile.name)
	name, ext = os.path.splitext(basename)
	image_name = name + ".png"
 	save_image(g, image_name)
	print "Create " + image_name
	return

def save_image(g, filepath):
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	
	layout = nx.circular_layout(g)
#	nx.draw(g, with_labels=False, node_size=50, linewidths=0, alpha=0.5, node_color='#3399ff', edge_color='#666666', pos=layout)
	nx.draw(g, with_labels=True, node_size=50, linewidths=0, alpha=0.5, node_color='#3399ff', edge_color='#666666', pos=layout)
	plt.draw()
	plt.savefig(filepath)
	return

if __name__ == '__main__':
	main(argumentparser.parse_args())
