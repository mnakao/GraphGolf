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

import argparse
import networkx as nx
import sys
import os
import subprocess

def main(infile, is_num, centers):
#        g = nx.read_edgelist(infile)
	g          = nx.read_edgelist(infile, nodetype=int)
	basename   = os.path.basename(infile.name)
	name, ext  = os.path.splitext(basename)
	image_name = name + ".png"
 	save_image(g, image_name, is_num, centers)
	print "Create " + image_name
	cmd = "open " + image_name
	subprocess.call(cmd, shell=True)
	return

def my_circular_layout(G, centers):
    import numpy as np

    if len(G) == 0:
        return {}

    twopi = 2.0*np.pi
    theta = np.arange(0, twopi, twopi/(len(G)-centers))
    pos   = np.column_stack([np.cos(theta), np.sin(theta)])
    for x in range(centers):
        pos   = np.insert(pos, len(G)-centers, [0, 0], axis=0)
    return dict(zip(G, pos))

def save_image(g, filepath, is_num, centers):
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	
	if centers != 0:
		layout = my_circular_layout(g, centers)
	else:
		layout = nx.circular_layout(g)

	nx.draw(g, with_labels=is_num, node_size=50, linewidths=0, alpha=0.5, node_color='#3399ff', edge_color='#666666', pos=layout)
	plt.draw()
	plt.savefig(filepath)
	return

def parser():
	usage = 'Usage: python {} FILE [--num] [--center] [--help]'.format(__file__)
	argparser = argparse.ArgumentParser(usage=usage)
	argparser.add_argument('infile', type=argparse.FileType('r'), default=sys.stdin, help='Input file')
	argparser.add_argument('-n', '--num',    action='store_true', help='write node number')
	argparser.add_argument('-c', '--center', type=int, help='Number of centers')
	return argparser.parse_args()

if __name__ == '__main__':
	args = parser()
	main(args.infile, args.num, args.center)
