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

def main(infile, is_num, is_center):
#        g = nx.read_edgelist(infile)
	g          = nx.read_edgelist(infile, nodetype=int)
	basename   = os.path.basename(infile.name)
	name, ext  = os.path.splitext(basename)
	image_name = name + ".png"
 	save_image(g, image_name, is_num, is_center)
	print "Create " + image_name
	cmd = "open " + image_name
	subprocess.call(cmd, shell=True)
	return

def my_circular_layout(G):
    import numpy as np

    if len(G) == 0:
        return {}

    twopi = 2.0*np.pi
    theta = np.arange(0, twopi, twopi/(len(G)-1))
    pos   = np.column_stack([np.cos(theta), np.sin(theta)])
    pos   = np.insert(pos, len(G)-1, [0, 0], axis=0)
    return dict(zip(G, pos))

def save_image(g, filepath, is_num, is_center):
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	
	if is_center:
		layout = my_circular_layout(g)
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
	argparser.add_argument('-c', '--center', action='store_true', help='Node number 0 is located at the center')
	return argparser.parse_args()

if __name__ == '__main__':
	args = parser()
	main(args.infile, args.num, args.center)
