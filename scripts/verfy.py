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
import time
argumentparser = argparse.ArgumentParser()
argumentparser.add_argument('infile', type=argparse.FileType('r'), default=sys.stdin)

def main(args):
	infile = args.infile
	g = nx.read_edgelist(infile)
	if nx.is_connected(g):
		start = time.time()
		hops = nx.shortest_path_length(g, weight=None)
		diam, aspl = max_avg_for_matrix(hops)
		elapsed_time = time.time() - start
		print ("elapsed_time:{0}".format(elapsed_time)) + "[sec]"
	
	print("Diam. k = {}\t ASPL l = {}".format(diam, aspl))
	return

def max_avg_for_matrix(data):
	cnt = 0
	sum = 0.0
	max = 0.0
	for i in data:
		for j in data[i]:
			if i != j:
				cnt += 1
				sum += data[i][j]
				if max < data[i][j]:
					max = data[i][j]
	return max, sum / cnt

if __name__ == '__main__':
	main(argumentparser.parse_args())
