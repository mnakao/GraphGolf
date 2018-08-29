#!/usr/bin/python2.7
# coding: utf-8

import networkx as nx
import sys

def main(filename):
  infile = open(filename, "r")

  num_of_lines   = 0
  max_vertex_num = 0
  for line in infile:
    num_of_lines += 1
    itemList = line[:-1].split(' ')
    max_item = int(max(itemList))
    if max_item > max_vertex_num:
      max_vertex_num = max_item

  nnodes = max_vertex_num + 1
  assert((num_of_lines*2) % nnodes == 0)
  degree = num_of_lines * 2 / nnodes
  print "N = ", nnodes, " Degree = ", degree

  infile.seek(0)
  g = nx.read_edgelist(infile)
  if nx.is_connected(g):
    hops = nx.shortest_path_length(g, weight=None)
    diam, aspl = max_avg_for_matrix(hops)
  
  print("Diam. = {}\t ASPL = {}".format(diam, aspl))

  low_diam, low_aspl = lower_bound_of_diam_aspl(nnodes, degree)
  print(": Lower Diam. = {}\t ASPL      = {}".format(low_diam, low_aspl))
  print(": Diam. Gap.  = {}\t ASPL Gap. = {}".format(diam-low_diam, aspl-low_aspl))

def lower_bound_of_diam_aspl(nnodes, degree):
  diam = -1
  aspl = 0.0
  n = 1
  r = 1
  while True:
    tmp = n + degree * pow(degree - 1, r - 1)
    if tmp >= nnodes:
      break
    n = tmp
    aspl += r * degree * pow(degree - 1, r - 1)
    diam = r
    r += 1
  diam += 1
  aspl += diam * (nnodes - n)
  aspl /= (nnodes - 1)
  return diam, aspl

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

  print "SUM =", sum / 2
  return max, sum / cnt

if __name__ == '__main__':
  params = sys.argv
  if len(params) != 2:
    print 'Usage: python %s parameter' % params[0]
    quit()
  
  filename = params[1]
  main(filename)
