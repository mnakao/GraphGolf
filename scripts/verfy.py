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
  closed_circuit_3 = count_closed_circuit(nnodes, degree, 3, infile)
  closed_circuit_4 = count_closed_circuit(nnodes, degree, 4, infile)
  print "Num. of closed curcuit (length is 3) : ", closed_circuit_3
  print "Num. of closed curcuit (length is 4) : ", closed_circuit_4

  infile.seek(0)
  g = nx.read_edgelist(infile)
  if nx.is_connected(g):
    hops = nx.shortest_path_length(g, weight=None)
    diam, aspl = max_avg_for_matrix(hops)
  
  print("Diam. = {}\t ASPL = {}".format(diam, aspl))

  low_diam, low_aspl = lower_bound_of_diam_aspl(nnodes, degree)
  print(": Lower Diam. = {}\t ASPL      = {}".format(low_diam, low_aspl))
  print(": Diam. Gap.  = {}\t ASPL Gap. = {}".format(diam-low_diam, aspl-low_aspl))

def count_closed_circuit(nnodes, degree, distance, infile):
  adjacency = [[0 for i in range(degree)] for j in range(nnodes)]
  count     = [0 for i in range(nnodes)]

  infile.seek(0)
  for line in infile:
    itemList = line[:-1].split(' ')
    n1 = int(itemList[0])
    n2 = int(itemList[1])
    adjacency[n1][count[n1]] = n2
    adjacency[n2][count[n2]] = n1
    count[n1]+=1
    count[n2]+=1
  
  count = 0
  if distance == 3:
    for v1 in range(nnodes):
      for i in range(degree):
        v2 = adjacency[v1][i]
        for j in range(i+1, degree):
          v3 = adjacency[v1][j]
          for d in range(degree):
            if v3 == adjacency[v2][d]:
              count += 1

    return count/3
  elif distance == 4:
    for v1 in range(nnodes):
      for i in range(degree):
        v2 = adjacency[v1][i]
        for j in range(i+1, degree):
          v3 = adjacency[v1][j]
          for v4 in range(v1+1, nnodes):
            if(v4 == v2 or v4 == v3):
              continue
            for d1 in range(degree):
              for d2 in range(d1+1,degree):
                if((adjacency[v4][d1] == v2 and adjacency[v4][d2] == v3) or
                   (adjacency[v4][d1] == v3 and adjacency[v4][d2] == v2)):
                  count += 1
    return count/2
  else:
    return -1

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

  return max, sum / cnt

if __name__ == '__main__':
  params = sys.argv
  if len(params) != 2:
    print 'Usage: python %s parameter' % params[0]
    quit()
  
  filename = params[1]
  main(filename)
