from collections import Counter

import pandas as pd
import os

import toyplot
import debruijn
from tqdm import tqdm, trange




def get_kmer_count_from_sequence(sequence, k=3):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # dict to store kmers
    kmers = {}

    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]

        if len(kmer) != k:
            break

        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    return kmers


def get_debruijn_edges_from_kmers(kmers):
    """
    Every possible (k-1)mer (n-1 suffix and prefix of kmers) is assigned
    to a node, and we connect one node to another if the (k-1)mer overlaps
    another. Nodes are (k-1)mers, edges are kmers.
    """
    # store edges as tuples in a set
    edges = set()

    # compare each (k-1)mer
    for k1 in kmers:
        for k2 in kmers:
            if k1 != k2:
                # if they overlap then add to edges
                if k1[1:] == k2[:-1]:
                    edges.add((k1[:-1], k2[:-1]))
                if k1[:-1] == k2[1:]:
                    edges.add((k2[:-1], k1[:-1]))

    return edges

def plot_debruijn_graph(edges, width=500, height=500):
    "returns a toyplot graph from an input of edges"
    graph = toyplot.graph(
        [i[0] for i in edges],
        [i[1] for i in edges],
        width=width,
        height=height,
        tmarker=">",
        vsize=25,
        vstyle={"stroke": "black", "stroke-width": 2, "fill": "none"},
        vlstyle={"font-size": "11px"},
        estyle={"stroke": "black", "stroke-width": 2},
        layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
    return graph

k = 3
sequences = []
sequences_scores = []
for root, dir, files in os.walk('avastin/avastin'):
    root = root + '/'
    for file in files:
        filename = root + file
        data = pd.read_csv(filename, delimiter='\t')
        temp = data[data['Score']>0.5]
        temp = temp[-50<temp['PPM Difference']]
        temp = temp[temp['PPM Difference']<50]
        temp.reset_index(inplace=True)
        sequences.extend(temp['DENOVO'].values)
        for i in range(len(temp)):
            sequences_scores.append([temp['DENOVO'][i],temp['Score'][i]])

inputs1 = sequences_scores
delList = []

for i in inputs1:
    if 'LLLLLL' in i[0]:
        delList.append(i)
    elif 'PPPPPP' in i[0]:
        delList.append(i)
    elif 'GGGGGG' in i[0]:
        delList.append(i)
    elif 'SSSSSS' in i[0]:
        delList.append(i)
    elif 'KKKKKK' in i[0]:
        delList.append(i)

for i in delList:
    inputs1.remove(i)

inputsBefore = []

for i in inputs1:
    if len(i[0]) > 4:
        inputsBefore.append(i)
    # print(222222222222,inputsBefore)
    # graphkmer = debruijn.DeBruijnGraph(inputsBefore, 5)
    # pathkmer = graphkmer.longestPath()
    # removeList = []
    # inputsBefore = []

    # for i in inputsBefore:
    #     for j in pathkmer:
    #         if i[0] in j[1][0]:
    #             removeList.append(i)

    # for i in removeList:
    #     inputsBefore.remove(i)

pathkmer1 = []
lengthList = []
print(len(inputsBefore))
for a in trange(8):

    inputskmer = []

    for i in inputsBefore:
        if len(i[0]) > (3+a):
            inputskmer.append(i)
    # print(33333333333,inputskmer)
    graphkmer = debruijn.DeBruijnGraph(inputskmer, 3+a)
    pathkmer = graphkmer.longestPath()

    removeList = []
    inputsBefore = []

    for i in inputskmer:
        for j in pathkmer:
            if i[0] in j[1][0]:
                removeList.append(i)

    for i in removeList:
        inputskmer.remove(i)

    for i in inputskmer:
        inputsBefore.append(i)

    for i in pathkmer:
        inputsBefore.append(i[1])

    for i in pathkmer:
        lengthList.append(i[0])

    for i in pathkmer:
        pathkmer1.append([i[0],i[1][0]])


pathkmer1.sort(key=lambda x:x[0])
print(pathkmer1)
# kmers_dic = {}
# for sequence in sequences:
#     kmer = get_kmer_count_from_sequence(sequence, k=k)
#     kmers_dic = dict(Counter(kmers_dic) + Counter(kmer))
# print('共有 {} 个kmer(k = {})！！！！'.format(len(kmers_dic),k))





