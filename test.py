from collections import Counter

import pandas as pd
import os

from tqdm import tqdm, trange


class node:
    def __init__(self, kmer):
        self.sequence = kmer
        self.k = len(kmer)
        self.next = []
        self.previous = []

    def getNextMatch(self):
        return self.sequence[1:]

    def getPreviousMatch(self):
        return self.sequence[:-1]

    def addNext(self,node):
        self.next.append(node)

    def addPrevious(self,node):
        self.previous.append(node)

    def getKmer(self):
        return self.sequence

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



k = 4
sequences = []
for root, dir, files in os.walk('avastin/avastin'):
    root = root + '/'
    for file in files:
        filename = root + file
        data = pd.read_csv(filename, delimiter='\t')
        temp = data[data['Score']>0.5]
        temp = temp[-50<temp['PPM Difference']]
        temp = temp[temp['PPM Difference']<50]
        sequences.extend(temp['DENOVO'].values)

kmers_dic = {}
for sequence in sequences:
    kmer = get_kmer_count_from_sequence(sequence, k=k)
    kmers_dic = dict(Counter(kmers_dic) + Counter(kmer))
print('共有 {} 个kmer(k = {})！！！！'.format(len(kmers_dic),k))


nodes = []
for kmer in kmers_dic:
    nodes.append(node(kmer))

for i in trange(len(nodes)):
    for node_1 in nodes:
        for node_2 in nodes:
            if node_1.sequence != node_2.sequence:
                if node_1.sequence[1:] == node_2.sequence[:-1]:
                    node_1.addNext(node_2)
                if node_1.sequence[:-1] == node_2.sequence[1:]:
                    node_1.addPrevious(node_2)




