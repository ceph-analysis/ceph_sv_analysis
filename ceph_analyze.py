#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse

parser = argparse.ArgumentParser(description="analyze the curation of CEPH de novo variants")
parser.add_argument("-s", "--summary", help="summary file of the variants scored",required=True)
args = parser.parse_args()

variants = {}

with open(args.summary, 'r') as summary_file:
    for line in summary_file:
        if line[0] != '#':
            fields = line.strip().split()
            position = "_".join(fields[:3])
            if position not in variants:
                variants[position] = []
            variants[position].append(fields)

def all_same(scores_list):
    if len(scores_list) == 0: 
        return True
    first_score = scores_list[0]
    for score in scores_list:
        if score != first_score:
            return False
    return True


#quality: good (de novo), bad (not de novo), unknown (mix)
for position in variants:
    quality = 0
    scores = []
    for entry in variants[position]:
        score = float(entry[3])
        scores.append(score)
    if all_same(scores):
        quality = scores[0]
    else:
        quality = -1
    result = position.replace("_", "\t") + "\t" + str(int(quality))
    print (result)
