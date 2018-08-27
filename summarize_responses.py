#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-r", "--responses", help="tab-separated file of responses (from retrieval.py",required=True)
args = parser.parse_args()

imgs_scored = {}
keys = {}
with open(args.responses, 'r') as responses_file:
    question = responses_file.readline().strip()
    answers = responses_file.readline().strip().split('\t')
    keys_list = responses_file.readline().strip().split('\t')
    for i in range(len(keys_list)):
        keys[keys_list[i].replace(" ", "_")] = i

    for line in responses_file:
        fields = line.strip().split("\t")
        img = fields[keys['IMAGE']]
        if img not in imgs_scored:
            imgs_scored[img] = []
        imgs_scored[img].append(fields)

#replace this with arg
score_map = {'Supports':1,
        "Does_not_support":0,
        "Review_(IGV)": 0}
print ("#CHROM\tSTART\tEND\tSCORE\tCOUNT\tSVTYPE\tIMAGE")
for img in imgs_scored:
    summary = []
    score = 0.0
    count = len(imgs_scored[img])
    for entry in imgs_scored[img]:
        if len(summary) == 0:
            summary = [entry[keys['chrom']],
                        entry[keys['start']],
                        entry[keys['end']],
                        0,0,
                        entry[keys['sv_type']],
                        entry[keys['IMAGE']]]
        text_score = (entry[keys['SCORE']]).replace(" ", "_")
        score += score_map[text_score]
    summary[3] = str(score/count)
    summary[4] = str(count)
    print ("\t".join(summary))
