#! /usr/bin/env python
from __future__ import print_function
import sys

class Family:
    def __init__(self, samples, pedigree):
        self.family_id = pedigree
        self.generations = Family.getGenerations(samples)

        self.f1_male = '????'
        self.f1_female = '????'
        self.p0_male_paternal = '????'
        self.p0_female_paternal = '????'
        self.p0_male_maternal = '????'
        self.p0_female_maternal = '????'

        for f1 in self.generations['F1']:
            #if paternal, set paternal p0s
            if f1[4] == '1':
                self.f1_male = f1[1]
                if f1[2] != '0':
                    self.p0_male_paternal = f1[2]
                if f1[3] != '0':
                    self.p0_female_paternal = f1[3]
            #if maternal, set maternal p0s
            else:
                self.f1_female = f1[1]
                if f1[2] != '0':
                    self.p0_male_maternal = f1[2]
                if f1[3] != '0':
                    self.p0_female_maternal = f1[3]


    @staticmethod
    def getGenerations(samples):
        generations = {'P0':[], 'F1':[], 'F2':[]}
        
        #scan for grandparents (P0s)
        #using the assumption that the grandparents 
        #will be the only samples without IDs for either parent
        for sample in samples:
            if sample[2] == '0' and sample[3] == '0':
                generations['P0'].append(sample)
        #scan for parents (F1s)
        #using the assumption that they will have parent IDs in the P0s list
        for sample in samples:
            for p0 in generations['P0']:
                #if the sample's dad or mom is the current p0, it's an f1
                if (sample[2] == p0[1] or sample[3] == p0[1]) and sample not in generations['F1']:
                    generations['F1'].append(sample)
        #scan for kids (F2s)
        #using the assumption that they will have parent IDs in the F1s list 
        for sample in samples:
            for f1 in generations['F1']:
                if (sample[2] == f1[1] or sample[3] == f1[1]) and sample not in generations['F2']:
                    generations['F2'].append(sample)
        return generations

    @staticmethod
    def extractFamiliesFromPedigree(samples):
        generations = Family.getGenerations(samples)
        
        families = {}
        for f2 in generations['F2']:
            #family_id is the ped id and both F1s
            family_id = f2[0] + "_" + f2[2]+ "_" + f2[3]
            if family_id not in families:
                families[family_id] = []
            families[family_id].append(f2)

        for family_id in families:
            parent_ids = family_id.split('_')[1:]
            parents = []
            for f1 in generations['F1']:
                if f1[1] in parent_ids:
                    parents.append(f1)
            for f1 in parents:
                for p0 in generations['P0']:
                    if p0[1] in f1[2:4]:
                        families[family_id].append(p0)
            families[family_id] += parents
        return families

    def __str__(self):
        output = self.family_id + "\n"
        output += "Generations:\n"
        for gen in self.generations:
            output += gen + ":\n"
            for sample in self.generations[gen]:
                output += "\t".join(sample) + "\n"
        return output



    def pretty_print(self):
        overall_length = 8*len(self.generations['F2'])-4

        top_start = self.p0_male_paternal + '    ' + self.p0_female_paternal
        top_end = self.p0_male_maternal + '    ' + self.p0_female_maternal
        gap_len = max(overall_length - (len(top_start) + len(top_end)),4)
        top_str = top_start + " "*gap_len + top_end

        spacer1 = "\n  |______|" + " "*(gap_len+4) + "|______|  "
        spacer1 += "\n      |" + " "*(gap_len+10) + "|      \n"

        middle_start = "    " + self.f1_male
        middle_end = self.f1_female + "    "
        gap_len = overall_length - (len(middle_start) + len(middle_end))
        middle_str = middle_start + " "*gap_len + middle_end
        
        spacer2 = "\n  ____|" + "_"*(gap_len+2) + "|____  \n  "
        for i in range(len(self.generations['F2'])):
            if i%2==0:
                spacer2 += "|       " 
            else:
                spacer2 += "|       "
        spacer2 += "\n"

        
        return "Pedigree " + self.family_id + ":\n" \
                + top_str + spacer1 \
                + middle_str + spacer2 \
                + "    ".join([f2[1] for f2 in self.generations['F2']]) + "\n"


def CreateDictFromPED(ped, remove_failed = False):
    pedigrees = {}
    with open(ped, 'r') as ped_file:
        for line in ped_file:
            fields = line.strip().split()
            pedigree = fields[0]
            if pedigree not in pedigrees:
                pedigrees[fields[0]] = []
            pedigrees[pedigree].append(fields)

    families = {}
    for pedigree in pedigrees:
        sub_pedigrees = Family.extractFamiliesFromPedigree(pedigrees[pedigree])
        for family_id in sub_pedigrees:
            family = Family(sub_pedigrees[family_id], pedigree)
            families[family_id] = family
    return families
