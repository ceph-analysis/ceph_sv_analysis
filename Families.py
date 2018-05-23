#! /usr/bin/env python
from __future__ import print_function
import sys

class Family:
#    family_id = ''
#    mother = []
#    father = []
#    mgrandmother = []
#    mgrandfather = []
#    pgrandmother = []
#    pgrandfather = []
#    kids = []
#    all_members = []

    def __init__(self, samples, pedigree):
        #add to family
        self.all_members = samples
        
        
        #scan for grandparents, using the assumption that the grandparents 
        #will be the only samples without parent IDs
        grandparents = []
        for sample in samples:
            if sample[2] == '0' and sample[3] == '0':
                grandparents.append(sample)
        #scan for parents
        parents = []
        for sample in samples:
            for grandparent in grandparents:
                if sample[2] == grandparent[1]:
                    parents.append(sample)
#        if len(grandparents) != 4:
#            #print ("missing grandparent in pedigree: " + pedigree)
#            return None
#        if len(parents) != 2:
#            #print ("missing parent in pedigree: " + pedigree)
#            return None
        try:
            if parents[0][4] == '1':
                self.father = parents[0]
                self.mother = parents[1]
            if parents[1][4] == '1':
                self.father = parents[1]
                self.mother = parents[0]
        except:
            print ("failed to add a parent")
        
        try:
            for gparent in grandparents:
                if gparent[1] == self.father[2]:
                    self.pgrandfather = gparent
                elif gparent[1] == self.father[3]:
                    self.pgrandmother = gparent
                elif gparent[1] == self.mother[2]:
                    self.mgrandfather = gparent
                elif gparent[1] == self.mother[3]:
                    self.mgrandmother = gparent
        except:
            print ("failed to add a grandparent")

        self.kids = []
        for sample in samples:
            if sample not in grandparents and sample not in parents:
                self.kids.append(sample)
        self.family_id = pedigree

#    def __init__(self, samples, pedigree):
#        #add to family
#        self.all_members = samples        
#        
#        #scan for grandparents, using the assumption that the grandparents 
#        #will be the only samples without parent IDs
#        self.generations = []
#        self.couples = []
#        self.siblings_by_mom = {}
#        p0 = []
#        for sample in samples:
#            sample_dad = sample[2]
#            sample_mom = sample[3]
#            if sample_mom not in self.siblings_by_mom:
#                self.siblings_by_mom[sample_mom] = []
#            self.siblings_by_mom[sample_mom].append(sample[1])
#            self.couples.append([sample_dad, sample_mom])
#
#            if sample[2] == '0' and sample[3] == '0':
#                print (sample)
#                p0.append(sample)
#        self.generations.append(p0)
#
#        #scan to place samples in generations
#        count_included = len(p0)
#        current_gen_pos = 0
#        while count_included < len(self.all_members):
#            current_gen_samples = []
#            current_gen_pos += 1
#            for sample in samples:
#                for prev_gen_sample in self.generations[current_gen_pos-1]:
#                    if sample[2] == prev_gen_sample[1]:
#                        current_gen_samples.append(sample)
#                        count_included += 1
#            self.generations.append(current_gen_samples)
#        self.family_id = pedigree

    def get_couple_by_sample_id(sample_id):
        for couple in self.couples:
            if sample_id in couple:
                return couple
        return sample_id

    def couplize_generation(generation):
        couples = set()
        for sample in generation:
            couples.add(self.get_couple_by_sample_id(sample))
        return [couples]
    
    def get_gen_by_sample_id(sample_id):
        for i in range(len(self.generations)):
            if sample_id in self.generations[i]:
                return i
        return -1

#    def pretty_print(self):
#        if len(self.generations > 1):
#            print ("prety-printing of less than 2-generation pedigrees is not yet supported")
#            sys.exit(1)
#        if len(self.generations > 3):
#            print ("prety-printing of more than 3-generation pedigrees is not yet supported")
#            sys.exit(1)
#
#        for generation  in self.generations:
#            gen_couples = self.coupleize_generation(generation)
#            for couple in gen_couples:
#                
#
#        first_gen_string = 
#        if len(self.generations) == 1:
#            print ("\t".join(sel.generations[0]))
#        elif len(self.generations) == 2:
#
#        
#        for sample in 
#
#        def create_line(generation):
#            spacer = ''
#            line = ''
#            odd = False
#            for sample in generation:
#                spacer += "  |   "
#                if odd:
#                    spacer += "  "
#                    odd = !odd
#                line += sample[1] + "    "
#        overall_length = 8*len(self.generations[-1])-4
#
#        top_start = self.pgrandfather[1] + '    ' + self.pgrandmother[1]
#        top_end = self.mgrandfather[1] + '    ' + self.mgrandmother[1]
#        gap_len = max(overall_length - (len(top_start) + len(top_end)),4)
#        top_str = top_start + " "*gap_len + top_end
#
#        spacer1 = "\n  |______|" + " "*(gap_len+4) + "|______|  "
#        spacer1 += "\n      |" + " "*(gap_len+10) + "|      \n"
#
#        middle_start = "    " + self.father[1]
#        middle_end = self.mother[1] + "    "
#        gap_len = overall_length - (len(middle_start) + len(middle_end))
#        middle_str = middle_start + " "*gap_len + middle_end
#        
#        spacer2 = "\n  ____|" + "_"*(gap_len+2) + "|____  \n  "
#        
#        for i in range(len(self.kids)):
#            if i%2==0:
#                spacer2 += "|       " 
#            else:
#                spacer2 += "|       "
#        spacer2 += "\n"
#
#        
#        return "Pedigree " + self.family_id + ":\n" \
#                + top_str + spacer1 \
#                + middle_str + spacer2 \
#                + "    ".join([kid[1] for kid in self.kids]) + "\n"

    def test(self):
        passed_test = True
        error_msg = []
        # step one: make sure generations are there
        if not len(self.mgrandmother) > 0:
            error_msg.append("Missing maternal grandmother")
            passed_test = False
        if not len(self.mgrandfather) > 0:
            error_msg.append("Missing maternal grandfather")
            passed_test = False
        if not len(self.pgrandmother) > 0:
            error_msg.append("Missing paternal grandmother")
            passed_test = False
        if not len(self.pgrandfather) > 0:
            error_msg.append("Missing paternal grandfather")
            passed_test = False
        if not len(self.mother) > 0:
            error_msg.append("Missing mother")
            passed_test = False
        if not len(self.father) > 0:
            error_msg.append("Missing father")
            passed_test = False

        if not passed_test:
            error_msg.append("Family " + self.family_id + " failed test")
            return error_msg

        # step two: make sure relationships match
        for kid in self.kids:
            if kid[2] != self.father[1]:
                passed_test = False
                error_msg.append("Kid " + kid[1] +  " doesn't match father")
            if kid[3] != self.mother[1]:
                passed_test = False
                error_msg.append ("Kid " + kid[1] +  " doesn't match mother")
        if self.mother[2] != self.mgrandfather[1]:
            error_msg.append ("Mother doesn't match maternal grandfather")
            passed_test = False
        if self.mother[3] != self.mgrandmother[1]:
            error_msg.append ("Mother doesn't match maternal grandmother")
            passed_test = False
        if self.father[2] != self.pgrandfather[1]:
            error_msg.append ("Father doesn't match maternal grandfather")
            passed_test = False
        if self.father[3] != self.pgrandmother[1]:
            error_msg.append ("Father doesn't match maternal grandmother")
            passed_test = False

        if not passed_test:
            error_msg.append ("Family " + self.family_id + " failed test")
            return error_msg
        return 1


    def __str__(self):
        if self.test() != 1:
            return "\nFailed test\n"
        ret_str = "Paternal Grandfather: " + self.pgrandfather[1] + "\n" + \
            "Paternal Grandmother: " + self.pgrandmother[1] + "\n" + \
            "Maternal Grandfather: " + self.mgrandfather[1] + "\n" + \
            "Maternal Grandmother: " + self.mgrandmother[1] + "\n" + \
            "Father: " + self.father[1] + "\n" + \
            "Mother: " + self.mother[1] + "\n" + \
            "\n".join([("Child: " + kid[1]) for kid in self.kids])
        return ret_str

    def pretty_print(self):
        if self.test() != 1:
            return "\nFailed test\n"
        overall_length = 8*len(self.kids)-4

        top_start = self.pgrandfather[1] + '    ' + self.pgrandmother[1]
        top_end = self.mgrandfather[1] + '    ' + self.mgrandmother[1]
        gap_len = max(overall_length - (len(top_start) + len(top_end)),4)
        top_str = top_start + " "*gap_len + top_end

        spacer1 = "\n  |______|" + " "*(gap_len+4) + "|______|  "
        spacer1 += "\n      |" + " "*(gap_len+10) + "|      \n"

        middle_start = "    " + self.father[1]
        middle_end = self.mother[1] + "    "
        gap_len = overall_length - (len(middle_start) + len(middle_end))
        middle_str = middle_start + " "*gap_len + middle_end
        
        spacer2 = "\n  ____|" + "_"*(gap_len+2) + "|____  \n  "
        for i in range(len(self.kids)):
            if i%2==0:
                spacer2 += "|       " 
            else:
                spacer2 += "|       "
        spacer2 += "\n"

        
        return "Pedigree " + self.family_id + ":\n" \
                + top_str + spacer1 \
                + middle_str + spacer2 \
                + "    ".join([kid[1] for kid in self.kids]) + "\n"

    def has_sample(self, sample_id):
        try:
            sample_id = str(sample_id)
        except:
            return "Invalid sample ID"

        for sample in self.all_members:
            if sample_id == sample[1]:
                return True
        return False

    def check_generation(self, sample_id):
        try:
            sample_id = str(sample_id)
        except:
            return "Invalid sample ID"
        if not self.has_sample(sample_id):
            return False
        for sample in self.kids:
            if sample_id == sample[1]:
                return "F2"
        if sample_id == self.father[1] or sample_id == self.mother[1]:
            return "F1"
        else:
            return "P0"


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
        family = Family(pedigrees[pedigree], pedigree)
#        if remove_failed and family.test() != 1:
#            continue
        families[pedigree] = family
    return families

def CreateListFromPED(ped, remove_failed = False):
    pedigrees = {}
    with open(ped, 'r') as ped_file:
        for line in ped_file:
            fields = line.strip().split()
            pedigree = fields[0]
            if pedigree not in pedigrees:
                pedigrees[fields[0]] = []
            pedigrees[pedigree].append(fields)

    families = []
    for pedigree in pedigrees:
        family = Family(pedigrees[pedigree], pedigree)
        if remove_failed and family.test() != 1:
            continue
        families.append(family)
    return families
