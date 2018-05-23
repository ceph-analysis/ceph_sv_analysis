#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse
import cyvcf2
import os
import families
import variants

ceph_vcf = "sv.pruned.re_svtyped.reclassed.filtered.gt50bp.allchr.CEPH.no_sname.vcf.gz"
ceph_ped = "16-08-06_WashU-Yandell-CEPH.ped"
sexes = {
    '1' : 'male',
    '2' : 'female',
    '3' : 'both'
}
# check this 
generations = ['P0', 'F1', 'F2']
exclude_pedigrees = []
pedigree_lookup = {}

parser = argparse.ArgumentParser(description="script to find de novo SVs in the CEPH dataset")
parser.add_argument("-v", "--vcf", help="vcf file of called SVs for the whole sample set",default=ceph_vcf)
parser.add_argument("-p", "--ped", help="ped file CEPH families",default=ceph_ped)
parser.add_argument("-f", '--family', help="report only de novos for family specified (by pedigree ID)")
args = parser.parse_args()

families = families.CreateDictFromPED(args.ped, True)

variant_info = {}
vcf = cyvcf2.VCF(os.path.expanduser(args.vcf))
vcf_samples = [int(x.replace("H_VT-", '')) for x in vcf.samples]


#step 1: make sure variants don't exist in multiple families (more than 1 for now)
filtered1_variants_by_ped = {}
variant_count = 0
for family_id in families:
    filtered1_variants_by_ped[family_id] = []
number = 0
for variant in vcf:
    affected_counts = []
#    number += 1
#    if number % 10 != 0:
#        continue
    affected_peds = set()
    for sample in vcf_samples:
        if variants.sample_has_variant(sample, variant, vcf_samples):
            for family_id in families:
                if families[family_id].has_sample(sample):
                    affected_peds.add(family_id)
    affected_counts.append(len(affected_peds))
    #if ((len(affected_peds)) == 2) or ((len(affected_peds)) == 1):
    if (len(affected_peds) == 1):
        for affected_ped in affected_peds:
            filtered1_variants_by_ped[affected_ped].append(variant)
            variant_count += 1
variant_count = 0
#step 2: make sure variants follow basic de novo pattern (not in P0, in F1, in F2)
#filtered2_variants_by_ped = {}
#filtered2_variants_by_ped[ped_id] = []
#the info needed for a samplot of the variant will be:
#CHROM POS END TYPE PO-0 P0-1 F1 F2
#with potentially a list of F2s
filtered_variants2 = []

for ped_id in filtered1_variants_by_ped:
    for variant in filtered1_variants_by_ped[ped_id]:
        affected_gens = {
            "P0" : [],
            "F1" : [],
            "F2" : []
        }
        for sample_list in families[ped_id].all_members:
            sample = int(sample_list[1])
            #malformed sample IDs (I hope)
            if variants.sample_has_variant(sample, variant, vcf_samples):
                gen = families[ped_id].check_generation(sample)
                affected_gens[gen].append(sample)
#                if gen == "P0":dd
#                    #a grandparent has this variant, so it's not de novo
#                    in_P0 = False
#                elif gen == "F1":
#                    #a parent has this variant, so it could be de novo in the F1
#                    in_F1 = True
#                elif gen == "F2":
#                    in_F2 = True
#                    #a child has this variant, so it could be de novo in the F1 

        # if there are no affected grandparents, exactly one affected parent, and at least 1 affected child, it could be de novo
        if (len(affected_gens["P0"]) == 0) and \
                (len(affected_gens["F1"]) == 1) and \
                (len(affected_gens["F2"]) > 0):
            F1 = affected_gens["F1"][0]
            #get gramma and grampaA
            f1_ismom = (int(families[ped_id].mother[1]) == F1)
            if f1_ismom:
                P0s = families[ped_id].mgrandfather[1],families[ped_id].mgrandmother[1]
            else:
                P0s = families[ped_id].pgrandfather[1],families[ped_id].pgrandmother[1]

            variant_entry = [
                variant.CHROM,
                variant.POS,
            ]
            if variant.INFO['SVTYPE'] == "BND":
                variant_entry.append("NA")
            else:
                variant_entry.append(variant.INFO['END'])
            
            variant_entry.append(variant.INFO['SVTYPE'])
            variant_entry.append(P0s[0])
            variant_entry.append(P0s[1])
            variant_entry.append(F1)
            for F2 in affected_gens['F2']:
                variant_entry.append(F2)
            filtered_variants2.append(variant_entry)
            #filtered2_variants_by_ped[ped_id].append(variant)
            variant_count += 1
for variant in filtered_variants2:
    print("\t".join([str(x) for x in variant]))
