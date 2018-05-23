#! /usr/bin/env python
import sys
import argparse
import cyvcf2
import os

ceph_vcf = "sv.pruned.re_svtyped.reclassed.filtered.gt50bp.allchr.CEPH.no_sname.vcf.gz"
ceph_ped = "16-08-06_WashU-Yandell-CEPH.ped"
sexes = {
    '1' : 'male',
    '2' : 'female',
    '3' : 'both'
}
# check this 
REF = 0
generations = ['P0', 'F1', 'F2']
exclude_pedigrees = []
pedigree_lookup = {}

parser = argparse.ArgumentParser(description="one-off script to find de novo SVs in the CEPH dataset")
parser.add_argument("-v", "--vcf", help="vcf file of called SVs for the whole sample set",default=ceph_vcf)
parser.add_argument("-p", "--ped", help="ped file CEPH families",default=ceph_ped)
parser.add_argument("-f", '--family', help="report only de novos for family specified (by pedigree ID)")
args = parser.parse_args()
pedigrees = {}
with open (args.ped, 'r') as pedfile:
    for line in pedfile:
        fields = line.strip().split()
        sample = {
            'pedigree_id'   : int(fields[0]),
            'id'            : int(fields[1]),
            'dad_id'        : int(fields[2]),
            'mom_id'        : int(fields[3]),
            'sample_sex'    : sexes[fields[4]],
            'offspring'     : [],
            'generation'    : ''
        }
        pedigree_lookup[sample['id']] = sample['pedigree_id']
        sample_id = sample['id']
        ped_id = sample['pedigree_id']
        if ped_id not in pedigrees:
            pedigrees[ped_id] = {}
        pedigrees[ped_id][sample_id] = sample
#recursively set the generation for a sample's offspring, starting from P0
def set_generations(pedigree,sample_id, generation):
    # if the generations are already set for this branch, move on
    if pedigree[sample_id]['generation'] != '':
        return pedigree

    #set generation for the current sample
    pedigree[sample_id]['generation'] = generations[generation]
    generation += 1
    
    #set generations for the kids
    if len(pedigree[sample_id]['offspring']) > 0:
        for offspring_id in pedigree[sample_id]['offspring']:
            pedigree = set_generations(pedigree, offspring_id, generation)
    return pedigree

# graph the relatedness (bidirectionally)  
for ped_id in pedigrees:
    for id1 in pedigrees[ped_id]:
        for id2 in pedigrees[ped_id]:
            sample2 = pedigrees[ped_id][id2]
            
            # if the first id is that of a parent, add the second id to the parent's offspring list
            if sample2['dad_id'] == id1 or sample2['mom_id'] == id1:
                pedigrees[ped_id][id1]['offspring'].append(id2)
    
    # add generation info 
    for id1 in pedigrees[ped_id]:
        # P0s should have children and not parents
        if (len(pedigrees[ped_id][id1]['offspring']) > 0) and \
                (pedigrees[ped_id][id1]['mom_id'] == 0) and \
                (pedigrees[ped_id][id1]['dad_id'] == 0):
            pedigrees[ped_id] = set_generations(pedigrees[ped_id], id1, 0)

    # look for f2s in each pedigree and print exclude the pedigree if there are none
    found = False
    for sample_id in pedigrees[ped_id]:
        sample = pedigrees[ped_id][sample_id]
        if "F2" == sample['generation']:
            found = True
    if not found:
        exclude_pedigrees.append(ped_id)

variant_info = {}
vcf = cyvcf2.VCF(os.path.expanduser(args.vcf))
vcf_samples = [int(x.replace("H_VT-", '')) for x in vcf.samples]

def is_homref(gt):
    if gt[0] == REF and gt[1] == REF:
        return True
    return False

def is_het(gt):
    if gt[0] != gt[1]:
        return True
    return False

def sample_has_variant(sample_id, variant):
    for i in range(len(variant.genotypes)):
        gt_sample_id = vcf_samples[i]
        if gt_sample_id == sample_id:
            if not is_homref(variant.genotypes[i]):
                return True

def parents_have_variant(sample_id, variant):
    return (sample_has_variant(sample['dad_id'], variant)) or (sample_has_variant(sample['mom_id'], variant))

def offspring_have_variant(sample_id, variant):
    sample = pedigrees[pedigree_lookup[sample_id]][sample_id]

    # check all offspring
    for offspring_id in sample['offspring']:
        if sample_has_variant(offspring_id, variant):
            return True
    return False

def is_common(variant):
    variant_ped_ids = set()
    for i in range(len(variant.genotypes)):
        if (not is_homref(variant.genotypes[i])):
            ped_id = pedigree_lookup[vcf_samples[i]]
            variant_ped_ids.add(ped_id)
    if len(variant_ped_ids) > 1:
        return True
    return False

denovo_count = 0
for variant in vcf:
    if variant.is_sv and not (is_common(variant)):
        for i in range(len(variant.genotypes)):
            gt = variant.genotypes[i]
            if is_het(gt):
                # this sample is het for this variant, so we add it to the variant info
                # if the sample has parents, consider them before adding
                sample_id = vcf_samples[i]
                if sample_id not in variant_info:
                    variant_info[sample_id] = []
                sample = pedigrees[pedigree_lookup[sample_id]][sample_id]
                
                # only report scores for one family if specified
                if args.family and (sample['pedigree_id'] != int(args.family)):
                    continue

                # if both of the sample's parents are known (an F1 missing parents is no use)
                if sample['dad_id'] == 0 or sample['mom_id'] == 0:
                    continue
                # if the sample has fewer than 5 offspring
                if len(sample['offspring']) < 5:
                    continue

                # if the total number of supporting reads for most common ALT is less than 5
                su = variant.format('AP')
                for i in range (len(variant.format('AS'))):
                    if len(su) > i:
                        su[i] += variant.format('AS')[i]
                    else:
                        su.append(variant.format('AS')[i])
                    
                if max(su) < 5:
                    continue

                # if the allele balance is less than .25 for ALT
#                if variant.format('AB')[i] < 0.1:
#                    continue

                #if the sample is a P0, add the variant
                if sample['generation'] == generations[0]:
                    variant_info[sample_id].append(variant)
                # else if the sample is an F1 and the variant isn't in the parents, but is in the F2, add it
                elif sample['generation'] == generations[1]:
                    
                            if not parents_have_variant(sample_id, variant) and offspring_have_variant(sample_id, variant):
                                variant_info[sample_id].append(variant)
                                denovo_count += 1
                else:
                    if parents_have_variant(sample_id, variant):
                        variant_info[sample_id].append(variant)
print (denovo_count)
