#! /usr/bin/env python
from __future__ import print_function
import sys

REF = 0
def is_homref(gt):
    if gt[0] == REF and gt[1] == REF:
        return True
    return False

def is_het(gt):
    if gt[0] != gt[1]:
        return True
    return False

def sample_has_variant(sample_id, variant, vcf_samples):
    for i in range(len(variant.genotypes)):
        gt_sample_id = vcf_samples[i]
        if gt_sample_id == sample_id:
            if not is_homref(variant.genotypes[i]):
                return True

def is_common(variant, vcf_samples):
    variant_ped_ids = set()
    for i in range(len(variant.genotypes)):
        if (not is_homref(variant.genotypes[i])):
            ped_id = pedigree_lookup[vcf_samples[i]]
            variant_ped_ids.add(ped_id)
    if len(variant_ped_ids) > 1:
        return True
    return False
