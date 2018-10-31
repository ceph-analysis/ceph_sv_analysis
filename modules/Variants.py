#! /usr/bin/env python
from __future__ import print_function
import sys

REF = [0,True,-1]
REF = [0]
ALT = [1,2]
def is_ref(gt):
    if gt in REF:
        return True
    if gt in ALT:
        return False
    else:
        return False

def is_homref(gt):
    #if the genotype is ref/ref, it's homref
    #if either is missing, assume it's ref
    if is_ref(gt[0]) and is_ref(gt[1]):
        return True
    return False

def is_het(gt):
    if gt[0] != gt[1]:
        return True
    return False

def sample_get_var_support(sample_idx, variant):
    AO = variant.format("AO")
    if AO is None:
        AO_len = 0
    else:
        AO_len = len(AO)

    if AO_len > 0 and sample_idx in AO:
        AO_support = variant.format('AO')[sample_idx]
        return int(AO_support)
    else:
        return 0

def sample_has_variant(sample_id, variant, vcf_samples, cutoff):
    sample_id = int(sample_id)
    for i in range(len(variant.genotypes)):
        gt_sample_id = int(vcf_samples[i])
        if gt_sample_id == sample_id:
            if (not is_homref(variant.genotypes[i])) | sample_get_var_support(i, variant) > cutoff:
                return True
            return False

def getVariantSupport(sample_id, variant,vcf_samples):
    sample_id = int(sample_id)
    for i in range(len(variant.genotypes)):
        gt_sample_id = int(vcf_samples[i])
        if gt_sample_id == sample_id:
            return sample_get_var_support(i, variant)

def is_common(variant, vcf_samples):
    variant_ped_ids = set()
    for i in range(len(variant.genotypes)):
        if (not is_homref(variant.genotypes[i])):
            ped_id = pedigree_lookup[vcf_samples[i]]
            variant_ped_ids.add(ped_id)
    if len(variant_ped_ids) > 1:
        return True
    return False
