#! /usr/bin/env python
from __future__ import print_function
import sys
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', RuntimeWarning)
    import pandas as pd
sys.path.append('/uufs/chpc.utah.edu/common/home/u1072557/ceph_denovosv/ceph_scripts/modules/')
import argparse
import cyvcf2
import os
import Families
import Variants as Variants 
import argparse
import numpy as np
import pandas as pd


#########################################################################################################
# define functions
#########################################################################################################

#step 1: get variants with sample data (id, support, ped, sex, generation) 
def vcf_to_matrix(vcf_name):
    vcf = cyvcf2.VCF(os.path.expanduser(vcf_name))
    variant_support_dict = {}

    for variant in vcf:
        #we can't do much with either of these types
        if variant.INFO['SVTYPE'] == "BND" or variant.INFO['SVTYPE'] == "MEI":
            continue
        # get essential information about variant
        chrom,pos,end,svtype = variant.CHROM,variant.POS,variant.INFO['END'], variant.INFO['SVTYPE']
        variant_name = "_".join([str(chrom),str(pos),str(end),svtype])
        
        support_array = []
        # need to get the support from either AO, PR+SR, or genotype
        if 'AO' in vcf:    #Smoove
            support_array = np.array([x[0] for x in variant.format("AO")])
            strong_support = 4
            weak_support = 0
        elif 'PR' in vcf or 'SR' in vcf:  #manta
            split_read_su = np.zeros(len(vcf.samples))
            if 'SR' in vcf and variant.format("SR") is not None:
                split_read_su = np.array([x[0] for x in variant.format("SR")])

            paired_end_su = np.zeros(len(vcf.samples))
            if 'PR' in vcf and variant.format("PR") is not None:
                split_read_su = np.array([x[0] for x in variant.format("PR")])
            support_array = np.add(split_read_su, paired_end_su).astype(int)
            strong_support = 4
            weak_support = 0
        elif variant.genotypes:  #cnvnator
            support_array = [1 in x[:2] for x in variant.genotypes]
            strong_support = True
            weak_support = True
        else:
            print("No valid support field in VCF (checked for AO, PR, SR, GT)", file=sys.stderr)
            sys.exit()
        variant_support_dict[variant_name] = support_array
    variant_support = pd.DataFrame.from_dict(variant_support_dict, orient='index', columns=vcf.samples)
    return variant_support, strong_support, weak_support

# returns only the variants where at least one of the given samples has at least greater than min_support support
def extract_variants_by_samples(variant_support, sample_ids, min_support):
    if type(min_support) == int:
        extracted_variant_support = variant_support.loc[(variant_support[sample_ids] > min_support).any(axis="columns")]
    elif type(min_support) == bool:
        extracted_variant_support = variant_support.loc[(variant_support[sample_ids] == min_support).any(axis="columns")]
    return extracted_variant_support


# returns only the variants where none of the given samples has min_support support 
# basically, exclude if P0s have variant
def exclude_variants_by_samples(variant_support, sample_ids, min_support):
    #select and keep all regions where every one of the given samples has less than the min_support
    if type(min_support) == int:
        unexcluded_variant_support = variant_support.loc[(variant_support[sample_ids] <= min_support).all(axis="columns")]
    # select and keep regions where the given samples are false (ie, do not have variant)
    elif type(min_support) == bool:
        unexcluded_variant_support = variant_support.loc[(variant_support[sample_ids] != min_support).all(axis="columns")]
    return unexcluded_variant_support

def get_sample_generations(ped, family, verbose):
    families = Families.CreateFamilies(ped, True)
    generations = {'P0':[], 'F1':[], 'F2':[]}
    for family_id in families:
        if family in family_id:
            if verbose:
                print(families[family_id].pretty_print())
            for gen in generations:
                generations[gen] +=  [x[1] for x in families[family_id].generations[gen]]
    return generations

#########################################################################################################
# define args
#########################################################################################################
parser = argparse.ArgumentParser(description="script to find de novo SVs in a family of the CEPH dataset")
parser.add_argument("-v", 
    "--vcf", 
    help="vcf files of SVs"
)
parser.add_argument("-p", 
    "--ped", 
    help="ped file CEPH families"
)
parser.add_argument("-f", 
    '--family', 
    help="family specified (by pedigree ID)",
    required=True
)
parser.add_argument("-c", 
    '--support_cutoff', 
    help="support cutoff (AO). If fewer reads support a call and GT is homref, no variant.", 
    type=int, 
    default=0
)
parser.add_argument('--verbose', 
    help="print extra stuff out, like the ascii representation of the family", 
    action="store_true", 
    default=False
)
args = parser.parse_args()

#########################################################################################################
# main 
#########################################################################################################
# find F1 de novos that are transmitted to F2 generation
#########################################################################################################
# convert VCF inputs to dataframes with the support or genotype information for all samples
variant_support, strong_support, weak_support = vcf_to_matrix(args.vcf)

# find the variants that only appear in F1 and F2 samples
generations = get_sample_generations(args.ped, args.family, args.verbose)

# remove variants where P0s have even weak support
f1_support = exclude_variants_by_samples(variant_support, generations['P0'], weak_support)

# remove variants where F2s lack even weak support
f1_support = extract_variants_by_samples(f1_support, generations['F2'], weak_support)

# remove variants where F1s lack even weak  support
f1_support = extract_variants_by_samples(f1_support, generations['F1'], weak_support)

#remove variants where no F1 or F2 sample has strong support
f1_support = extract_variants_by_samples(f1_support, generations['F1'] + generations['F2'], strong_support)
print (f1_support)
print()

#########################################################################################################
# find F2 de novos
#########################################################################################################

# remove variants that have even weak support in the P0 or F1 generations
f2_support = exclude_variants_by_samples(variant_support, generations['P0']+generations['F1'], weak_support)

# remove variants that don't have even weak support in the F2 generation (this shouldn't be changing anything probably)
f2_support = extract_variants_by_samples(f2_support, generations['F2'], weak_support)

print (f2_support)
