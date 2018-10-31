#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse
import cyvcf2
import os
import Families
import Variants_smv as Variants 
import argparse

ceph_vcf = "/uufs/chpc.utah.edu/common/home/u1072557/ceph_denovosv/sv.pruned.re_svtyped.reclassed.filtered.gt50bp.allchr.CEPH.no_sname.vcf.gz"
ceph_vcf = "/uufs/chpc.utah.edu/common/home/u1072557/ceph_denovosv/ceph-smoove-0.1.11.smoove.square.vcf.gz"
ceph_ped = "/uufs/chpc.utah.edu/common/home/u1072557/ceph_denovosv/16-08-06_WashU-Yandell-CEPH.ped"
vcfs = ["/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1330/1330-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1331/1331-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1334/1334-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1340/1340-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1341/1341-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1344/1344-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1345/1345-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1346/1346-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1347/1347-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1349/1349-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1350/1350-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1353/1353-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1354/1354-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1355/1355-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1356/1356-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1357/1357-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1358/1358-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1362/1362-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1377/1377-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1408/1408-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1418/1418-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1420/1420-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1424/1424-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1444/1444-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1447/1447-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1451/1451-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1454/1454-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1458/1458-smoove.genotyped.vcf.gz",
#    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1459/1459-smoove.genotyped.vcf.gz",
    "/uufs/chpc.utah.edu/common/home/u6000771/Projects/2018/smoove-CEPH/results-dn-by-fam/1463/1463-smoove.genotyped.vcf.gz"
]
vcfs = [ceph_vcf]

parser = argparse.ArgumentParser(description="script to find de novo SVs in the CEPH dataset")
#parser.add_argument("-v", "--vcf", help="vcf file of called SVs for the whole sample set",default=ceph_vcf)
parser.add_argument("-p", "--ped", help="ped file CEPH families",default=ceph_ped)
parser.add_argument("-f", '--family', help="report only de novos for family specified (by pedigree ID)")
parser.add_argument("-c", '--support_cutoff', help="support cutoff (AO). If fewer reads support a call and GT is homref, no variant. ", type=int, default=0)
args = parser.parse_args()


sexes = {
    '1' : 'male',
    '2' : 'female',
    '3' : 'both'
}
# check this 
generations = ['P0', 'F1', 'F2']
families = Families.CreateFamilies(args.ped, True)
#counts = {'P0':0,'F1':0,'F2':0}
#valid_f1_fam_count = 0
#valid_f2_fam_count = 0
#for family_id in families:
##    print (families[family_id].pretty_print())
#    f1s = families[family_id].getValidF1Count()
#    counts['F2'] += families[family_id].getValidF2Count()
#    counts['P0'] += families[family_id].getValidP0Count()
#    if f1s == 2:
#        valid_f1_fam_count += 1
#    if f1s > 0:
#        valid_f2_fam_count += 1
#    counts['F1'] += f1s
#print(counts)
#print(valid_f1_fam_count)
#print(valid_f2_fam_count)
#sys.exit()
variant_sample_states = {}
#step 1: get variants with sample data (id, support, ped, sex, generation) 
filtered_vars = {}
for family_id in families:
    filtered_vars[family_id] = []
vcf_handles = {}
for vcf_name in vcfs:
    vcf = cyvcf2.VCF(os.path.expanduser(vcf_name))
    #vcf handles will contain the name of the vcf as key, and as a pair handle and the samples in it 
    vcf_handles[vcf_name] = vcf
    for variant in vcf:
        if variant.INFO['SVTYPE'] == "BND" or variant.INFO['SVTYPE'] == "MEI":
            continue
        chrom,pos,end = variant.CHROM,variant.POS,variant.INFO['END']


        #if not (chrom == "19" and str(pos) == "38459850" and str(end) == "38460147"):
            #continue
        variant_name = "-".join([str(chrom),str(pos),str(end)])
        if variant_name not in variant_sample_states:
            variant_sample_states[variant_name] = {}
        for i in range(len(vcf.samples)):
            sample = vcf.samples[i]
            family_id = 0
            for fam_id in families:
                if families[fam_id].hasSample(sample):
                    family_id = fam_id
                    break
            if family_id == 0:
                #print("Sample not found. ID: " + sample)
                continue

            AO_field = variant.format("AO")
            AO = 0
            if AO_field is not None and len(AO_field) > i:
                AO = variant.format("AO")[i][0]

            variant_sample_states[variant_name][sample] = {
                "AO": int(AO),
                "chrom": str(chrom),
                "pos":pos, 
                "end":end,
                "type":variant.INFO["SVTYPE"],
                "sex":families[family_id].checkSex(sample),
                "generation":families[family_id].checkGeneration(sample),
                "fam_id": fam_id,
                "vcf": vcf_name
            }

AO_cutoff = 0
#step 1.5: remove variants that appear in multiple families
variants_to_remove = []
for variant in variant_sample_states:
    affected_peds = []
    for sample in variant_sample_states[variant]:
        if variant_sample_states[variant][sample]['AO'] > AO_cutoff:
            affected_peds.append(variant)
    if len(list(set(affected_peds))) > 1:
        variants_to_remove.append(variant)
for variant in variants_to_remove:
    variant_sample_states.pop(variant)

#step 2 for each variant, make sure that it only exists in either f1+f2s or f2s only
#add to appropriate list
f1_denovos = {}
f2_denovos = {}
f2_mosaics = {}
for variant in variant_sample_states:
    affected_gens = []
    for sample in variant_sample_states[variant]:
        if variant_sample_states[variant][sample]["AO"] > AO_cutoff:
            affected_gens.append(variant_sample_states[variant][sample]["generation"])
    if "P0" in affected_gens:
        continue
    elif affected_gens.count("F1") > 1:
        continue
    elif "F1" not in affected_gens:
        #need to make sure the F2s' parents are both available
        for sample in variant_sample_states[variant]:
            if variant_sample_states[variant][sample]["AO"] <= AO_cutoff:
                continue
            family_id = variant_sample_states[variant][sample]["fam_id"]
            dad = families[family_id].checkDad(sample)
            mom = families[family_id].checkMom(sample)
            if str(dad) != "0" and str(mom) != "0":
                if affected_gens.count("F2") > 1:
                    if variant not in f2_mosaics:
                        f2_mosaics[variant] = {}
                    f2_mosaics[variant][sample] = variant_sample_states[variant][sample]
                else:
                    if variant not in f2_denovos:
                        f2_denovos[variant] = {}
                    f2_denovos[variant][sample] = variant_sample_states[variant][sample]
    elif "F1" in affected_gens and "F2" in affected_gens:
        #if we get to this point, there are no grandparents with the variant and only one parent with the variant
        #but at least one child
        #still need to know if the parent is the right one
        f1 = ''
        family_id = 0
        for sample in variant_sample_states[variant]:
            if str(sample) != "0" and variant_sample_states[variant][sample]["generation"] == "F1" and variant_sample_states[variant][sample]["AO"] > AO_cutoff:
                family_id = variant_sample_states[variant][sample]["fam_id"]
                dad_id = families[family_id].checkDad(sample)
                mom_id = families[family_id].checkMom(sample)
                if str(dad_id) != '0' and str(mom_id) != '0':
                    f1 = sample
                    if variant not in f1_denovos:
                        f1_denovos[variant] = {}
                    f1_denovos[variant][sample] = variant_sample_states[variant][sample] 
                    break
        if f1 != "":
            for sample in variant_sample_states[variant]:
                if variant_sample_states[variant][sample]["AO"] <= AO_cutoff:
                    continue
                if variant_sample_states[variant][sample]["generation"] == "F2":
                    #if the f1 is the parent of the f2 who has the variant, it's a valid possible
                    #as long as the parent of the f1 has a valid ID
                    if (families[family_id].checkDad(sample) == f1) or (families[family_id].checkMom(sample) == f1):
                        f1_denovos[variant][sample] = variant_sample_states[variant][sample] 
                    #gotta remove the break once things are narrowed down
print ("f1: " + str(len(f1_denovos)))
print ("f2: " + str(len(f2_denovos)))
print ("f2 mosaics: " + str(len(f2_mosaics)))

def get_AO(sample_id, chrom,start,end, vcf_name):
    vcf = vcf_handles[vcf_name]
    sample_idx = 0
    for i in range(len(vcf.samples)):
        if vcf.samples[i] == sample_id:
            sample_idx = i
            break
    for variant in vcf:
        if variant.chrom == chrom and variant.pos == start and variant.INFO['END'] == end:
            return variant.format("AO")[sample_idx][0]
    return 0


for variant in f1_denovos:
    variant_entry = []
    samples_to_add ={ 
        "P0": [],
        "F1": [],
        "F2": []
    }
    for sample in f1_denovos[variant]:
        sample_info = f1_denovos[variant][sample]
        if len(variant_entry) == 0:
            variant_entry = [
                sample_info['chrom'], 
                sample_info['pos'],
                sample_info['end'],
                sample_info['type']
            ]
        samples_to_add[sample_info['generation']].append(",".join([sample, sample_info['generation'], sample_info['sex'], str(sample_info['AO'])]))
        if sample_info['generation'] == "P0":
            print ("sometins rong: " + sample_info)
            sys.exit()
        elif sample_info['generation'] == "F1":
            dad_id = families[sample_info['fam_id']].checkDad(sample)
            dad_AO = get_AO(dad_id,sample_info['chrom'],sample_info['pos'],sample_info['end'],sample_info['vcf'])
            mom_id = families[sample_info['fam_id']].checkMom(sample)
            mom_AO = get_AO(mom_id,sample_info['chrom'],sample_info['pos'],sample_info['end'],sample_info['vcf'])
            
            samples_to_add['P0'].append(",".join([dad_id, "P0", '1', str(dad_AO)]))
            samples_to_add['P0'].append(",".join([mom_id, "P0", '2', str(mom_AO)]))
    variant_entry = variant_entry + samples_to_add['P0'] + samples_to_add['F1'] + samples_to_add['F2']
    print ("\t".join([str(x) for x in variant_entry]) )

for variant in f2_denovos:
    variant_entry = []
    samples_to_add ={ 
        "F1": [],
        "F2": []
    }
    for sample in f2_denovos[variant]:
        sample_info = f2_denovos[variant][sample]
        if len(variant_entry) == 0:
            variant_entry = [
                sample_info['chrom'], 
                sample_info['pos'],
                sample_info['end'],
                sample_info['type']
            ]
        samples_to_add[sample_info['generation']].append(",".join([sample, sample_info['generation'], sample_info['sex'], str(sample_info['AO'])]))
        if sample_info['generation'] == "P0" or sample_info['generation'] == "F1":
            print ("sometins rong: " + sample_info)
            sys.exit()
            
        dad_id = families[sample_info['fam_id']].checkDad(sample)
        mom_id = families[sample_info['fam_id']].checkMom(sample)
        if str(dad_id) == "0" or str(mom_id) == "0":
            continue
        dad_AO = get_AO(dad_id,sample_info['chrom'],sample_info['pos'],sample_info['end'],sample_info['vcf'])
        mom_AO = get_AO(mom_id,sample_info['chrom'],sample_info['pos'],sample_info['end'],sample_info['vcf'])
        
        samples_to_add['F1'].append(",".join([dad_id, "F1", '1', str(dad_AO)]))
        samples_to_add['F1'].append(",".join([mom_id, "F1", '2', str(mom_AO)]))
    variant_entry = variant_entry + samples_to_add['F1'] + samples_to_add['F2']
    print ("\t".join([str(x) for x in variant_entry]) )

for variant in f2_mosaics:
    variant_entry = []
    samples_to_add ={ 
        "F1": [],
        "F2": []
    }
    for sample in f2_mosaics[variant]:
        sample_info = f2_mosaics[variant][sample]
        if len(variant_entry) == 0:
            variant_entry = [
                sample_info['chrom'], 
                sample_info['pos'],
                sample_info['end'],
                sample_info['type']
            ]
        samples_to_add[sample_info['generation']].append(",".join([sample, sample_info['generation'], sample_info['sex'], str(sample_info['AO'])]))
        if sample_info['generation'] == "P0" or sample_info['generation'] == "F1":
            print ("sometins rong: " + sample_info)
            sys.exit()
            
        dad_id = families[sample_info['fam_id']].checkDad(sample)
        mom_id = families[sample_info['fam_id']].checkMom(sample)
        if str(dad_id) == "0" or str(mom_id) == "0":
            continue
        dad_AO = get_AO(dad_id,sample_info['chrom'],sample_info['pos'],sample_info['end'],sample_info['vcf'])
        mom_AO = get_AO(mom_id,sample_info['chrom'],sample_info['pos'],sample_info['end'],sample_info['vcf'])
        
        samples_to_add['F1'].append(",".join([dad_id, "F1", '1', str(dad_AO)]))
        samples_to_add['F1'].append(",".join([mom_id, "F1", '2', str(mom_AO)]))
    variant_entry = variant_entry + samples_to_add['F1'] + samples_to_add['F2']
    print ("\t".join([str(x) for x in variant_entry]) )




     










#
#variant_count = 0
##step 2: make sure variants follow basic de novo pattern (not in P0, in F1, in F2)
##filtered2_variants_by_ped = {}
##filtered2_variants_by_ped[ped_id] = []
##the info needed for a samplot of the variant will be:
##CHROM POS END TYPE PO-0 P0-1 F1 F2
##with potentially a list of F2s
#filtered_variants_f1 = []
#filtered_variants_f2mislabeled = []
#filtered_variants_f2_germline = []
#filtered_variants_somatic_mosaic = []
#
#for ped_id in filtered_vars:
#    for variant in filtered_vars[ped_id]:
#        sample_variant_support = {}
#        affected_gens = {
#            "P0" : [],
#            "F1" : [],
#            "F2" : []
#        }
#        for sample in families[ped_id].sample_ids:
#            if Variants.sample_has_variant(sample, variant, vcf_samples, args.support_cutoff):
#                gen = families[ped_id].checkGeneration(sample)
#                affected_gens[gen].append(sample)
#                sample_variant_support[sample] = Variants.getVariantSupport(sample, variant,vcf_samples)
#
#        # if there are no affected grandparents, exactly one affected parent, and at least 1 affected child, it could be de novo f1
#        if (len(affected_gens["P0"]) == 0) and \
#                (len(affected_gens["F1"]) == 1) and \
#                (len(affected_gens["F2"]) > 0):
#            F1 = affected_gens["F1"][0]
#            #get gramma and grampa
#            try:
#                f1_ismom = (int(families[ped_id].f1_female) == int(F1))
#            except:
#                print ("error finding F1 " + F1 )
#                continue
#
#            if f1_ismom:
#                P0s = families[ped_id].p0_male_maternal,families[ped_id].p0_female_maternal
#            else:
#                P0s = families[ped_id].p0_male_paternal,families[ped_id].p0_female_paternal
#
#            #if we're missing a P0, assume the variant came from that sample
#            unknown_p0 = False
#            for P0 in P0s:
#                if P0 == "????":
#                    unknown_p0 = True
#            if unknown_p0:
#                continue
#
#            variant_entry = [
#                variant.CHROM,
#                variant.POS,
#            ]
#            if variant.INFO['SVTYPE'] == "BND":
#                variant_entry.append("NA")
#            else:
#                variant_entry.append(variant.INFO['END'])
#            
#            variant_entry.append(variant.INFO['SVTYPE'])
#            samples_to_add = list(P0s) +[F1] + [x for x in affected_gens["F2"]]
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_to_add)
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_variant_support[sample_to_add])
#
#            #variant_entry.append(P0s[0])
#            #variant_entry.append(P0s[1])
#            #variant_entry.append(F1)
#            #for F2 in affected_gens['F2']:
#            #    variant_entry.append(F2)
#            
#            filtered_variants_f1.append(variant_entry)
#
#            #trying weird stuff
#            variant_entry = [
#                variant.CHROM,
#                variant.POS,
#            ]
#            if variant.INFO['SVTYPE'] == "BND":
#                variant_entry.append("NA")
#            else:
#                variant_entry.append(variant.INFO['END'])
#            variant_entry.append(variant.INFO['SVTYPE'])
#            
#            samples_to_add = [x for x in affected_gens["F2"]] + [families[ped_id].f1_male, families[ped_id].f1_male]
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_to_add)
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_variant_support[sample_to_add])
#
#
#            #for F2 in affected_gens['F2']:
#            #    variant_entry.append(F2)
#            #variant_entry.append(families[ped_id].f1_male)
#            #variant_entry.append(families[ped_id].f1_female)
#            filtered_variants_f2mislabeled.append(variant_entry)
#
#
#            variant_count += 1
#        # if there are no affected grandparents, no affected parent, and exactly 1 affected child, it could be de novo germline f2
#        elif (len(affected_gens["P0"]) == 0) and \
#                (len(affected_gens["F1"]) == 0) and \
#                (len(affected_gens["F2"]) == 1):
#            variant_entry = [
#                variant.CHROM,
#                variant.POS,
#            ]
#            if variant.INFO['SVTYPE'] == "BND":
#                variant_entry.append("NA")
#            else:
#                variant_entry.append(variant.INFO['END'])
#            #if we're missing F1s, assume the variant came from that sample
#            if families[ped_id].f1_male == "????" or families[ped_id].f1_female == "????":
#                continue
#            variant_entry.append(variant.INFO['SVTYPE'])
#
#            samples_to_add = [affected_gens["F2"][0]] + [families[ped_id].f1_male, families[ped_id].f1_male]
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_to_add)
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_variant_support[sample_to_add])
#
#            #variant_entry.append(affected_gens["F2"][0])
#            #variant_entry.append(families[ped_id].f1_male)
#            #variant_entry.append(families[ped_id].f1_female)
#            filtered_variants_f2_germline.append(variant_entry)
#            variant_count += 1
#        # if there are no affected grandparents, no affected parent, and at least 2 affected children, it could be de novo somatic mosaic f1->f2
#        elif (len(affected_gens["P0"]) == 0) and \
#                (len(affected_gens["F1"]) == 0) and \
#                (len(affected_gens["F2"]) > 1):
#            variant_entry = [
#                variant.CHROM,
#                variant.POS,
#            ]
#            if variant.INFO['SVTYPE'] == "BND":
#                variant_entry.append("NA")
#            else:
#                variant_entry.append(variant.INFO['END'])
#            
#            #if we're missing F1s, assume the variant came from that sample
#            if families[ped_id].f1_male == "????" or families[ped_id].f1_female == "????":
#                continue
#            variant_entry.append(variant.INFO['SVTYPE'])
#
#            samples_to_add = [families[ped_id].f1_male, families[ped_id].f1_male] + [x for x in affected_gens["F2"]]
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_to_add)
#            for sample_to_add in samples_to_add:
#                variant_entry.append(sample_variant_support[sample_to_add])
#
#
#
#            #variant_entry.append(families[ped_id].f1_male)
#            #variant_entry.append(families[ped_id].f1_female)
#            #for F2 in affected_gens['F2']:
#            #    variant_entry.append(F2)
#            filtered_variants_somatic_mosaic.append(variant_entry)
#            variant_count += 1
#
#print("#F1 de novos: ")
#print ("# CHROM\tPOS\tEND\tTYPE\tP0-1\tP0-2\tF1\tF2+")
#for variant in filtered_variants_f1:
#    print("\t".join([str(x) for x in variant]))
#
#print("#F2 mislabeled de novos:")
#print ("# CHROM\tPOS\tEND\tTYPE\tF2\tF1-1\tF1-2")
#for variant in filtered_variants_f2mislabeled:
#    print("\t".join([str(x) for x in variant]))
#
#print("#F2 germline de novos:")
#print ("# CHROM\tPOS\tEND\tTYPE\tF2\tF1-1\tF1-2")
#for variant in filtered_variants_f2_germline:
#    print("\t".join([str(x) for x in variant]))
#
#print("#F2 somatic mosaics:")
#print ("# CHROM\tPOS\tEND\tTYPE\tF1-1\tF1-2\tF2")
#for variant in filtered_variants_somatic_mosaic:
#    print("\t".join([str(x) for x in variant]))
