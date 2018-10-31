from __future__ import print_function
import sys
import argparse
import os

parser = argparse.ArgumentParser(description="creates a bash script to run cnvnator for the CEPH samples")
parser.add_argument("-d", "--dir", help="directory of bams to call CNVs from",required=False, default="/scratch/ucgd/lustre/u1006375/ceph/ceph-bams/")
args = parser.parse_args()


root_dir = "/uufs/chpc.utah.edu/common/home/u1072557/ceph_denovosv/cnvnator_calling/cnvnator_results"
outfile_dir = "/uufs/chpc.utah.edu/common/home/u1072557/ceph_denovosv/cnvnator_calling/running_files/"
for sample in os.listdir(args.dir):
    sample_name,ext = os.path.splitext(sample)
    if ext != ".bam":
        continue
    sample_root_name = os.path.join(root_dir,sample_name +".root")
    sample_running_name = os.path.join(outfile_dir, sample_name+".sh")
    with open(sample_running_name, 'w') as outfile:
        outfile.write("#!/usr/bin/bash\n"+ "set -eo\n")
        outfile.write('export ROOTSYS="/uufs/chpc.utah.edu/common/home/u0055382/root"\n')
        outfile.write('export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ROOTSYS/lib"\n')
        outfile.write('export cnvnator="/uufs/chpc.utah.edu/common/home/u0055382/CNVnator_v0.3.3/src/cnvnator"\n')


        command = "$cnvnator hg19 -root "+sample_root_name +" -tree " + os.path.join(args.dir, sample) + "\n"
        command += "$cnvnator hg19 -root "+sample_root_name+" -his 100 -d /uufs/chpc.utah.edu/common/home/u1072557/na12878_cnvnator/reference"+ "\n"
        command += "$cnvnator -root "+sample_root_name+" -stat 100"+ "\n"
        command += "$cnvnator -root "+sample_root_name+" -partition 100"+ "\n"
        sample_calls_name = sample_root_name.replace("root", "calls")
        
        command += "$cnvnator -root "+sample_root_name+" -call 100 > " + sample_calls_name + "\n"

        sample_bvcf_name = sample_root_name.replace("root", "vcf.gz")

        command += "/uufs/chpc.utah.edu/common/home/u0055382/CNVnator_v0.3.3/cnvnator2VCF.pl " + sample_calls_name + " | bgzip > " + sample_bvcf_name+ "\n"
        command += "echo " + sample_name
        outfile.write (command)
        print ("bash " + sample_running_name)
