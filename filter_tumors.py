"""
Second step of the Tuba-seq pipeline- takes as input the output of count_reads.py
This script performs 3 functions:
1) identication and removal of tumors that don't map to any project
2) identification and removal of tumors that map to other projects (i.e., likely contamination)
3) clustering of tumors with highlyg similar barcodes (small tumors with barcodes very similar to those of much larger tumors are collapsed into them)
"""


########################################################################################
# 1. Load packages
########################################################################################
import sys
import os
import regex as re
import gzip
import getopt
import time
import datetime
import numpy as np
from datetime import timedelta
from subprocess import call
from helper_functions import hamming_distance, check_master_sgid_file, read_parameter_info, make_sgIDDict

now = datetime.datetime.now()
start_time= time.monotonic()

########################################################################################
# 2. Take arguments and build necessary resources
########################################################################################
try:
    opts, args = getopt.getopt(sys.argv[1:], "s:", ["sample=", 
        "sample_path=", 
        "barcode_length=", 
        "spikein_barcode_lengths=", 
        "sgRNAs=", 
        "sgids=", 
        "bc_dist=", 
        "R1_regex_looseness=", 
        "R2_regex_looseness=", 
        "root=", 
        "project_name=", 
        "parameter_name="])
except getopt.GetoptError:
    print("no arguments recognized\n")
    sys.exit(2)

for opt, arg in opts:
    if opt in ("--sample"):
        sample = arg
    elif opt in ("--sample_path"):
        sample_path = arg
    elif opt in ("--barcode_length"):
        barcode_length = arg
    elif opt in ("--spikein_barcode_lengths"):
        spikein_barcode_lengths = arg.strip().split(",")
    elif opt in ("--sgRNAs"):
        sgRNAs = arg.strip().split(",")
    elif opt in ("--sgids"):
        sgids = arg.strip().split(",")
    elif opt in ("--bc_dist"):
        bc_dist = arg
    elif opt in ("--R1_regex_looseness"):
        R1_regex_looseness = arg
    elif opt in ("--R2_regex_looseness"):
        R2_regex_looseness = arg
    elif opt in ("--root"):
        root = arg
    elif opt in ("--project_name"):
        project_name = arg
    elif opt in ("--parameter_name"):
        parameter_name = arg
    else:
        assert False, "unhandled option"


def check_sgid_against_other_projects(sgids_other_projects, current_sgid):
    if current_sgid in sgids_other_projects.keys():
        flag = True
    else:
        flag = False
    return(flag)

sgids_other_projects = check_master_sgid_file("/labs/mwinslow/Emily/tubaseq_pipeline/egs_all_sgRNAs.txt", sgRNAs, project_name)

sgIDDict = make_sgIDDict(sgids, sgRNAs)

########################################################################################
# 3. Set up input/output Files
########################################################################################

naming_stub = project_name + "_" + parameter_name + "_" + sample
input_name = root + "/" + project_name + "/" + parameter_name + "/raw_counts/" + naming_stub + "_counts.txt"
record_name = root + "/" + project_name + "/" + parameter_name + "/records/" + naming_stub + "_record.txt"
no_sgid_name = root + "/" + project_name + "/" + parameter_name + "/rejected/" + naming_stub + "_sgid_mismatchs.txt"
output_name = root + project_name + "/" + parameter_name + "/filtered/" + naming_stub + "_clustered.txt"
record = open(record_name, 'a+')
output = open(output_name, 'wt')
_input = open(input_name, 'rt')
_input.readline()

wrong_sgid = open(no_sgid_name, 'wt')
wrong_sgid.write("Sample,SgID_seq,OtherProj,BC,RC\n")

output.write("sgid,bc,rc\n")


parameter_info_file = root + "/tubaseq_parameter_files/" + parameter_name + "_parameter_file.txt"
parameter_info = read_parameter_info(parameter_info_file)
threshold = float(parameter_info.read_error_size_threshold)

#########################################################################################
# 4. read in data
########################################################################################

tumors_in = 0
reads_in = 0
none_tumors = 0
none_tumor_reads = 0
tumors_match_other = 0
reads_match_other = 0
tumors_dont_match_other = 0
reads_dont_match_other = 0

tumors_out = 0
reads_out = 0


tumors_all = {}
tumors_all_in = {}
reads_all_in = {}

for l in _input:
    fields = l.strip().split(",")
    sgid = fields[0]

    bc = fields[1]
    sgid_seq = fields[3]
    rc = int(fields[2])
    tumors_in += 1
    reads_in += rc

    if sgid == "None":
        none_tumors +=1
        none_tumor_reads+=rc
        if check_sgid_against_other_projects(sgids_other_projects, sgid_seq) == True:
            tumors_match_other += 1
            reads_match_other += rc
            wrong_sgid.write("{},{},{},{},{}\n".format(sample,sgid_seq,"-".join(sgids_other_projects[sgid_seq]),bc,rc))
        else:
            tumors_dont_match_other += 1
            reads_dont_match_other += rc
            wrong_sgid.write("{},{},None,{},{}\n".format(sample,sgid_seq,bc,rc))
    else:

        if sgid in tumors_all.keys():
            tumors_all_in[sgid]+=1
            tumors_all[sgid][0].append(bc)
            tumors_all[sgid][1].append(rc)
            tumors_all[sgid][2].append(1)
            reads_all_in[sgid]+=rc
        else:
            tumors_all[sgid] = [[bc],[rc],[1]]
            tumors_all_in[sgid]=1
            reads_all_in[sgid] = rc

record.write("###########CLUSTERING NOW############\n")
record.write("Going into clustering, there are a total of {} tumors and {} reads\n".format(tumors_in, reads_in))
record.write("Tumors not mapping to project: {}\n".format(none_tumors))
record.write("Reads not mapping to this project: {}\n".format(none_tumor_reads))
record.write("Tumors mapping to other projects: {}\n".format(tumors_match_other))
record.write("Reads mapping to other projects: {}\n".format(reads_match_other))
record.write("Tumors NOT mapping to other projects: {}\n".format(tumors_dont_match_other))
record.write("Reads NOT mapping to other projects: {}\n".format(reads_dont_match_other))
record.write("Reads excluding none/other projects: {}\n".format((int(reads_in)-int(none_tumor_reads))))

 
collapsed_per_sgid = {}
collapsed_per_sgid_reads = {}
other_counter = {}
collect_collapses = {}

tumors_all_clustered = {}

for sgid in tumors_all.keys():
    #processing occurs by SGID. This is all within one sample.
    tumors_all_clustered[sgid]=[[],[]]
    collapsed_per_sgid[sgid]=0
    other_counter[sgid]=0
    collect_collapses[sgid]=[]
    bc_rc = zip(tumors_all[sgid][0],tumors_all[sgid][1], tumors_all[sgid][2])
    bc_rc=sorted(bc_rc, key=lambda x: x[1], reverse=True) #sorting barcodes by descending read count

    bc,rc,status = zip(*bc_rc)
    bc=list(bc)
    rc=list(rc)
    status=list(status)
    # run through your barcodes... for each barcode, check barcodes AFTER IT (in terms of readcount), within that sgid.
    print("working on barcodes associated with sgid {}\n".format(sgid))
    for i,val in enumerate(bc):
        for j in range(i+1,len(bc)):

            if rc[i]>20:
                if status[j] == 1: #if this tumor hasn't already been collapsed into another tumor
                    if rc[j] <threshold*rc[i]: #I might want to collapse j into i because it is enough smaller
                        if len(val) == len(bc[j]):
                            if hamming_distance(val,bc[j])<=int(parameter_info.bc_similarity_threshold): #if these barcodes are very similar
                                other_counter[sgid]+=1
                                collect_collapses[sgid].append(rc[j])
                                print("big tumor is {},{}. Little tumor getting collapsed is {},{}\n".format(val,rc[i],bc[j],rc[j]))
                    #I want to collapse j into i
                                rc[i] += rc[j]
                                status[j] = 0
    collapsed_per_sgid[sgid] = len(rc) - sum(status)
    for i,bc in enumerate(bc):
        if status[i] ==1:
            tumors_all_clustered[sgid][0].append(bc)
            tumors_all_clustered[sgid][1].append(rc[i])


for sgid in tumors_all_clustered.keys():
    for i in range(len(tumors_all_clustered[sgid][0])):
        output.write("{},{},{}\n".format(sgid, tumors_all_clustered[sgid][0][i], tumors_all_clustered[sgid][1][i]))

end_time = time.monotonic()
td = timedelta(seconds=end_time - start_time)
record.write("Clustering took {}\n".format(td))


cluster_summmary_name = root + project_name + "/" + parameter_name + "/summary/clustering_summaries/" + naming_stub + "_clustering_summary.txt"
cluster_out = open(cluster_summmary_name, 'wt')
cluster_out.write("Sample,sgid,tumors_in,tumors_collapsed,percent_tumors_collapsed,reads_in,reads_collapsed,percent_reads_collapsed\n")
for sgid, value in tumors_all_in.items():
    cluster_out.write("{},{},{},{},{},{},{},{}\n".format(sample, sgid, value, collapsed_per_sgid[sgid], float(collapsed_per_sgid[sgid])/float(value), reads_all_in[sgid], sum(collect_collapses[sgid]), float(sum(collect_collapses[sgid]))/float(reads_all_in[sgid])))


record.write("Going into clustering, the breakdown of tumors by sgid was:\n")

for sgid,val in tumors_all_in.items():
    record.write("{}: {}\n".format(sgid,val))
record.write("There were this many tumors collapsed per sgid:\n")
for sgid,collapsed in collapsed_per_sgid.items():
    record.write("{}, # collapsed: {}\n".format(sgid,collapsed))
record.close()
