
'''
First step of Tuba-seq pipeline: processes fastq files to count reads that map to each unique sgID-barcode combination.

'''
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
from helper_functions import unique, make_regexes, make_sgIDDict, getsgID, hamming_distance, revcom, check_indices, get_bad_indices, read_project_info

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

sgIDDict = make_sgIDDict(sgids, sgRNAs) #keys are sequences, values are genes targeted
sgIDDict_invert = {value : key for (key, value) in sgIDDict.items()}

########################################################################################
# 3. Set up output Files
########################################################################################

project_info_file = root + "/tubaseq_project_files/" + project_name + "_project_file.txt"

naming_stub = project_name + "_" + parameter_name + "_" + sample
record_name = root + "/" + project_name + "/" + parameter_name + "/records/" + naming_stub + "_record.txt"
reject_read_name = root + "/" + project_name + "/" + parameter_name + "/rejected/" + naming_stub + "_rejects.txt"
output_name = root + "/" + project_name + "/" + parameter_name + "/raw_counts/" + naming_stub + "_counts.txt"
record = open(record_name, 'wt')
reject = open(reject_read_name, 'wt')
output = open(output_name, 'wt')

project_info = read_project_info(project_info_file)

bad_indices = get_bad_indices(sample, project_info.sample_to_indices, project_info.sample_to_lane, project_info.lane_to_samples)

reject.write("Reason,line_seq1,line_seq2_comp\n")

########################################################################################
# 4. Set up regular expressions for testing
########################################################################################

bc_lengths=spikein_barcode_lengths
bc_lengths.append(barcode_length)
unique_bc_lengths = unique(bc_lengths)
regexes_to_check = make_regexes(unique_bc_lengths, R1_regex_looseness, R2_regex_looseness)


#########################################################################################
# 5. main loop, running through lines of file
########################################################################################

File1 = sample_path
File2 =File1.replace('_R1','_R2')

record.write("Tuba seq pipeline\n")
record.write("Beginning Part 1 of Tuba seq pipeline: counting reads\n ")
record.write("Timestamp: {}\n".format(now.strftime("%Y-%m-%d %H:%M")))
record.write("Project: {}, Parameters: {}\n".format(project_name, parameter_name))
record.write("Sample: {},\n\n".format(sample))
record.write("Input files:\nRead1: {}\nRead2: {}\n\n".format(File1, File2))
record.write("Output files:\nPrimary: {}\nRejected Reads: {}\n Record: {}\n".format(output_name,reject_read_name, record_name))
record.write("Regexes being searched for in each read:\n")
for i in range(len(regexes_to_check[0])):
    record.write("R1: {}, R2: {}\n".format(regexes_to_check[0][i], regexes_to_check[1][i]))

f1 = gzip.open(File1,'rt')
f2 = gzip.open(File2,'rt')

None_seq_dict = {}
sgIDBCdict = {}
indexdict = {}
line_seq1 = 1

total_reads = 0
n_tumors = 0
accepted_reads = 0
rejected_reads = 0

R1_didnt_match = 0
R2_didnt_match = 0
regex_rejects = 0
BC_match = 0
BC_match_fail = 0

prelim_none_counts = 0
prelim_good_counts = 0


index_tracker = {"Exact match": 0, "Index hop" : 0, "Inexact match" : 0}
index_tracker_full = {}

while (line_seq1):
    line1 = f1.readline().rstrip() # skip the first line
    print("line 1 is {}\n".format(line1))
    indices = line1.split(":")[-1]
    print("indices are {}\n".format(indices))
    if len(indices)>1:
        try:
            index_tracker[check_indices(indices, project_info.sample_to_indices[sample], bad_indices, index_encoding=project_info.index_encoding)] += 1
            if indices in index_tracker_full.keys():
                index_tracker_full[indices]+=1
            else:
                index_tracker_full[indices] = 1
        except:
            print("there is some problem with indices. Likely indices are not encoded as expeceted. You may be working with an older file from Illumina\n")


    line_seq1 = f1.readline().rstrip() 
    line1 = f1.readline().rstrip() # skip the third line
    line_qua1 = f1.readline().rstrip() # get the sequencing quality
    line2 = f2.readline().rstrip()
    line_seq2 = f2.readline().rstrip()
    line2 = f2.readline().rstrip()
    line_qua2 = f2.readline().rstrip()
    line_seq2_comp=revcom(line_seq2)
    total_reads+=1
    match = 0
    for i in range(len(regexes_to_check[0])): #CHECKING ALL THE regexes (R1 and R2)
        if (match<1):
            regexR1 = re.compile(regexes_to_check[0][i])
            regexR2 = re.compile(regexes_to_check[1][i])
            if regexR1.search(line_seq1) and regexR2.search(line_seq2_comp):
                    match += 1
                    k1 = regexR1.search(line_seq1) # align R1
                    k2 = regexR2.search(line_seq2_comp) # align R2
    if match==0:
        rejected_reads += 1
        regex_rejects += 1
        reject.write("R1orR2_didnt_match,{},{}\n".format(line_seq1,line_seq2_comp))

    if match>1:
        print("More than one of the possible regexes matched this read. This was not expected\n")
        sys.exit(2)
    if match ==1:
        R1sgID = k1.group(1) # upstream of R1
        R1BC = k1.group(2) # downstream of R1
        sgID_1 = getsgID(sgIDDict, R1sgID, 2) #note to self: getsgID will assign none to anything that is not in the expected sgids for this project.
        R2sgID = k2.group(1) # upstream of R1
        R2BC = k2.group(2) # downstream of R2_RC
        sgID_2 = getsgID(sgIDDict, R2sgID, 2)
        
        BC_dist=hamming_distance(R1BC,R2BC)

        if BC_dist <= int(bc_dist) and ('N' not in R1BC): #barcodes match
            BC_match += 1
            myKey = sgID_1 + "," + R1BC
            accepted_reads += 1
            if myKey in sgIDBCdict:
                sgIDBCdict[myKey] += 1
            else:
                n_tumors += 1
                sgIDBCdict[myKey] = 1

            if sgID_1 == "None":
                None_seq_dict[myKey]=R1sgID

        else:
            BC_match_fail += 1
            rejected_reads += 1
            reject.write("BC_match_fail,{},{}\n".format(line_seq1,line_seq2_comp))


f1.close()
end_time = time.monotonic()
td = timedelta(seconds=end_time - start_time)
percent_accepted_including_none = float(accepted_reads)/float(total_reads)
percent_rejected = float(rejected_reads)/float(total_reads)


record.write("Counting reads took: {}\n\n".format(td))
record.write("##########\n")
record.write("Reads processed: {}\n".format(total_reads))
record.write("Reads accepted including None tumors: {}\n".format(accepted_reads))
record.write("Accepted_including_none as percent: {:.2%}\n".format(percent_accepted_including_none))
record.write("Reads rejected (total): {}\n".format(rejected_reads))
record.write("Rejected_as_percent: {:.2%}\n\n".format(percent_rejected))

record.write("Reads discarded at regex matching stage: {}\n".format(regex_rejects))
record.write("Reads discarded at barcode matching stage: {}\n".format(BC_match_fail))
record.write("Reads where R1-R2 barcodes match: {}\n".format(BC_match))
record.write("Looking at indices now...\n")
record.write("For Sample {} the correct indices are {}\n".format(sample,project_info.sample_to_indices[sample]))
for i in index_tracker.keys():
   record.write("{}:{}\n".format(i,index_tracker[i]))

record.write("Counts for every index observed:\n")
indices_sorted = sorted(index_tracker_full.items(), reverse=True, key=lambda x: x[1])
for elem in indices_sorted:
   record.write("{},{}\n".format(elem[0],elem[1]))

record.write("##########\n")

output.write("sgID,BC,Count,sgIDseq\n")
for k,v in sorted(sgIDBCdict.items()):
    if k.strip().split(",")[0]=="None":
        prelim_none_counts+=int(v)
        output.write("{},{},{}\n".format(k,v,None_seq_dict[k]))
    else:
        prelim_good_counts+=int(v)
        output.write("{},{},{}\n".format(k,v,sgIDDict_invert[k.strip().split(",")[0]]))
record.write("Prelim None read count: {}\n".format(prelim_none_counts))
record.write("Prelim good read count: {}\n".format(prelim_good_counts))
