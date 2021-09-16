'''
Third step of the Tuba-seq pipeline- takes as input the output of filter_tumors.py
Converts the read number associated with each tumor (i.e., sgID-barcode combination) to an approximate
cell number based on the number of reads in the spike-in cell lines.
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
from helper_functions import hamming_distance, read_project_info, unique
import numpy as np


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
        print("Found sample, {}\n".format(sample))
    elif opt in ("--sample_path"):
        sample_path = arg
        print("Found sample_path, {}\n".format(sample_path))
    elif opt in ("--barcode_length"):
        barcode_length = arg
        print("Found barcode_length, {}\n".format(barcode_length))
    elif opt in ("--spikein_barcode_lengths"):
        spikein_barcode_lengths = arg.strip().split(",")
        print("Found spikein_barcode_lengths: {}\n".format(spikein_barcode_lengths))
    elif opt in ("--sgRNAs"):
        sgRNAs = arg.strip().split(",")
        print("found sgRNAs, {}\n".format(sgRNAs))
    elif opt in ("--sgids"):
        sgids = arg.strip().split(",")
        print("found sgids, {}\n".format(sgids))
    elif opt in ("--bc_dist"):
        bc_dist = arg
        print("found bc_dist, {}\n".format(bc_dist))
    elif opt in ("--R1_regex_looseness"):
        R1_regex_looseness = arg
        print("found R1_regex_looseness, {}\n".format(R1_regex_looseness))
    elif opt in ("--R2_regex_looseness"):
        R2_regex_looseness = arg
        print("found R2_regex_looseness, {}\n".format(R2_regex_looseness))
    elif opt in ("--root"):
        root = arg
        print("found root, {}\n".format(root))
    elif opt in ("--project_name"):
        project_name = arg
        print("found project_name, {}\n".format(project_name))
    elif opt in ("--parameter_name"):
        parameter_name = arg
        print("found parameter_name, {}\n".format(parameter_name))
    else:
        assert False, "unhandled option"


#get all the info for this project:
project_info_file = root + "/tubaseq_project_files/" + project_name + "_project_file.txt"
project_info = read_project_info(project_info_file) 
########################################################################################
# 3. Set up input/output Files
########################################################################################

naming_stub = project_name + "_" + parameter_name + "_" + sample
input_name = root + "/" + project_name + "/" + parameter_name + "/filtered/" + naming_stub + "_clustered.txt"
record_name = root + "/" + project_name + "/" + parameter_name + "/records/" + naming_stub + "_record.txt"
output_name = root + project_name + "/" + parameter_name + "/filtered/" + naming_stub + "_final.txt"
record = open(record_name, 'a+')
output = open(output_name, 'wt')
_input = open(input_name, 'rt')
_input.readline()

output.write("sgid,bc,rc,cellnum\n")

def confirm_spike_ins(largest_spike_ins, expected_spi_bc):
    '''
    This function confirms that the spike-ins with the largest read counts in this sample have the expected spike-in barcodes.
    '''
    for i in largest_spike_ins:
        fields = i.strip().split(",")
        bc = fields[1]
        if bc not in expected_spi_bc:
            print("UNEXPECTED SPIKE-IN BARCODE! {} is top 3 spike in at size {}, but expected barcodes were {}\n".format(bc, fields[2], expected_spi_bc))
            sys.exit(2)

def get_conversion_factor_unweighted(spi_bc_sizes, largest_spike_ins, expected_spi_bc):
    '''
    This function returns a conversion factor to go from read count to cell number for each tumor within a sample. 
    This is an unweighted average (i.e. each spike-in contributes equally to the factor), which I think is most sensible for cases where spike-ins of equal size were added.
    This code is flexible and should work for cases where spike-ins of different size were added (i.e. a ladder)
    '''
    bc_to_rc = {}
    sizes = []
    print("largest_spike_ins is {}\n".format(largest_spike_ins))
    for i,v in enumerate(largest_spike_ins):
        z = v.strip().split(",")
        bc_to_rc[z[1]] = z[2] # spike-in barcodes linked to read counts in the data
    print("expected_spi_bc is {}\n".format(expected_spi_bc))
    for i in expected_spi_bc:
        sizes.append(int(bc_to_rc[i])) #sizes holds the readcounts for the spike-ins IN THE ORDER THAT CORRESPONDS TO SPI_BC_SIZES
    to_average = []
    print("sizes is: {}\n".format(sizes))
    for i,val in enumerate(sizes):
        to_average.append(val/float(spi_bc_sizes[i]))

    print("to average is {}\n".format(to_average))
    print("not inverted factor is {}\n".format(np.mean(to_average)))

    factor = 1/np.mean(to_average)
    print("factor is: {}\n".format(factor))
    return(factor)

#########################################################################################
# 4. read in data
########################################################################################

tumors_in = 0
tumors_all = {}

for l in _input:
    fields = l.strip().split(",")
    sgid = fields[0]
    bc = fields[1]
    rc = fields[2]
    tumors_in += 1
    key = sgid + "," + bc + "," + rc
    if sgid in tumors_all.keys():
        tumors_all[sgid].append(key)
    else:
        tumors_all[sgid] = [key]

#########################################################################################
# 5. Analyze the spike-ins (get read count --> cell number conversion factor)
########################################################################################
print("project_info.sample_to_spike_bcs is {}\n".format(project_info.sample_to_spike_bcs[sample]))
print("project_info.sample_to_spike_sizes is {}\n".format(project_info.sample_to_spike_sizes[sample]))
expected_spi_bc = project_info.sample_to_spike_bcs[sample]
spi_bc_sizes =project_info.sample_to_spike_sizes[sample]
if len(expected_spi_bc) != len(spi_bc_sizes):
    print("ERROR WITH SPIKE-INS! you expect {} spike-ins, but only have sizes for {}\n".format(len(expected_spi_bc), len(spi_bc_sizes)))
    sys.exit(2)

#Spike in will automatically be in descending read count order
largest_spike_ins = []
for i in range(len(spi_bc_sizes)):
    largest_spike_ins.append(tumors_all["Spi"][i])

print("The largest spike-ins are: {}\n".format(largest_spike_ins))

confirm_spike_ins(largest_spike_ins, expected_spi_bc)
factor = get_conversion_factor_unweighted(spi_bc_sizes,largest_spike_ins, expected_spi_bc)

#########################################################################################
# 6. Re-write file with cell number
########################################################################################
for sgid in tumors_all.keys():
    for tumor in tumors_all[sgid]:
        fields=tumor.strip().split(",")
        cellnum = int(fields[2])*factor
        output.write("{},{}\n".format(tumor,cellnum))
