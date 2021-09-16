'''
Helper functions for Tuba-seq pipeline
'''

import numpy as np
import sys


def unique(list1):
    x = np.array(list1)
    return(np.unique(x))

def revcom(text):
    complement={'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    return ''.join(complement[base] for base in text[::-1])

def hamming_distance(str1, str2):
    if len(str1)!=len(str2):
        print("ALERT! strings must be of same length to calculate hamming distance, but str1 has length {} and str2 has length {}\n".format(len(str1), len(str2)))
    return sum(str1[i]!=str2[i] for i in range(len(str1)))

def get_quality_string(sample_str):
    list_of_quals = [int(ord(c))-33 for c in sample_str]
    return list_of_quals

def get_mean_quality_read(list_of_quals):
    mean_qual = np.average(list_of_quals)
    return mean_qual

def get_gc(seq):
    A=float(seq.count('A')/len(seq))
    G=float(seq.count('G')/len(seq))
    T=float(seq.count('T')/len(seq))
    C=float(seq.count('C')/len(seq))
    GC=float((seq.count('G')+seq.count('C'))/len(seq))
    #return([A,T,G,C,GC])
    return(GC)

def read_project_info(project_info_file):
    print("in read_project_info function. Project_info_file is {}\n".format(project_info_file))
    class project_info():
        def __init__(self, _project_name = None, _project_type = None, _sample_ids = None, _sample_paths = None, _sample_to_sample_paths = None, _genotypes = None, _gt_to_sample = None, _sample_to_gt = None, _gender_to_sample = None, _sample_to_gender= None, _sample_to_lungweight = None, _sample_to_titer = None, _control_gt = None, _sample_to_indices = None, _sample_to_injection_route = None, 
            _sgids = None, _genes = None, _sgRNAs = None, _inerts = None, _sgids_to_gene = None, _sgrnas_to_sgids = None,
            _sample_to_bc_length = None,
            _sample_to_n_spike = None, _sample_to_spike_bcs = None, _sample_to_spike_sizes = None, _sample_to_spike_bc_lengths = None,
            _sample_groupings = None,
            _mouse = None, _organ = None, _source = None, _mouse_to_sample = None, _sample_to_mouse = None, _organ_to_sample = None, _sample_to_organ = None, _source_to_sample = None, _sample_to_source = None,
            _sample_to_treatment = None, _treatment_to_sample = None,
            _sample_to_lane = None, _lane_to_samples = None,
            _donors = None,
            _index_encoding = None,
            _initial_pooling = None, _desired_rel_depth = None,
            _donor_mixing = None, _core_donors = None,
            _inteded_rel_depths = None, _sample_to_depth=None):
            
            self.project_name = _project_name
            self.project_type = _project_type

            self.sample_ids = _sample_ids
            self.sample_paths = _sample_paths
            self.sample_to_sample_paths = _sample_to_sample_paths
            self.genotypes = _genotypes
            self.gt_to_sample = _gt_to_sample
            self.sample_to_gt = _sample_to_gt
            self.gender_to_sample = _gender_to_sample
            self.sample_to_gender = _sample_to_gender
            self.sample_to_lungweight = _sample_to_lungweight
            self.sample_to_titer = _sample_to_titer
            self.control_gt = _control_gt
            self.sample_to_indices = _sample_to_indices
            self.sample_to_injection_route = _sample_to_injection_route
            self.sample_to_treatment = _sample_to_treatment
            self.treatment_to_sample = _treatment_to_sample

            self.sgids = _sgids
            self.genes = _genes
            self.sgRNAs = _sgRNAs
            self.inerts = _inerts
            self.sgids_to_gene = _sgids_to_gene
            self.sgRNAs_to_sgids = _sgrnas_to_sgids

            self.sample_to_bc_length = _sample_to_bc_length

            self.sample_to_n_spike = _sample_to_n_spike
            self.sample_to_spike_bcs = _sample_to_spike_bcs
            self.sample_to_spike_sizes = _sample_to_spike_sizes
            self.sample_to_spike_bc_lengths = _sample_to_spike_bc_lengths

            self.sample_groupings = _sample_groupings

            self.mouse = _mouse
            self.organ = _organ
            self.source = _source
            self.mouse_to_sample = _mouse_to_sample
            self.sample_to_mouse =_sample_to_mouse
            self.organ_to_sample = _organ_to_sample
            self.sample_to_organ = _sample_to_organ
            self.source_to_sample = _source_to_sample
            self.sample_to_source = _sample_to_source

            self.sample_to_lane = _sample_to_lane
            self.lane_to_samples = _lane_to_samples

            self.donors = _donors
            self.index_encoding = _index_encoding

            self.initial_pooling = _initial_pooling
            self.desired_rel_depth = _desired_rel_depth

            self.donor_mixing = _donor_mixing
            self.core_donors = _core_donors

            self.intended_rel_depths = _inteded_rel_depths
            self.sample_to_depth = _sample_to_depth


            self.populate_project_info(project_info_file)

        def populate_project_info(self, project_info_file):
            project_info={}
            f = open(project_info_file, 'rt')
            for l in f:

                #small addition 6/28/2021 to skip comments
                if l[0]!="#":

                    fields = l.strip().split("\t")
                    project_info[fields[0]]=fields[1]

            self.project_name = project_info["PROJECT"].strip()
            try:
                self.project_type = project_info["PROJECT_TYPE"].strip()
            except:
                print("project type not specified in project file!!!")
                sys.exit()

            self.sample_ids = project_info["SAMPLES"].strip().split(",")
            self.sample_paths = project_info["SAMPLES_PATHS"].strip().split(",")
            

            genotypes_raw = project_info["GENOTYPES"].strip().split(",")
            genders_raw = project_info["GENDERS"].strip().split(",")
            lungweights_raw = project_info["LUNGWEIGHTS"].strip().split(",")
            barcode_lengths_raw = project_info["BARCODE_LENGTH"].strip().split(",")
            nspike_raw = project_info["NSPIKE"].strip().split(",")
            indices_raw = project_info["INDICES"].strip().split(",")
            injection_routes_raw = project_info["INJECTION"].strip().split(",")
            lanes = project_info["LANES"].strip().split(",")
            treatments_raw = project_info["TREATMENT"].strip().split(",")
            depths_raw = project_info["DESIRED_RELATIVE_DEPTHS"].strip().split(",")

            spike_sizes_raw = project_info["SIZE_SPIKE"].strip().split(",")
            spike_bc_raw = project_info["SPIKE_BC"].strip().split(",")
            spike_bc_length_raw = project_info["SPIKE_BC_LENGTHS"].strip().split(",")
            titers_raw = project_info["TITERS"].strip().split(",")

            if self.project_type == "Transplant" or self.project_type == "Transplant_mixed": 
                self.core_donors = project_info["CORE_DONORS"].strip().split(",")
                donor_mixing_dict = {}
                if project_info["MIX_INFO"]!="NA":
                    donor_mixing_raw = project_info["MIX_INFO"].strip().split(";") 
                    for i in donor_mixing_raw:
                        x = i.strip().split(":")
                        identifier=x[0]
                        samps = x[1].strip().split(",")
                        donor_mixing_dict[identifier] = samps
                self.donor_mixing = donor_mixing_dict

            sample_groupings_dict = {}
            if project_info["SAMPLE_GROUPINGS"]!="NA":
                sample_groupings_raw = project_info["SAMPLE_GROUPINGS"].strip().split(";") 
                for i in sample_groupings_raw:
                    x = i.strip().split(":")
                    identifier = x[0]
                    samps = x[1].strip().split(",")
                    sample_groupings_dict[identifier] = samps
            self.sample_groupings = sample_groupings_dict

            self.control_gt = project_info["CAS9NEG_GT"].strip()

            self.sample_to_sample_paths = {}
            self.sample_to_gt = {}
            self.sample_to_gender = {}
            self.sample_to_lungweight = {}
            self.sample_to_bc_length = {}
            self.sample_to_n_spike = {}
            self.sample_to_titer = {}
            self.sample_to_indices = {}
            self.sample_to_injection_route = {}
            self.sample_to_treatment = {}

            self.sample_to_spike_bcs = {}
            self.sample_to_spike_sizes = {}
            self.sample_to_spike_bc_lengths = {}

            self.sample_to_depth = {}

            self.sample_to_lane = {}
            self.index_encoding = project_info["INDEX_ENCODING"].strip()

            for i in range(len(self.sample_ids)):
                self.sample_to_sample_paths[self.sample_ids[i]] = self.sample_paths[i]
                self.sample_to_gt[self.sample_ids[i]]=genotypes_raw[i]
                self.sample_to_gender[self.sample_ids[i]] = genders_raw[i]
                self.sample_to_lungweight[self.sample_ids[i]] = lungweights_raw[i]
                self.sample_to_bc_length[self.sample_ids[i]] = barcode_lengths_raw[i]
                self.sample_to_n_spike[self.sample_ids[i]] = nspike_raw[i]
                self.sample_to_titer[self.sample_ids[i]] = titers_raw[i]
                self.sample_to_indices[self.sample_ids[i]] = indices_raw[i]
                self.sample_to_injection_route[self.sample_ids[i]] = injection_routes_raw[i]
                self.sample_to_lane[self.sample_ids[i]] = lanes[i]
                self.sample_to_treatment[self.sample_ids[i]] = treatments_raw[i]
                self.sample_to_depth[self.sample_ids[i]] = depths_raw[i]

                self.sample_to_spike_bcs[self.sample_ids[i]] = spike_bc_raw[i].strip().split(":")
                self.sample_to_spike_sizes[self.sample_ids[i]] = list(map(int,spike_sizes_raw[i].strip().split(":")))
                self.sample_to_spike_bc_lengths[self.sample_ids[i]] = spike_bc_length_raw[i].strip().split(":")



            self.gt_to_sample = {}
            self.lane_to_samples = {}
            self.treatment_to_sample= {}

            for i in range(len(lanes)):
                if lanes[i] in self.lane_to_samples.keys():
                    self.lane_to_samples[lanes[i]].append(self.sample_ids[i])
                else:
                    self.lane_to_samples[lanes[i]] = [self.sample_ids[i]]

            for i in range(len(genotypes_raw)):
                if genotypes_raw[i] in self.gt_to_sample.keys():
                    self.gt_to_sample[genotypes_raw[i]].append(self.sample_ids[i])
                else:
                    self.gt_to_sample[genotypes_raw[i]]=[self.sample_ids[i]]  
            for i in range(len(treatments_raw)):
                if treatments_raw[i] in self.treatment_to_sample.keys():
                    self.treatment_to_sample[treatments_raw[i]].append(self.sample_ids[i])
                else:
                    self.treatment_to_sample[treatments_raw[i]]=[self.sample_ids[i]]  

            self.genotypes = unique(genotypes_raw)

            self.sgids = project_info["SGIDS"].strip().split(",")
            self.sgRNAs = project_info["SGRNAS"].strip().split(",")
            genes_raw = project_info["GENES"].strip().split(",")

            self.sgids_to_gene = {}
            for i in range(len(self.sgids)):
                self.sgids_to_gene[self.sgRNAs[i]] = genes_raw[i]
            self.genes = unique(genes_raw)

            self.inerts = project_info["INERT"].strip().split(",")
            self.sgRNAs_to_sgids = {}

            for i in range(len(self.sgRNAs)):
                self.sgRNAs_to_sgids[self.sgRNAs[i]] = self.sgids[i]

            if self.project_type == "Nano":
                self.initial_pooling = project_info["INITIAL_POOLING_PROP"].strip().split(",")
                self.desired_rel_depth = project_info["DESIRED_RELATIVE_DEPTHS"].strip().split(",")

            if self.project_type == "Transplant" or self.project_type == "Transplant_mixed": ####these parameters are only relevant to transplantation experiments; will just be default None for other exp.
                self.sample_to_mouse = {}
                self.sample_to_organ = {}
                self.sample_to_source = {}


                mice_raw = project_info["MOUSE"].strip().split(",")
                organs_raw = project_info["ORGAN"].strip().split(",")
                sources_raw = project_info["SOURCE"].strip().split(",")
                self.source =sources_raw
                self.donors = project_info["DONORS"].strip().split(",")
                #print("mice_raw {} organs raw {}, sources_raw {}\n".format(mice_raw, organs_raw, sources_raw))
                for i in range(len(self.sample_ids)):
                    self.sample_to_mouse[self.sample_ids[i]]=mice_raw[i]
                    self.sample_to_organ[self.sample_ids[i]]=organs_raw[i]
                    self.sample_to_source[self.sample_ids[i]]=sources_raw[i]

                self.mouse_to_sample = {}
                self.organ_to_sample = {}
                self.source_to_sample = {}
                for i in range(len(mice_raw)):
                    if mice_raw[i] in self.mouse_to_sample.keys():
                        self.mouse_to_sample[mice_raw[i]].append(self.sample_ids[i])
                    else:
                        self.mouse_to_sample[mice_raw[i]]=[self.sample_ids[i]] 
                    if organs_raw[i] in self.organ_to_sample.keys():
                        self.organ_to_sample[organs_raw[i]].append(self.sample_ids[i])
                    else:
                        self.organ_to_sample[organs_raw[i]]=[self.sample_ids[i]]
                    if sources_raw[i] in self.source_to_sample.keys():
                        self.source_to_sample[sources_raw[i]].append(self.sample_ids[i])
                    else:
                        self.source_to_sample[sources_raw[i]]=[self.sample_ids[i]] 

    return(project_info())


def read_parameter_info(parameter_info_file):
    class parameter_info():
        def __init__(self, _parameter_set_name = None, _dist_between_bc_for_read = None, _R1_regex_looseness = None, _R2_regex_looseness = None, _read_error_size_threshold = None, _bc_similarity_threshold = None, _gc_bias_correction = None, _min_t_size = None,
            _bundle_dist = None, _n_bundle= None, _contam_corr = None, _contam_removal_threshold = None):
            self.dist_between_bc_for_read = _dist_between_bc_for_read
            self.parameter_set_name = _parameter_set_name
            self.R1_regex_looseness = _R1_regex_looseness
            self.R2_regex_looseness = _R2_regex_looseness
            self.read_error_size_threshold = _read_error_size_threshold
            self.bc_similarity_threshold = _bc_similarity_threshold 
            self.gc_bias_correction = _gc_bias_correction
            self.min_t_size = _min_t_size
            self.bundle_dist = _bundle_dist
            self.n_bundle = _n_bundle
            self.contam_corr = _contam_corr
            self.contam_removal_threshold = _contam_removal_threshold
            self.populate_parameter_info(parameter_info_file)

        def populate_parameter_info(self, parameter_info_file):
            parameter_info={}
            f = open(parameter_info_file, 'rt')
            for l in f:
                fields = l.strip().split("\t")
                parameter_info[fields[0]]=fields[1]

            self.parameter_set_name = parameter_info["PARAMETER_SET"].strip()
            self.dist_between_bc_for_read = parameter_info["DIST_BETWEEN_BC_FOR_READ"].strip()
            self.R1_regex_looseness = parameter_info["R1_REGEX_LOOSENESS"].strip()
            self.R2_regex_looseness = parameter_info["R2_REGEX_LOOSENESS"].strip()
            self.read_error_size_threshold = parameter_info["READ_ERROR_SIZE_THRESHOLD"].strip()
            self.bc_similarity_threshold = parameter_info["BC_SIMILARITY_THRESHOLD"].strip()
            self.gc_bias_correction = parameter_info["GC_BIAS"].strip()
            self.min_t_size = float(parameter_info["MIN_SIZE_CUTOFF"].strip())
            self.bundle_dist = float(parameter_info["BUNDLE_DIST"].strip())
            self.n_bundle = int(parameter_info["N_BUNDLE".strip()])
            self.contam_corr = parameter_info["CONTAMINATION_CORRECT"].strip()
            self.contam_removal_threshold = parameter_info["CONTAMINATION_REMOVAL_THRESHOLD"].strip()
    return(parameter_info())

def write_input_file(root, project_info, parameter_info):
    fname = root + "/tubaseq_inp_files/" + project_info.project_name + "_" + parameter_info.parameter_set_name + ".inp"
    f = open(fname, 'wt')
    for sample in project_info.sample_ids:
        f.write("--sample={} --sample_path={} --barcode_length={} --spikein_barcode_lengths={} --sgRNAs={} --sgids={} --bc_dist={} --R1_regex_looseness={} --R2_regex_looseness={} --root={} --project_name={} --parameter_name={}\n".format(sample, project_info.sample_to_sample_paths[sample], project_info.sample_to_bc_length[sample], ",".join(project_info.sample_to_spike_bc_lengths[sample]), ",".join(project_info.sgRNAs), ",".join(project_info.sgids), parameter_info.dist_between_bc_for_read, parameter_info.R1_regex_looseness,parameter_info.R2_regex_looseness,  root, project_info.project_name, parameter_info.parameter_set_name))

def make_regexes(bc_lengths, R1_regex_looseness, R2_regex_looseness):
    forward_structure_20 = "GA" + "(........)" + "(GC.....TA.....GC.....TA.....GC)" + "(ATGCCCA){e<"+str(R1_regex_looseness) +"}"
    forward_structure_15 = 'GA' + '(........)' + '(AA.....TT.....AA.....)' + '(ATGCCCA){e<' + str(R1_regex_looseness) + "}"
    reverse_structure_20 ='(........)' +'(GC.....TA.....GC.....TA.....GC)' + '(ATGCCCA){e<'+ str(R2_regex_looseness) +"}"
    reverse_structure_15 = '(........)' +'(AA.....TT.....AA.....)' + '(ATGCCCA){e<'+ str(R2_regex_looseness) +"}"
    
    R1_regex_structures = {20: forward_structure_20, "20":forward_structure_20, 15:forward_structure_15, "15":forward_structure_15}
    R2_regex_structures = {20:reverse_structure_20, "20":reverse_structure_20, 15:reverse_structure_15, "15":reverse_structure_15}
    regex_R1 = []
    regex_R2 = []
    for bc in bc_lengths:
        if bc ==20 or bc=="20":
            regex_R1.append(R1_regex_structures[20])
            regex_R2.append(R2_regex_structures[20])
        elif bc ==15 or bc =="15":
            regex_R1.append(R1_regex_structures[15])
            regex_R2.append(R2_regex_structures[15])
        else:
            "NEED TO DEFINE ADDITIONAL BC STRUCTURES- don't have anything for length {}\n".format(bc)
    return([regex_R1,regex_R2])

def make_sgIDDict(sgids, sgRNAs):
    sgIDdict = {}
    for i, sgid in enumerate(sgids):
        sgIDdict[sgid] = sgRNAs[i]
    return(sgIDdict)

def getsgID(sgIDdict, sgid, sgiddist):
    '''flexbile getsgID returns match if sgid is within sgiddist
    from valid sgid'''
    for valid_sgid in sgIDdict.keys():
        if hamming_distance(valid_sgid, sgid)<sgiddist:
            return(sgIDdict[valid_sgid])

    return("None")

def check_master_sgid_file(master_file, sgRNAs_this_project, project):

    sgids_other_projects = {}
    f = open(master_file, 'rt')
    f.readline()
    for l in f:
        fields = l.strip().split("\t")
        #print("fields is {}\n".format(fields))
        project_name = fields[0]
        target_gene = fields[1]
        #guide_version = fields[2]
        guide_id = fields[2]
        sgRNA = fields[3]
        if project_name != project and sgRNA not in sgRNAs_this_project:
            key = project_name + "_" + target_gene
            if sgRNA in sgids_other_projects.keys():
                sgids_other_projects[sgRNA].append(key)
            else:
                sgids_other_projects[sgRNA] = [key]
    return(sgids_other_projects)

def get_bad_indices(sample, index_dict, sample_to_lane, lane_to_sample):
    '''get a list of indices from other samples in the lane to use to detect index hopping'''
    bad_indices = [] #these are indices that are come from the other samples in your lane
    samples_in_lane = lane_to_sample[sample_to_lane[sample]]
    for s in samples_in_lane:
        if s != sample:
            bad_indices = bad_indices + index_dict[s].strip().split(":")
    return bad_indices

def check_indices(index_info_fastq, correct_index, bad_indices, index_encoding):
    '''function is passed index info for a given read (in fastq formatting)
    and decides whether it matches the correct pair of indices for that sample
    The function will return Exact match if it is an exact match for the correct index,
    Inexact match if it is not an exact match, and Index hop if the pair contains indices
    that belong to multiple samples'''
    #ATCACG+AGCGCTAG
    info = index_info_fastq.strip().split("+")
    print("info is {}\n".format(info))
    if index_encoding == "revcomp:revcomp":
        reverse = revcom(info[0])
        forward=revcom(info[1])
    elif index_encoding == "norm:revcomp":
        reverse = revcom(info[0])
        forward = info[1]
    elif index_encoding == "revcomp:norm":
        reverse = info[0]
        forward = revcom(info[1])
    elif index_encoding == "norm:norm":
        reverse = info[0]
        forward = info[1]
    else:
        print("ERROR! the index encoding you specified ({}) is not a recognized option.\n".format(index_encoding))
        sys.exit(1)
    test = forward + ":" + reverse
    print("correct index is {}\n".format(correct_index))
    print("test is {}\n".format(test))
    if test == correct_index:
        print("found a match\n")
        return "Exact match"
    else:
        if forward in bad_indices:
            return "Index hop"
        elif reverse in bad_indices:
            return "Index hop"
        else:
            return "Inexact match"

