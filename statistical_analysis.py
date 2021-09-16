import numpy as np
import random
import sys

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2)) 

def LN_mean(t_sizes):
	'''Calculates the log-normal mean over a list of tumor sizes.

	Arguments:
	-t_sizes: list of tumor sizes.
	'''
	LN_x = np.log(t_sizes)
	m = np.mean(LN_x)
	v = np.var(LN_x, ddof=1)
	output = np.exp(m + 0.5*v) 
	return(output)

def fdr_correct_pval_by_sgid(pval_dict):
	'''FDR-correct p-values.

	Arguments:
	-pval_dict: dictionary where keys are sgIDs and values are uncorrected p-values.
	'''
	import statsmodels.stats.multitest as smm
	fdr_corr_dict = {}
	sgids = []
	p = []

	for sgid in pval_dict.keys():
		sgids.append(sgid)
		p.append(pval_dict[sgid])
	pvals_corr = list(smm.fdrcorrection(np.array(p))[1])

	for i, sgid in enumerate(sgids):
		fdr_corr_dict[sgid] = pvals_corr[i]
	return(fdr_corr_dict)


def input_data(infname, list_of_samples, sgids_to_exclude):
	'''Reads in tumor file (output of tumor_processing.py) and a list of samples and 
	returns a nested dictionary where keys are sample --> sgid --> list of tumor sizes (cell #)

	Arguments:
	-infname is the tumor file (output of tumor_processing.py)
	-list_of_samples: list of samples to include in the analysis (usually, all samples).
=	-sgids_to_exclude: list of sgIDs to exclude from analysis (usually, none).
	'''
	sgids_to_exclude.append("Spi")
	mouse_dict = {}
	f = open(infname, 'rt')
	f.readline()
	for l in f:
		fields = l.strip().split(",")
		mouse = fields[5]
		sgid = fields[0]
		if sgid not in sgids_to_exclude and mouse in list_of_samples:
			GT = fields[6]
			CN = float(fields[3])
			if mouse in mouse_dict.keys():
				if sgid in mouse_dict[mouse].keys():
					mouse_dict[mouse][sgid].append(CN)
				else:
					mouse_dict[mouse][sgid] = [CN]
			else:
				mouse_dict[mouse] = {}
				mouse_dict[mouse][sgid] = [CN]
	f.close()
		
	return(mouse_dict)


def merge_mice(mouse_dict):
	'''Collects tumors across mice into a single dataset. Returns a dictionary where the keys are sgids and values are lists of tumor sizes.

	Arguments:
	-mouse_dict: not-merged tumor dataset. (This is a nested dictionary with the structure samples >> sgIDs >> list of cell numbers).

	'''
	merged_mouse_dict = {}

	for mouse in mouse_dict.keys():
		for sgid in mouse_dict[mouse].keys():
			if sgid in merged_mouse_dict.keys():
				merged_mouse_dict[sgid].extend(mouse_dict[mouse][sgid])
			else:
				merged_mouse_dict[sgid] = mouse_dict[mouse][sgid]

	return(merged_mouse_dict)


def getP(stats):
	'''Calculates a two-sided p-value by comparing the bootstrapped values of an inert-adjusted statistic 
	to its value under the null hypothesis of no effect of tumor genotype (i.e., 1).

	'''
	n_boot = len(stats)
	boolean_stat = []
	for i in stats:
		boolean_stat.append(i>1)
	return(2*min((sum(boolean_stat)/n_boot), (1-(sum(boolean_stat)/n_boot))))


def adj_percentiles(data, percentiles, inerts):
	'''Calculates relative (inert-adjusted) percentiles of the tumor size distribution for each sgID.

	Arguments:
	-data: tumor dataset (This is a dictionary where the keys are sgids and values are lists of tumor sizes)
	-percentiles: list of percentiles to evaluate.
	-inerts = list of inert sgIDs.

	'''
	perc_out = {}
	inert_tumors = []
	for i in inerts:
		if i in data.keys():
			inert_tumors.extend(data[i])
	print("in adj_percentiles, it has {} inert tumors\n".format(len(inert_tumors)))
	for sgid in data.keys():

		perc_out[sgid] = {}
		tumors = data[sgid]
		print("working on sgid {}, it has {} tumors, {} inert tumors\n".format(sgid,len(tumors), len(inert_tumors)))
		for p in percentiles:
			perc = p*100
			perc_out[sgid][p] = np.percentile(tumors, perc) / np.percentile(inert_tumors,perc)
		perc_out[sgid]["LNmean"] = LN_mean(tumors) / LN_mean(inert_tumors)

	return(perc_out)


def adj_burden(g_data, kt_data, inerts):
	'''Calculates relative (inert-adjusted) tumor burden statistic.

	Arguments:
	-g_data: tumor data for focal genotype. (This is a dictionary where the keys are sgids and values are lists of tumor sizes)
	-kt_data: tumor data for KT mice (these mice lack Cas9, therefore they are used to assess the make-up of the viral pool)
	-inerts: list of inert sgIDs.

	'''
	burden_out = {}
	inert_tumors_g = []
	inert_tumors_KT = []
	for i in inerts:
		if i in g_data.keys():
			inert_tumors_g.extend(g_data[i])
		if i in kt_data.keys():
			inert_tumors_KT.extend(kt_data[i])

	for sgid in intersection(list(g_data.keys()),list(kt_data.keys())):
		burden_out[sgid] = {}
		tumors_g = g_data[sgid]
		tumors_KT = kt_data[sgid]
		ts_stat = sum(tumors_g) / sum(tumors_KT)
		inert_stat = sum(inert_tumors_g) / sum(inert_tumors_KT)
		burden_out[sgid] = ts_stat / inert_stat
	return(burden_out)


def adj_TN(g_data, kt_data, inerts):
	'''Calculates relative (inert-adjusted) tumor number statistic.

	Arguments:
	-g_data: tumor data for focal genotype.
	-kt_data: tumor data for KT mice (these mice lack Cas9, therefore they are used to assess the make-up of the viral pool).
	-inerts: list of inert sgIDs.

	'''
	TN_out = {}
	inert_tumors_g = []
	inert_tumors_KT = []
	for i in inerts:
		if i in g_data.keys():
			inert_tumors_g.extend(g_data[i])
		if i in kt_data.keys():
			inert_tumors_KT.extend(kt_data[i])
	for sgid in intersection(list(g_data.keys()),list(kt_data.keys())):
		TN_out[sgid] = {}
		tumors_g = g_data[sgid]
		tumors_KT = kt_data[sgid]
		ts_stat = len(tumors_g) / len(tumors_KT)
		inert_stat = len(inert_tumors_g) / len(inert_tumors_KT)
		TN_out[sgid] = ts_stat / inert_stat
	return(TN_out)


def bootstrap_mice(mouse_dict, n_mice):
	'''Bootstraps mice.

	Arguments:
	-mouse_dict: tumor_dataset (this is a nested dictionary with the structure samples >> sgIDs >> list of cell numbers).
	-n_mice: number of mice to sample (for bootstrapping, set equal to the original number of mice).

	'''
	bs_mouse_dict = {}
	mouse_list = list(mouse_dict.keys())
	bs_mouse_list = np.random.choice(mouse_list,n_mice, replace=True)

	for i,m in enumerate(bs_mouse_list):
		bs_mouse_dict[i] = mouse_dict[m]
	return(bs_mouse_dict)


def bootstrap_tumors(data):
	'''Bootstraps tumors within each mouse, by sgID.

	Arguments:
	-data: tumor dataset (this is a nested dictionary with the structure samples >> sgIDs >> list of cell numbers).

	'''
	new_data = {}
	for m in data.keys(): #for each mouse
		new_data[m] = {}
		for s in data[m].keys():
			t_indices = range(len(data[m][s]))
			new_t_indices = np.random.choice(t_indices,len(t_indices), replace=True)
			new_t = [data[m][s][i] for i in new_t_indices]
			new_data[m][s] = new_t
	return(new_data)

