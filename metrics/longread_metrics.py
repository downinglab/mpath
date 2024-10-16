import os
import sys
import argparse

from tqdm import tqdm

import pandas as pd
import numpy as np

from scipy.stats import pearsonr
from scipy.stats import beta

import matplotlib.pyplot as plt

def calc_smc(pairs_1, pairs_2):
	# match results in 1, not match results in -1
	return np.sum(((pairs_1 * pairs_2) + 1) / 2) / len(pairs_1)
	
def calc_inverse_smc(pairs_1, pairs_2):
	# match results in 1, not match results in -1
	return np.abs(np.sum(((pairs_1 * pairs_2) - 1) / 2)) / len(pairs_1)
	
def calc_pearson_fast(pairs_1, pairs_2):
	
	mean_pairs_1 = np.mean(pairs_1)
	mean_pairs_2 = np.mean(pairs_2)
	
	pearson_numerator = np.sum((pairs_1 - mean_pairs_1) * (pairs_2 - mean_pairs_2))
	pearson_denominator = np.sqrt(np.sum(np.square(pairs_1 - mean_pairs_1)) * np.sum(np.square(pairs_2 - mean_pairs_2)))
	
	pearson_r = pearson_numerator / pearson_denominator
	
	return pearson_r

# python /home/data/nlaszik/nanopore/code/longread_metrics.py -path_input_bed /home/data/atrinh6/nanopore/231208_BrdU500_PGvM/readlevel_meth/intersect_avgBrdU02/readlevel_meth_avgBrdU02ONLY/codeinput/uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq.bed -path_output_csv /home/data/Shared/shared_datasets/longread/nlaszik/uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq_full_matrix.csv --use_full_matrix

# python /home/data/nlaszik/nanopore/code/longread_metrics.py -path_input_bed /home/data/atrinh6/nanopore/231208_BrdU500_PGvM/readlevel_meth/intersect_avgBrdU02/readlevel_meth_avgBrdU02ONLY/codeinput/DS1000_uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq.bed -path_output_csv /home/data/Shared/shared_datasets/longread/nlaszik/DS1000_uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq.csv --use_full_matrix

parser = argparse.ArgumentParser()

parser.add_argument("-path_input_bed", required = True)
parser.add_argument("-path_output_csv", required = True)
parser.add_argument("--use_full_matrix", action='store_true')

# get arguments
args = parser.parse_args()

path_input_bed = args.path_input_bed
path_output_csv = args.path_output_csv
use_full_matrix = args.use_full_matrix

distance_bins = [[0, 100], [100, 1000], [1000, 5000], [5000, 10000]]

print('loading long read data...')
df_cpg = pd.read_csv(path_input_bed, sep='\t', names=['chrom', 'start', 'stop', 'strand', 'read_id', 'methylation', 'wgbs'])
df_cpg_records = df_cpg.to_dict('records')

print(df_cpg)

# construct first
print('grouping cpgs by read...')
read_ids = list(set(list(df_cpg['read_id'])))
read_dict = {read_id:{'meth_values':[], 'positions':[], 'wgbs_values':[]} for read_id in read_ids}
for record in tqdm(df_cpg_records):
	# to make multiplication work to do correlation, we must change 0 to -1, but keep 1
	read_dict[record['read_id']]['meth_values'].append(record['methylation'])
	read_dict[record['read_id']]['positions'].append(record['start'])
	read_dict[record['read_id']]['wgbs_values'].append(record['wgbs'])
	
print('correlating...')
df_corr_records = []
for read_id in tqdm(read_dict):
	
	if len(read_dict[read_id]['meth_values']) < 3:
		continue
	
	meth_values = np.array(read_dict[read_id]['meth_values'])
	wgbs_values = np.array(read_dict[read_id]['wgbs_values'])
	positions = np.array(read_dict[read_id]['positions'])
	
	meth_ratio = np.mean(meth_values)
	wgbs_distance = np.sum(np.square(wgbs_values - meth_values)) / len(meth_values)
	
	# make 0 into -1, 1 stays as 1. makes smc calculation easier
	# and does not influence pearson
	meth_values = meth_values * 2 - 1
	
	meth_2d = np.repeat(meth_values[None,:], len(meth_values), axis=0)
	meth_2d_t = meth_2d.T
	
	if use_full_matrix:
		filter = np.full((len(meth_values), len(meth_values)), True)
	else:
		# real meth_values will never be 0 (since -1 & 1) so this keeps all upper triangle
		filter = np.logical_not(np.triu(meth_2d, k=1) == 0)
	
	pairs_1 = meth_2d[filter]
	pairs_2 = meth_2d_t[filter]
	
	pearson_r = calc_pearson_fast(pairs_1, pairs_2)
	smc = calc_smc(pairs_1, pairs_2)
	inverse_smc = calc_inverse_smc(pairs_1, pairs_2)
	
	df_corr_records.append({'read_id': read_id, 'read_wgbs_distance': wgbs_distance, 'read_meth_ratio': meth_ratio, 'distance_bin': 'bin_all', 'num_pairs': len(pairs_1), 'smc': smc, 'uniformity': smc - inverse_smc, 'pearson_r': pearson_r})
	
	# get matrix of distance between pairs of cpgs
	if use_full_matrix:
		distances = np.abs(np.subtract.outer(positions, positions))
	else:
		distances = np.triu(np.abs(np.subtract.outer(positions, positions)), k=1)
	
	# distance bins
	for i_bin, distance_bin in enumerate(distance_bins):
		
		distance_filter = np.logical_and(distances > distance_bin[0], distances <= distance_bin[1])
		
		pairs_1 = meth_2d[distance_filter]
		pairs_2 = meth_2d_t[distance_filter]
		
		if len(pairs_1) > 1:
		
			pearson_r = calc_pearson_fast(pairs_1, pairs_2)
			smc = calc_smc(pairs_1, pairs_2)
			inverse_smc = calc_inverse_smc(pairs_1, pairs_2)
			
			df_corr_records.append({'read_id': read_id, 'read_wgbs_distance': wgbs_distance, 'read_meth_ratio': meth_ratio, 'distance_bin': f'bin_{i_bin}', 'num_pairs': len(pairs_1), 'smc': smc, 'uniformity': smc - inverse_smc, 'pearson_r': pearson_r})
			
		else:
			
			df_corr_records.append({'read_id': read_id, 'read_wgbs_distance': wgbs_distance, 'read_meth_ratio': meth_ratio, 'distance_bin': f'bin_{i_bin}', 'num_pairs': len(pairs_1), 'smc': np.nan, 'uniformity': np.nan, 'pearson_r': np.nan})
			
	# closest cpgs
	if use_full_matrix:
		# real meth_values will never be 0 (since -1 & 1) so this keeps all upper triangle
		upper_diag = np.diag(np.diag(meth_2d, k=1), k=1)
		lower_diag = np.diag(np.diag(meth_2d, k=-1), k=-1)
		closest_filter = np.logical_not(upper_diag == 0) + np.logical_not(lower_diag == 0)
	else:
		# real meth_values will never be 0 (since -1 & 1) so this keeps all upper triangle
		upper_diag = np.diag(np.diag(meth_2d, k=1), k=1)
		closest_filter = np.logical_not(upper_diag == 0)
		
	pairs_1 = meth_2d[closest_filter]
	pairs_2 = meth_2d_t[closest_filter]
	
	pearson_r = calc_pearson_fast(pairs_1, pairs_2)
	smc = calc_smc(pairs_1, pairs_2)
	inverse_smc = calc_inverse_smc(pairs_1, pairs_2)
	
	df_corr_records.append({'read_id': read_id, 'read_wgbs_distance': wgbs_distance, 'read_meth_ratio': meth_ratio, 'distance_bin': 'bin_closest', 'num_pairs': len(pairs_1), 'smc': smc, 'uniformity': smc - inverse_smc, 'pearson_r': pearson_r})
	
df_corr = pd.DataFrame.from_records(df_corr_records)

list_num_pairs = list(df_corr['num_pairs'])
list_pearson_r = list(df_corr['pearson_r'])

unique_n_pairs = list(set(list_num_pairs))

print(len(unique_n_pairs))
print('getting beta distributions...')
dists = {n: beta(n/2 - 1, n/2 - 1, loc=-1, scale=2) for n in tqdm(unique_n_pairs)}

print('getting p values...')
pearson_p = [2*dists[n].cdf(-abs(r)) for r, n in tqdm(zip(list_pearson_r, list_num_pairs))]

df_corr['pearson_p'] = pearson_p

print(df_corr)
df_corr.to_csv(path_output_csv, index=False)

