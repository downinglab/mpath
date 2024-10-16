import os
import sys
import argparse

from tqdm import tqdm

import pandas as pd
import numpy as np

from scipy.stats import pearsonr
from scipy.stats import beta
from scipy.stats import t

import matplotlib.pyplot as plt

from joblib import Parallel, delayed, dump, load, parallel_backend

import subprocess

def blocks(files, size=65536):
	while True:
		b = files.read(size)
		if not b: break
		yield b

def create_sublists(data, M):
	# Create an empty list to store the sublists
	sublists = [data[i:i+M] for i in range(0, len(data), M)]
	return sublists

def calc_smc(pairs_1, pairs_2):
	# match results in 1, not match results in -1
	return np.sum(((pairs_1 * pairs_2) + 1) / 2) / len(pairs_1)
	
def calc_inverse_smc(pairs_1, pairs_2):
	# match results in 1, not match results in -1
	return np.abs(np.sum(((pairs_1 * pairs_2) - 1) / 2)) / len(pairs_1)

def pearson_p(r, dist):
	"""
	Excises p-value calculation from the scipy pearsonr function,
	essentially speeding it up
	"""
	return 2*dist.cdf(-abs(r))
	
def calculate_p_value_t_dist(r, n):
	"""
	Calculate the p-value for a Pearson correlation coefficient given the number of data points.
	
	Parameters:
	r (float): The Pearson correlation coefficient.
	n (int): The number of data points used to calculate the correlation coefficient.
	
	Returns:
	float: The p-value corresponding to the hypothesis that the correlation coefficient is zero.
	
	Attempts to do the same thing as Matlab corrcoef function
	
	"""
	# Calculate the t-statistic
	try:
		t_statistic = r * np.sqrt((n - 2) / (1 - r**2))
	except:
		# if this is an error, r = 1, so p = 0
		return 0.
	
	# Calculate two-tailed p-value
	p_value = 2 * (1 - t.cdf(np.abs(t_statistic), df=n-2))
	
	return p_value
	
def get_p_values(slice):
	
	#pearson_ps = [pearson_p(r, dists[n]) for r,n in slice]
	pearson_ps = [calculate_p_value_t_dist(r, n) for r,n in slice]
	
	return pearson_ps
	
def calc_pearson_fast(pairs_1, pairs_2):
	
	mean_pairs_1 = np.mean(pairs_1)
	mean_pairs_2 = np.mean(pairs_2)
	
	pearson_numerator = np.sum((pairs_1 - mean_pairs_1) * (pairs_2 - mean_pairs_2))
	pearson_denominator = np.sqrt(np.sum(np.square(pairs_1 - mean_pairs_1)) * np.sum(np.square(pairs_2 - mean_pairs_2)))
	
	pearson_r = pearson_numerator / pearson_denominator
	
	return pearson_r
	
def correlate_slice(slice, use_full_matrix):
	df_corr_records = []
	
	k = 0
	
	for read in slice:
		
		if len(read['meth_values']) < 3:
			continue
		
		read_id = read['read_id']
		meth_values = np.array(read['meth_values'])
		wgbs_values = np.array(read['wgbs_values'])
		positions = np.array(read['positions'])
		
		meth_ratio = np.mean(meth_values)
		wgbs_distance = np.sum(np.square(wgbs_values - meth_values)) / len(meth_values)
		
		# make 0 into -1, 1 stays as 1. makes smc calculation easier
		# and does not influence pearson
		meth_values = meth_values * 2 - 1
		
		# make the 1d methylation list into 2d, we just repeat the values
		meth_2d = np.repeat(meth_values[None,:], len(meth_values), axis=0)
		meth_2d_t = meth_2d.T
		
		if use_full_matrix:
			filter = np.full((len(meth_values), len(meth_values)), True)
		else:
			# real meth_values will never be 0 (since -1 & 1) so this keeps all upper triangle
			filter = np.logical_not(np.triu(meth_2d, k=k) == 0)
		
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
			distances = np.triu(np.abs(np.subtract.outer(positions, positions)), k=k)
		
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
	
	return df_corr_records


# python /home/data/nlaszik/nanopore/code/longread_metrics_parallel.py -path_input_bed /home/data/atrinh6/nanopore/240119_HUES64_100X/LIB1/20240119_1352_P2S-01599-A_PAU70119_f639b486/readlevelmeth_wgbs.bed -path_output_csv /home/data/Shared/shared_datasets/longread/nlaszik/100x_lib1.csv -p 8

# python /home/data/nlaszik/nanopore/code/longread_metrics_parallel.py -path_input_bed /home/data/atrinh6/nanopore/231208_BrdU500_PGvM/readlevel_meth/intersect_avgBrdU02/readlevel_meth_avgBrdU02ONLY/codeinput/uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq.bed -path_output_csv /home/data/Shared/shared_datasets/longread/nlaszik/uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq_full_matrix.csv -p 64 --use_full_matrix

# python /home/data/nlaszik/nanopore/code/longread_metrics_parallel.py -path_input_bed /home/data/atrinh6/nanopore/231208_BrdU500_PGvM/readlevel_meth/intersect_avgBrdU02/readlevel_meth_avgBrdU02ONLY/codeinput/uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq.bed -path_output_csv /home/data/Shared/shared_datasets/longread/nlaszik/uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq_full_matrix.csv -p 8 --use_full_matrix

# python /home/data/nlaszik/nanopore/code/longread_metrics_parallel.py -path_input_bed /home/data/atrinh6/nanopore/231208_BrdU500_PGvM/readlevel_meth/intersect_avgBrdU02/readlevel_meth_avgBrdU02ONLY/codeinput/DS1000_uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq.bed -path_output_csv /home/data/Shared/shared_datasets/longread/nlaszik/DS1000_uniq_sameStartEnd_PG_B500_16h_readlevelmeth_avgBrdU02ONLY_WGBS_uniq.csv -p 8 --use_full_matrix

parser = argparse.ArgumentParser()

parser.add_argument("-path_input_bed", required = True)
parser.add_argument("-path_output_csv", required = True)
parser.add_argument("-p", required = False)
parser.add_argument("--use_full_matrix", action='store_true')

# get arguments
args = parser.parse_args()

path_input_bed = args.path_input_bed
path_output_csv = args.path_output_csv
p = args.p
use_full_matrix = args.use_full_matrix

if p:
	num_processes = int(p)
else:
	num_processes = 1

distance_bins = [[0, 100], [100, 1000], [1000, 5000], [5000, 10000]]
read_dict = {}

print('counting cpgs...')
result = subprocess.run(["wc", "-l", path_input_bed], stdout=subprocess.PIPE, text=True, check=True)
num_cpgs = int(result.stdout.split()[0])

batch_size = min(num_cpgs, int(1e8))

num_batches = int(num_cpgs // batch_size)

# from here on we are batching
# keeping the file handle open allows us to restart
lines_remaining = True
with open(path_input_bed, 'r') as fh:
	
	i_batch = 0
	read_dict = {}
	last_read_id = ''
	
	while lines_remaining:
			
		print(f'loading data for batch {i_batch + 1} out of {num_batches}...')
		num_cpgs_in_batch = 0
		lines_remaining = False
		
		for line in tqdm(fh, total=batch_size):
			
			line_list = line.split('\t')
			read_id = line_list[4]
			
			if num_cpgs_in_batch > batch_size and read_id != last_read_id:
				lines_remaining = True
				break
			
			if last_read_id != read_id:
				read_dict[read_id] = {'meth_values': [], 'positions': [], 'wgbs_values': [], 'read_id': read_id}
			
			read_dict[read_id]['meth_values'].append(float(line_list[5]))
			read_dict[read_id]['positions'].append(int(line_list[1]))
			read_dict[read_id]['wgbs_values'].append(float(line_list[6]))
			
			last_read_id = read_id
			num_cpgs_in_batch += 1
			
		# load just this batch and slice
		read_list = [read_dict[read_id] for read_id in read_dict]
		num_lines = len(read_list)
		slice_size = num_lines // num_processes
		slices = [read_list[i:i+slice_size] for i in range(0, num_lines, slice_size)]
		print(f'processing, slice size: {slice_size}...')
		with parallel_backend("loky", inner_max_num_threads=2):
			batch_df_corr_records = Parallel(n_jobs=num_processes, verbose=100, pre_dispatch="all")(delayed(correlate_slice)(slice, use_full_matrix) for slice in slices)
			
		df_corr_records = [item for sublist in batch_df_corr_records for item in sublist]
		df_corr = pd.DataFrame.from_records(df_corr_records)
		
		list_pearson_r = list(df_corr['pearson_r'])
		list_num_pairs = list(df_corr['num_pairs'])
		
		unique_n_pairs = list(set(list_num_pairs))
		
		zip_list = list(zip(list_pearson_r, list_num_pairs))
		
		print(len(unique_n_pairs))
		#print('getting beta distributions...')
		#dists = {n: beta(n/2 - 1, n/2 - 1, loc=-1, scale=2) for n in tqdm(unique_n_pairs)}
		
		print('getting p values...')
		num_lines = len(zip_list)
		slice_size = num_lines // num_processes
		slices = [zip_list[i:i+slice_size] for i in range(0, num_lines, slice_size)]
		
		with parallel_backend("loky", inner_max_num_threads=2):
			all_p_values = Parallel(n_jobs=num_processes, verbose=100, pre_dispatch="all")(delayed(get_p_values)(slice) for slice in slices)
		
		df_corr['pearson_p'] = [item for sublist in all_p_values for item in sublist]
		
		print(df_corr)
		path_output_csv_batch = path_output_csv.replace('.csv', f'.batch{i_batch}.csv')
		df_corr.to_csv(path_output_csv_batch, index=False)
		
		# since file handle iterator cannot go back, we use the last value to make sure nothing is missed
		read_dict = {}
		read_dict[read_id] = {'meth_values': [], 'positions': [], 'wgbs_values': [], 'read_id': read_id}
		read_dict[read_id]['meth_values'].append(float(line_list[5]))
		read_dict[read_id]['positions'].append(int(line_list[1]))
		read_dict[read_id]['wgbs_values'].append(float(line_list[6]))
		
		del batch_df_corr_records
		del df_corr_records
		
		i_batch += 1

