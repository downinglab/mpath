#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Calculates metrics for Nanopore long read data.

@author: nglaszik

"""

import os
import sys
import argparse
import random
import string
import math

from tqdm import tqdm

import pandas as pd
import numpy as np

from scipy.stats import t

from joblib import Parallel, delayed, dump, load, parallel_backend

import subprocess

def generate_random_string(length):
	"""Generates a random string of letters and digits."""
	characters = string.ascii_letters + string.digits
	return ''.join(random.choice(characters) for i in range(length))

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
	
	pearson_ps = [calculate_p_value_t_dist(r, n) for r,n in slice]
	return pearson_ps
	
def calc_pearson_fast(pairs_1, pairs_2):
	
	mean_pairs_1 = np.mean(pairs_1)
	mean_pairs_2 = np.mean(pairs_2)
	
	pearson_numerator = np.sum((pairs_1 - mean_pairs_1) * (pairs_2 - mean_pairs_2))
	pearson_denominator = np.sqrt(np.sum(np.square(pairs_1 - mean_pairs_1)) * np.sum(np.square(pairs_2 - mean_pairs_2)))
	
	if pearson_denominator > 0:
		pearson_r = pearson_numerator / pearson_denominator
	else:
		pearson_r = np.nan
	
	return pearson_r
	
def correlate_slice(slice, use_full_matrix, min_cpgs):
	
	df_corr_records = []
	
	for read in slice:
		
		if len(read['meth_values']) < min_cpgs:
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
			filter = np.logical_not(np.triu(meth_2d) == 0)
		
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
			distances = np.triu(np.abs(np.subtract.outer(positions, positions)))
		
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

# python ./metrics/longread_metrics_parallel.py -path_input_bed ./data/input.bed -path_output_csv ./data/output.csv -p 8 -min_cpgs 3 -bin_limits 0,100,1000,5000,10000 --use_full_matrix

parser = argparse.ArgumentParser()

parser.add_argument("-path_input_bed", required = True)
parser.add_argument("-path_output_csv", required = True)
parser.add_argument("-min_cpgs", required = True)
parser.add_argument("-bin_limits", required = True)
parser.add_argument("-batch_size", required = False, default=1e8)
parser.add_argument("-p", required = False, default=1)
parser.add_argument("--use_full_matrix", action='store_true')

# get arguments
args = parser.parse_args()

path_input_bed = args.path_input_bed
path_output_csv = args.path_output_csv
min_cpgs = int(args.min_cpgs)
bin_limits = args.bin_limits
batch_size = int(args.batch_size)
num_processes = int(args.p)
use_full_matrix = args.use_full_matrix

# construct distance_bins from string
bin_limits_list = [int(limit) for limit in bin_limits.split(',')]
distance_bins = [[bin_limits_list[i], bin_limits_list[i + 1]] for i in range(len(bin_limits_list) - 1)]

for i_bin, distance_bin in enumerate(distance_bins):
	print(f'bin_{i_bin}: {distance_bin}')

print('counting cpgs...')
result = subprocess.run(["wc", "-l", path_input_bed], stdout=subprocess.PIPE, text=True, check=True)
num_cpgs = int(result.stdout.split()[0])
num_cpgs_remaining = num_cpgs
print(f'processing {num_cpgs} cpgs')

# minimize batch size if more than number of cpgs
batch_size = min(num_cpgs, batch_size)
num_batches = int(math.ceil(num_cpgs / batch_size))

# read file and process in batches
batch_csv_paths = []
lines_remaining = True
with open(path_input_bed, 'r') as fh:
	
	i_batch = 0
	i_line = 0
	read_dict = {}
	last_read_id = ''
	
	while lines_remaining:
			
		print(f'loading data for batch {i_batch + 1} out of {num_batches}...')
		num_cpgs_in_batch = 0
		running_batch_size = min(batch_size, num_cpgs_remaining)
		lines_remaining = False
		
		for line in tqdm(fh, total=running_batch_size):
			
			line_list = line.split('\t')
			
			# exit if enough lines aren't detected
			# this is fine for entire file since iterator will stop once we reach end of file
			if len(line_list) < 7:
				print(f'Please include all columns in input BED, only {len(line_list)} detected at row {i_line}')
				sys.exit()
			
			position = int(line_list[1])
			read_id = line_list[4]
			meth_value = float(line_list[5])
			wgbs_value = float(line_list[6])
			
			# will stop including reads in batch once batch size is full AND we have finished reading data for an entire read
			if num_cpgs_in_batch > batch_size and read_id != last_read_id:
				lines_remaining = True
				break
			
			if last_read_id != read_id:
				read_dict[read_id] = {'meth_values': [], 'positions': [], 'wgbs_values': [], 'read_id': read_id}
			
			read_dict[read_id]['meth_values'].append(meth_value)
			read_dict[read_id]['positions'].append(position)
			read_dict[read_id]['wgbs_values'].append(wgbs_value)
			
			last_read_id = read_id
			num_cpgs_in_batch += 1
			i_line += 1
			
		num_cpgs_remaining = num_cpgs_remaining - num_cpgs_in_batch
			
		# load just this batch and slice it
		read_list = [read_dict[read_id] for read_id in read_dict]
		num_lines = len(read_list)
		slice_size = num_lines // num_processes
		slices = [read_list[i:i+slice_size] for i in range(0, num_lines, slice_size)]
		print(f'processing, slice size: {slice_size}...')
		with parallel_backend("loky", inner_max_num_threads=2):
			batch_df_corr_records = Parallel(n_jobs=num_processes, verbose=10, pre_dispatch="all")(delayed(correlate_slice)(slice, use_full_matrix, min_cpgs) for slice in slices)
			
		df_corr_records = [item for sublist in batch_df_corr_records for item in sublist]
		df_corr = pd.DataFrame.from_records(df_corr_records)
		
		list_pearson_r = list(df_corr['pearson_r'])
		list_num_pairs = list(df_corr['num_pairs'])
		
		unique_n_pairs = list(set(list_num_pairs))
		
		zip_list = list(zip(list_pearson_r, list_num_pairs))
		
		print('getting p values...')
		num_lines = len(zip_list)
		slice_size = num_lines // num_processes
		slices = [zip_list[i:i+slice_size] for i in range(0, num_lines, slice_size)]
		
		with parallel_backend("loky", inner_max_num_threads=2):
			all_p_values = Parallel(n_jobs=num_processes, verbose=10, pre_dispatch="all")(delayed(get_p_values)(slice) for slice in slices)
		
		df_corr['pearson_p'] = [item for sublist in all_p_values for item in sublist]
		
		# output batch file
		random_string = generate_random_string(10)
		path_output_csv_batch = path_output_csv.replace('.csv', f'.{random_string}.csv')
		df_corr.to_csv(path_output_csv_batch, index=False)
		batch_csv_paths.append(path_output_csv_batch)
		
		# since file handle iterator cannot go back, we use the last value to make sure nothing is missed
		# this is triggered once a new read_id triggers execution of a batch
		# new read dict is required for batch, read_id should not be in it
		read_dict = {}
		read_dict[read_id] = {'meth_values': [], 'positions': [], 'wgbs_values': [], 'read_id': read_id}
		read_dict[read_id]['meth_values'].append(meth_value)
		read_dict[read_id]['positions'].append(position)
		read_dict[read_id]['wgbs_values'].append(wgbs_value)
		
		del batch_df_corr_records
		del df_corr_records
		
		i_batch += 1
		
print('Finished processing, stitching batch CSVs together...')
# stitch batch csvs together
df_corrs = [pd.read_csv(path) for path in batch_csv_paths]
df = pd.concat(df_corrs)
df.to_csv(path_output_csv, index=False)

# delete batch csvs
for path in batch_csv_paths:
	os.remove(path)

print('Done')