#!/usr/bin/python
from __future__ import division
import sys
import numpy as np
import scipy.io
import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os.path


def plot_results(log0):

	#
	# Count the succeeded simulations 
	#
	log0_suc = log0 
	avg_times_0 = 0
	avg_edges_0 = 0
	avg_iters_0 = 0
	avg_tries_0 = 0
	avg_times_suc_0 = 0
	avg_edges_suc_0 = 0
	avg_iters_suc_0 = 0
	avg_tries_suc_0 = 0
	num_0 = log0[log0.shape[0]-1,0] 
	num_suc_0 = log0.shape[0] 
	
	x = 0
	xs = 0
	while (x < log0.shape[0]):
		if(log0[x,3] == 0):
			num_suc_0 -= 1
			log0_suc = np.delete(log0_suc,xs,0)
		else:
			avg_times_suc_0 += log0[x,4]/1000
			avg_edges_suc_0 += log0[x,2]
			avg_iters_suc_0 += log0[x,1]
			if (log0.shape[1]>=6):
				avg_tries_suc_0 += log0[x,5]
			
			xs += 1

		avg_times_0 += log0[x,4]/1000
		avg_edges_0 += log0[x,2]
		avg_iters_0 += log0[x,1]
		if (log0.shape[1]>=6):
			avg_tries_0 += log0[x,5]
		x += 1

	
	print num_suc_0
	print num_0
	print 'With Manipulability:'
	print 'No. Successful: {0:.4f}'.format(num_suc_0/num_0)
	print 'Avg. No. Edges:  ', str(avg_edges_suc_0/num_suc_0) 
	print 'Avg. No. Iters:  ', str(avg_iters_suc_0/num_suc_0) 
	print 'Avg. No. Tries:  ', str(avg_tries_suc_0/num_suc_0) 
	print 'Avg. Comp. Time: ', str(avg_times_suc_0/num_suc_0)

	#
	# Plot succeeded scatter
	#
	fig = plt.figure(2)
	plt.subplot(311)
	plt.scatter(log0_suc[:,0], log0_suc[:,4]/1000, c='r', marker='.')	
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Comp. Time (s)')
	
	plt.subplot(312)
	plt.scatter(log0_suc[:,0], log0_suc[:,2], c='r', marker='.')	
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Edges')
	
	plt.subplot(313)
	plt.scatter(log0_suc[:,0], log0_suc[:,1], c='r', marker='.')	
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.xlabel('Simulations')
	plt.ylabel('Iterations')

	#
	# Histogram
	#
	fig = plt.figure(3)
	ax = fig.add_subplot(311)
	ax.hist(log0_suc[:,4]/1000, 50)
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Simulations')
	plt.xlabel('Comp. Time (s)')
	
	ax = fig.add_subplot(312)
	ax.hist(log0_suc[:,2], 50)
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Simulations')
	plt.xlabel('Edges')
	
	ax = fig.add_subplot(313)
	ax.hist(log0_suc[:,1], 50)
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Simulations')
	plt.xlabel('Iterations')


if __name__ == '__main__':

	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)
	
	log0 = []
	log0 = np.loadtxt(str(sys.argv[1]))
	
	plot_results(log0)
	plt.show(True)

