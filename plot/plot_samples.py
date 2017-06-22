#!/usr/bin/python
import sys
import numpy as np
import scipy.io
import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os.path

def plot_results():
	#
	# Get the data sizes
	#
	n_nodes=0
	samples = np.loadtxt("tree_samples_0.dat")
	if  (samples.ndim > 1):
		n_nodes=samples.shape[0]
	elif  (samples.size > 1):
		n_nodes=1
	

	n_failed=0
	if (os.path.isfile("failed_samples_0.dat")): 
		failed_samples = np.loadtxt("failed_samples_0.dat")
		if (failed_samples.ndim > 1):
			n_failed=failed_samples.shape[0]
		elif (failed_samples.size > 1):
			n_failed=1
	
	#
	# Load edge files and plot in 3d
	#
	fig = plt.figure(1)
	ax = fig.add_subplot(111,projection='3d')
	for x in range(0, n_failed):
		new_edge = np.loadtxt("failed_tcps_"+str(x)+".dat")
		# edges = np.append(edges, new_edge, axis=0)
		ax.plot(new_edge[:,1], new_edge[:,2], new_edge[:,3], c='r',
				linewidth='1.0', linestyle='-')
	
	for x in range(0, n_nodes-1):
		new_edge = np.loadtxt("tree_tcps_"+str(x)+".dat")
		# edges = np.append(edges, new_edge, axis=0)
		ax.plot(new_edge[:,1], new_edge[:,2], new_edge[:,3], c='b',
				linewidth='2.0', linestyle='-')
	
	#
	# Plot nodes as scatter
	#
	if  (n_failed > 1):
		ax.scatter(failed_samples[:,1], failed_samples[:,2], failed_samples[:,3], c='r', marker='o')
	elif (n_failed > 0):
		ax.scatter(failed_samples[1], failed_samples[2], failed_samples[3], c='r', marker='o')	

	if (n_nodes > 1):
		ax.scatter(samples[:,1], samples[:,2], samples[:,3], c='b', marker='^')
	elif (n_nodes > 0):
		ax.scatter(samples[1], samples[2], samples[3], c='b', marker='^')
		
	ax.grid(True)
	ax.autoscale(True)
	ax.axis('tight')
	ax.set_xlabel('X (m)')
	ax.set_ylabel('Y (m)')
	ax.set_zlabel('Z (m)')


	#
	# Plot Quaternions as scatter
	#
	fig = plt.figure(2)
	ax = fig.add_subplot(111,projection='3d')
	if  (n_failed > 1):
		ax.scatter(failed_samples[:,4], failed_samples[:,5], failed_samples[:,6], c='r', marker='o')
	elif (n_failed > 0):
		ax.scatter(failed_samples[4], failed_samples[5], failed_samples[6], c='r', marker='o')	

	if (n_nodes > 1):
		ax.scatter(samples[:,4], samples[:,5], samples[:,6], c='b', marker='^')
	elif (n_nodes > 0):
		ax.scatter(samples[4], samples[5], samples[6], c='b', marker='^')
	ax.grid(True)
	ax.autoscale(True)
	ax.axis('tight')
	ax.set_xlabel('Qx')
	ax.set_ylabel('Qy')
	ax.set_zlabel('Qz')
	


if __name__ == '__main__':

	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)

	plot_results()
	plt.show(True)
