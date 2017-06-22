#!/usr/bin/python
import sys
import numpy as np
import scipy.io
import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt


def plot_results(state,dstate):

	# make automatic for higher number of derivatives
	states = state[:, 1:8]
	dstates = dstate[:, 1:8]

	# Plot Positions
	fig = plt.figure(1)
	plt.subplot(211)
	plt.plot(state[:,0], states[:,0] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 0')
	plt.plot(state[:,0],states[:,1] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 1')
	plt.plot(state[:,0],states[:,2] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 2')
	plt.plot(state[:,0],states[:,3] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 3')
	plt.plot(state[:,0],states[:,4] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 4')
	plt.plot(state[:,0],states[:,5] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 5')
	plt.plot(state[:,0],states[:,6] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 6')	
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Joint Angles (deg)')
	
	plt.subplot(212)
	plt.plot(state[:,0],dstates[:,0] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 0')
	plt.plot(state[:,0],dstates[:,1] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 1')
	plt.plot(state[:,0],dstates[:,2] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 2')
	plt.plot(state[:,0],dstates[:,3] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 3')
	plt.plot(state[:,0],dstates[:,4] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 4')
	plt.plot(state[:,0],dstates[:,5] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 5')
	plt.plot(state[:,0],dstates[:,6] / np.pi * 180., linewidth=2.0, linestyle='-', marker='', label='joint 6')	
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Joint Velocities (deg/s)')
	plt.xlabel('Time')
	
if __name__ == '__main__':
	
	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)
	
	state = []
	state = np.loadtxt(str(sys.argv[1]))
	
	dstate = []
	dstate = np.loadtxt(str(sys.argv[2]))

	plot_results(state,dstate)
	plt.show(True)
