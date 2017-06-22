#!/usr/bin/python
import sys
import numpy as np
import scipy.io
import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt
import math


def plot_results(tcp,dtcp,ddtcp):

	

	# Plot Positions & Velocities
	fig = plt.figure(1)
	plt.subplot(321)
	plt.plot(tcp[:,0], tcp[:,1], linewidth=2.0, linestyle='-', color='r')
	plt.plot(tcp[:,0],tcp[:,2], linewidth=2.0, linestyle='-', color='g')
	plt.plot(tcp[:,0],tcp[:,3], linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('position m')

	plt.subplot(323)
	plt.plot(tcp[:,0],dtcp[:,1], linewidth=2.0, linestyle='-', color='r')
	plt.plot(tcp[:,0],dtcp[:,2], linewidth=2.0, linestyle='-', color='g')
	plt.plot(tcp[:,0],dtcp[:,3], linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('velocity m/s')

	plt.subplot(325)
	plt.plot(tcp[:,0],ddtcp[:,1], linewidth=2.0, linestyle='-', color='r')
	plt.plot(tcp[:,0],ddtcp[:,2], linewidth=2.0, linestyle='-', color='g')
	plt.plot(tcp[:,0],ddtcp[:,3], linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('acceleration m/s^2')
	
	# Plot Rotations & Angular Velocities
	plt.subplot(322)
	plt.plot(tcp[:,0],tcp[:, 4], linewidth=2.0, linestyle='-', color='r')
	plt.plot(tcp[:,0],tcp[:, 5], linewidth=2.0, linestyle='-', color='g')
	plt.plot(tcp[:,0],tcp[:, 6], linewidth=2.0, linestyle='-', color='b')
	plt.plot(tcp[:,0],tcp[:, 7], linewidth=2.0, linestyle='-', color='k')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Quaternion')

	plt.subplot(324)
	# / np.pi * 180.
	plt.plot(tcp[:,0],dtcp[:, 4], linewidth=2.0, linestyle='-', color='r')
	plt.plot(tcp[:,0],dtcp[:, 5], linewidth=2.0, linestyle='-', color='g')
	plt.plot(tcp[:,0],dtcp[:, 6], linewidth=2.0, linestyle='-', color='b')
	plt.plot(tcp[:,0],dtcp[:, 7], linewidth=2.0, linestyle='-', color='k')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Quaternion Rate 1/s')
	
	plt.subplot(326)
	plt.plot(tcp[:,0],ddtcp[:, 4], linewidth=2.0, linestyle='-', color='r')
	plt.plot(tcp[:,0],ddtcp[:, 5], linewidth=2.0, linestyle='-', color='g')
	plt.plot(tcp[:,0],ddtcp[:, 6], linewidth=2.0, linestyle='-', color='b')
	plt.plot(tcp[:,0],ddtcp[:, 7], linewidth=2.0, linestyle='-', color='k')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Quaternion Acceleration 1/s^2')

if __name__ == '__main__':

	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)

	tcp = []
	tcp = np.loadtxt(str(sys.argv[1]))

	dtcp = []
	dtcp = np.loadtxt(str(sys.argv[2]))

	ddtcp = []
	ddtcp = np.loadtxt(str(sys.argv[3]))

	plot_results(tcp,dtcp,ddtcp)
	plt.show(True)
