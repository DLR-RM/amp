#!/usr/bin/python
import sys
import numpy as np
import scipy.io
import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt
import math


#
# what is lambda? (c++ lambda?)
# filt?
#
def plot_results(Mm):

	# make automatic for higher number of derivatives
	time = Mm[:, 0]
	mm = Mm[:, 1]
	

	# Plot Positions
	fig = plt.figure(1)
	m = plt.plot(time, mm)
	plt.setp(m, linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('Manipulability Measure')
	plt.xlabel('Time')
	
	
if __name__ == '__main__':
	
	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)
	
	data = []
	data = np.loadtxt(str(sys.argv[1]))
	
	plot_results(data)
	plt.show(True)

