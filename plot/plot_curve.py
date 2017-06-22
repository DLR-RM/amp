#!/usr/bin/python
import sys
import numpy as np
import scipy.io
import matplotlib as mpl
import matplotlib.pyplot as plt
import math


def plot_results(fitted):

	# make automatic for higher number of derivatives
	timef = fitted[:, 0]
	posf = np.zeros(fitted.shape) 
	posf[:,0] =  fitted[:, 1+(0*3)]
	posf[:,1] =  fitted[:, 1+(1*3)]
	posf[:,2] =  fitted[:, 1+(2*3)]
	posf[:,3] =  fitted[:, 1+(3*3)]
	posf[:,4] =  fitted[:, 1+(4*3)]
	posf[:,5] =  fitted[:, 1+(5*3)]
	posf[:,6] =  fitted[:, 1+(6*3)]
	
	velf = np.zeros(fitted.shape) 
	velf[:,0] =  fitted[:, 2+(0*3)]
	velf[:,1] =  fitted[:, 2+(1*3)]
	velf[:,2] =  fitted[:, 2+(2*3)]
	velf[:,3] =  fitted[:, 2+(3*3)]
	velf[:,4] =  fitted[:, 2+(4*3)]
	velf[:,5] =  fitted[:, 2+(5*3)]
	velf[:,6] =  fitted[:, 2+(6*3)]
	
	accf = np.zeros(fitted.shape) 
	accf[:,0] =  fitted[:, 3+(0*3)]
	accf[:,1] =  fitted[:, 3+(1*3)]
	accf[:,2] =  fitted[:, 3+(2*3)]
	accf[:,3] =  fitted[:, 3+(3*3)]
	accf[:,4] =  fitted[:, 3+(4*3)]
	accf[:,5] =  fitted[:, 3+(5*3)]
	accf[:,6] =  fitted[:, 3+(6*3)]

	

	# Plot Positions
	fig = plt.figure(1)
	plt.subplot(311)
	pos_x = plt.plot(timef, posf[:,0])
	pos_y = plt.plot(timef, posf[:,1])
	pos_z = plt.plot(timef, posf[:,2])
	plt.setp(pos_x, linewidth=2.0, linestyle='-', color='r')
	plt.setp(pos_y, linewidth=2.0, linestyle='-', color='g')
	plt.setp(pos_z, linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('position m')
	
	
	plt.subplot(312)
	vel_x = plt.plot(timef, velf[:,0])
	vel_y = plt.plot(timef, velf[:,1])
	vel_z = plt.plot(timef, velf[:,2])
	plt.setp(vel_x, linewidth=2.0, linestyle='-', color='r')
	plt.setp(vel_y, linewidth=2.0, linestyle='-', color='g')
	plt.setp(vel_z, linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('velocity m/s')
	
	plt.subplot(313)
	acc_x = plt.plot(timef, accf[:,0])
	acc_y = plt.plot(timef, accf[:,1])
	acc_z = plt.plot(timef, accf[:,2])
	plt.setp(acc_x, linewidth=2.0, linestyle='-', color='r')
	plt.setp(acc_y, linewidth=2.0, linestyle='-', color='g')
	plt.setp(acc_z, linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('acceleration m/s^2')
	plt.xlabel('time')
    
if __name__ == '__main__':
	
	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)
	
	data = []
	data = np.loadtxt(str(sys.argv[1]))

	plot_results(data)
	plt.show(True)


	
