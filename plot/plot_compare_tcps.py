#!/usr/bin/python
import sys
import numpy as np
import scipy.io
import matplotlib as mpl
mpl.use('GTKAgg')
import matplotlib.pyplot as plt
import math


def plot_results(tcp,dtcp,tcp_a,dtcp_a):

	# make automatic for higher number of derivatives
	timef = tcp[:, 0]
	posf = np.zeros(tcp.shape) 
	posf[:,0] =  tcp[:, 1]
	posf[:,1] =  tcp[:, 2]
	posf[:,2] =  tcp[:, 3]
	posf[:,3] =  tcp[:, 4]
	posf[:,4] =  tcp[:, 5]
	posf[:,5] =  tcp[:, 6]
	posf[:,6] =  tcp[:, 7]
	
	velf = np.zeros(dtcp.shape) 
	velf[:,0] =  dtcp[:, 1]
	velf[:,1] =  dtcp[:, 2]
	velf[:,2] =  dtcp[:, 3]
	velf[:,3] =  dtcp[:, 4]
	velf[:,4] =  dtcp[:, 5]
	velf[:,5] =  dtcp[:, 6]
	
	# make automatic for higher number of derivatives
	time_curve = tcp_a[:, 0]
	pos = np.zeros(posf.shape) 
	pos[:,0] =  tcp_a[:, 1]
	pos[:,1] =  tcp_a[:, 2]
	pos[:,2] =  tcp_a[:, 3]
	pos[:,3] =  tcp_a[:, 4]
	pos[:,4] =  tcp_a[:, 5]
	pos[:,5] =  tcp_a[:, 6]
	pos[:,6] =  tcp_a[:, 7]
	
	vel = np.zeros(velf.shape) 
	vel[:,0] =  dtcp_a[:, 1]
	vel[:,1] =  dtcp_a[:, 2]
	vel[:,2] =  dtcp_a[:, 3]
	vel[:,3] =  dtcp_a[:, 4]
	vel[:,4] =  dtcp_a[:, 5]
	vel[:,5] =  dtcp_a[:, 6]

	#############################################################################
	# Plot Positions & Velocities
	fig = plt.figure(1)
	plt.subplot(221)
	posf_x = plt.plot(timef, posf[:,0], label='numerical')
	posf_y = plt.plot(timef, posf[:,1])
	posf_z = plt.plot(timef, posf[:,2])
	pos_x = plt.plot(time_curve, pos[:,0], label='desired')
	pos_y = plt.plot(time_curve, pos[:,1])
	pos_z = plt.plot(time_curve, pos[:,2])
	plt.setp(posf_x, linewidth=2.0, linestyle='-', color='r')
	plt.setp(posf_y, linewidth=2.0, linestyle='-', color='g')
	plt.setp(posf_z, linewidth=2.0, linestyle='-', color='b')
	plt.setp(pos_x, linewidth=1.0, linestyle='--', color='r')
	plt.setp(pos_y, linewidth=1.0, linestyle='--', color='g')
	plt.setp(pos_z, linewidth=1.0, linestyle='--', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('position m')
	
	plt.subplot(223)
	velf_x = plt.plot(timef, velf[:,0], label='numerical')
	velf_y = plt.plot(timef, velf[:,1])
	velf_z = plt.plot(timef, velf[:,2])
	vel_x = plt.plot(time_curve, vel[:,0], label='desired')
	vel_y = plt.plot(time_curve, vel[:,1])
	vel_z = plt.plot(time_curve, vel[:,2])
	plt.setp(velf_x, linewidth=2.0, linestyle='-', color='r')
	plt.setp(velf_y, linewidth=2.0, linestyle='-', color='g')
	plt.setp(velf_z, linewidth=2.0, linestyle='-', color='b')
	plt.setp(vel_x, linewidth=1.0, linestyle='--', color='r')
	plt.setp(vel_y, linewidth=1.0, linestyle='--', color='g')
	plt.setp(vel_z, linewidth=1.0, linestyle='--', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('velocity m/s')
	
	 
	# Plot Errors
	plt.subplot(222)
	pos_dx = plt.plot(timef, pos[:,0] - posf[:,0])
	pos_dy = plt.plot(timef, pos[:,1] - posf[:,1])
	pos_dz = plt.plot(timef, pos[:,2] - posf[:,2])
	plt.setp(pos_dx, linewidth=2.0, linestyle='-', color='r')
	plt.setp(pos_dy, linewidth=2.0, linestyle='-', color='g')
	plt.setp(pos_dz, linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('position error m')
		
	plt.subplot(224)
	vel_dx = plt.plot(timef, vel[:,0] - velf[:,0])
	vel_dy = plt.plot(timef, vel[:,1] - velf[:,1])
	vel_dz = plt.plot(timef, vel[:,2] - velf[:,2])
	plt.setp(vel_dx, linewidth=2.0, linestyle='-', color='r')
	plt.setp(vel_dy, linewidth=2.0, linestyle='-', color='g')
	plt.setp(vel_dz, linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('velocity error m/s')
   


	#############################################################################
	# Plot Quaternions
	fig = plt.figure(2)
	plt.subplot(221)
	posf_qx = plt.plot(timef, posf[:,3], label='numerical')
	posf_qy = plt.plot(timef, posf[:,4])
	posf_qz = plt.plot(timef, posf[:,5])
	posf_qw = plt.plot(timef, posf[:,6])
	pos_qx = plt.plot(timef, pos[:,3], label='desired')
	pos_qy = plt.plot(timef, pos[:,4])
	pos_qz = plt.plot(timef, pos[:,5])
	pos_qw = plt.plot(timef, pos[:,6])
	plt.setp(posf_qx, linewidth=2.0, linestyle='-', color='r')
	plt.setp(posf_qy, linewidth=2.0, linestyle='-', color='g')
	plt.setp(posf_qz, linewidth=2.0, linestyle='-', color='b')
	plt.setp(posf_qw, linewidth=2.0, linestyle='-', color='k')
	plt.setp(pos_qx, linewidth=1.0, linestyle='--', color='r')
	plt.setp(pos_qy, linewidth=1.0, linestyle='--', color='g')
	plt.setp(pos_qz, linewidth=1.0, linestyle='--', color='b')
	plt.setp(pos_qw, linewidth=1.0, linestyle='--', color='k')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('quaternions')
	
	plt.subplot(223)
	velf_wx = plt.plot(timef, velf[:,3], label='numerical')
	velf_wy = plt.plot(timef, velf[:,4])
	velf_wz = plt.plot(timef, velf[:,5])
	vel_wx = plt.plot(timef, vel[:,3], label='desired')
	vel_wy = plt.plot(timef, vel[:,4])
	vel_wz = plt.plot(timef, vel[:,5])
	plt.setp(velf_wx, linewidth=2.0, linestyle='-', color='r')
	plt.setp(velf_wy, linewidth=2.0, linestyle='-', color='g')
	plt.setp(velf_wz, linewidth=2.0, linestyle='-', color='b')
	plt.setp(vel_wx, linewidth=1.0, linestyle='--', color='r')
	plt.setp(vel_wy, linewidth=1.0, linestyle='--', color='g')
	plt.setp(vel_wz, linewidth=1.0, linestyle='--', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('angular vel. rad/s')
	
 	# Plot Errors
	plt.subplot(222)
	pos_dqx = plt.plot(timef, pos[:,3] - posf[:,3])
	pos_dqy = plt.plot(timef, pos[:,4] - posf[:,4])
	pos_dqz = plt.plot(timef, pos[:,5] - posf[:,5])
	pos_dqw = plt.plot(timef, pos[:,6] - posf[:,6])
	plt.setp(pos_dqx, linewidth=2.0, linestyle='-', color='r')
	plt.setp(pos_dqy, linewidth=2.0, linestyle='-', color='g')
	plt.setp(pos_dqz, linewidth=2.0, linestyle='-', color='b')
	plt.setp(pos_dqw, linewidth=2.0, linestyle='-', color='k')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('orientation error')
	
	
	plt.subplot(224)
	vel_dwx = plt.plot(timef, vel[:,3] - velf[:,3])
	vel_dwy = plt.plot(timef, vel[:,4] - velf[:,4])
	vel_dwz = plt.plot(timef, vel[:,5] - velf[:,5])
	plt.setp(vel_dwx, linewidth=2.0, linestyle='-', color='r')
	plt.setp(vel_dwy, linewidth=2.0, linestyle='-', color='g')
	plt.setp(vel_dwz, linewidth=2.0, linestyle='-', color='b')
	plt.legend(loc='best')
	plt.grid(True)
	plt.autoscale(True)
	plt.axis('tight')
	plt.ylabel('angular vel. error rad/s')




if __name__ == '__main__':
	
	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)
	
	tcp = []
	tcp = np.loadtxt(str(sys.argv[1]))
	
	dtcp = []
	dtcp = np.loadtxt(str(sys.argv[2]))
	
	tcp_a = []
	tcp_a = np.loadtxt(str(sys.argv[3]))
	
	dtcp_a = []
	dtcp_a = np.loadtxt(str(sys.argv[4]))

	plot_results(tcp,dtcp,tcp_a,dtcp_a)
	plt.show(True)
