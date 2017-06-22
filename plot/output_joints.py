#!/usr/bin/python
import sys
import numpy as np
import scipy.io


if __name__ == '__main__':
	
	
	state = []
	state = np.loadtxt(str(sys.argv[1]))
	
	for t in range(0, state.shape[0]) :
		for c in range(1, state.shape[1]) :
			sys.stdout.write(str(state[t,c]))
			sys.stdout.write(' ')
		print
