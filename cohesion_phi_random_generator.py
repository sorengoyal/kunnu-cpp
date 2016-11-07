#To run 
#python E_random_generator 10

import sys
import numpy as np

N = int(sys.argv[1])#Get the first command line argument
f = open('cohesion_phi_n1_n2.dat', 'w')

while (N > 0):
	n1 = np.random.randn(1)[0] #Genrate a single random number
	cohesion = 25000 + n1*2000 #Scale to correct mean and s.dev
	n2 = np.random.randn(1)[0]
	phi = 0.35 + n2*0.05 #Scale to correct mean and s.dev
	if(cohesion >= 0 and phi >= 0):
		print(cohesion, phi)
		f.write("%.3f, %.3f, %.3f, %.3f\n" % (cohesion, phi, n1, n2))
		N = N - 1
	
	
