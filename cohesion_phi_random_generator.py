#To run 
#python E_random_generator 10

import sys
import numpy as np

N = int(sys.argv[1])#Get the first command line argument

while (N > 0):
	n = np.random.randn(1)[0] #Genrate a single random number
	cohesion = 25000 + n*2000 #Scale to correct mean and s.dev
	n = np.random.randn(1)[0]
	phi = 0.35 + n*0.05 #Scale to correct mean and s.dev
	if(cohesion >= 0 and phi >= 0):
		print(cohesion, phi)
		N = N - 1
	
	
