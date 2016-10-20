#To run 
#python phi_random_generator 10

import sys
import numpy as np

N = int(sys.argv[1])#Get the first command line argument

while (N > 0):
	n = np.random.randn(1)[0] #Genrate a single random number
	final_no = 0.35 + n*0.05 #Scale to correct mean and s.dev
	if(final_no >= 0):
		print(final_no) 
		N = N - 1
	
