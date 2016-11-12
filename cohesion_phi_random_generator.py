import sys
import numpy as np
import os

N = 10
f = open('xi_n.dat', 'w')
g = open('cohesion_phi.dat','w')
while (N > 0):
	n1 = np.random.randn(1)[0] #Genrate a single random number
	cohesion = 25000 + n1*2000 #Scale to correct mean and s.dev
	n2 = np.random.randn(1)[0]
	phi = 0.35 + n2*0.05 #Scale to correct mean and s.dev
	if(cohesion >= 0 and phi >= 0):
		command = './demo_wavetrial_pulseVM ' + str(cohesion) + ' ' + str(phi)
		os.system(command)
		g.write("%.3f %.3f" % (cohesion, phi))
		f.write("%.3f, %.3f\n" % (n1, n2))
		N = N - 1	
	
