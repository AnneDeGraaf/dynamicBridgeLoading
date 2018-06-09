import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

EI = 90666667.0 	# [Nm^2]
density = 2400.0 	# kg/m^3
area = 0.4 			# m^2
span = 20.0 		# m
excitationModes = 15

lowestFreq = np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )

# w = som:  2/rhoAL * 1/eigenfreq * phi_m
# 	* int_0^t ( Q*phi_m(vtau)*heaviside * sin(eigenfreq(t-tau)) )dtau

def phi_m(m, x, L):
	phi = np.sin(m*np.pi*x / L)
	return(phi)

x = np.arange(0, span, 0.1)
t = np.arange(0, 2*np.pi/lowestFreq, 0.05/lowestFreq)

#intial conditions:
w = np.zeros((len(t)+1, len(x)))
Q = 20000.0 	# N




for iStep in range(len(t)):

	for m in np.arange(1, excitationModes+1):

		eigenfreq = m**2*np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )

		integral = Q * (1/eigenfreq) * (1-np.cos(eigenfreq*t[iStep])) * phi_m(m, span/4, span)
		w[iStep] += phi_m(m,x,span) * (2/(density*area*span*eigenfreq)) * integral
	



fig, ax = plt.subplots()
line, = ax.plot(w[0])
# maxW = max( abs(w) )


def animationFrames(i):
	line.set_ydata( -w[i] )
	return(line,)
	
ani = animation.FuncAnimation(fig, animationFrames, np.arange(0, len(w)), interval=33)

plt.show()