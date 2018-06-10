import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.animation as animation



## CALCULATES SITUATION FOR SINGLE STATIONARY CONSTANT LOAD


EI = 90666667.0 	# [Nm^2]
density = 2400.0 	# kg/m^3
area = 0.4 			# m^2
span = 20.0 		# m
excitationModes = 15
forceLocation = span/2

lowestFreq = np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )
print(lowestFreq)

def phi_m(m, x, L):
	phi = np.sin(m*np.pi*x / L)
	return(phi)

x = np.arange(0, span, 0.1)
t = np.arange(0, 2*np.pi/lowestFreq, 0.05/lowestFreq)


#intial conditions:
w = np.zeros((len(t)+1, len(x)))
Q = 20000.0 	# N

# Analytical w line
def analyticalW(F, x, L, EI, loc):
	W = np.zeros(len(x))
	bound = int(math.ceil(len(x) * loc/L))
	b = L-loc
	W[:bound+1] = F*b*x[:bound+1]*(L**2-b**2-(x[:bound+1])**2)/(6.0*L*EI)
	W[bound+1:] = F*b*x[bound+1:]*(L**2-b**2-(x[bound+1:])**2)/(6.0*L*EI) + F*(x[bound+1:]-loc)**3/(6.0*EI)
	return(W)

wAna = analyticalW(Q,x,span,EI, forceLocation)



for iStep in range(len(t)):

	for m in np.arange(1, excitationModes+1):

		eigenfreq = m**2*np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )

		integral = Q * (1/eigenfreq) * (1-np.cos(eigenfreq*t[iStep])) * phi_m(m, forceLocation, span)
		w[iStep] += phi_m(m,x,span) * (2/(density*area*span*eigenfreq)) * integral




fig, ax = plt.subplots()
line, = ax.plot(x, w[0])
line2, = ax.plot(x, -wAna, color='r')
line3, = ax.plot(x, -2*wAna, color = 'r')
forceLine = ax.axvline(forceLocation, color='b', linestyle='-')


# maxW = max( abs(w) )

def animationFrames(i):
	line.set_ydata( -w[i] )
	return(line,)
	
ani = animation.FuncAnimation(fig, animationFrames, np.arange(0, len(w)), interval=33)

plt.show()

## CALCULATED BY HAND AND RESULTS SEEM CORRECT.