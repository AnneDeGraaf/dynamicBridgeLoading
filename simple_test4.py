import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.animation as animation


## CALCULATES SITUATION FOR SINGLE CONSTANT MOVING POINT-LOAD WITH DAMPING AND SPRING


EI = 90666667.0 	# [Nm^2]
density = 2400.0 	# kg/m^3
area = 0.4 			# m^2
span = 20.0 		# m
velocity = 2.0 		# m/s
excitationModes = 15
Q = 20000.0 		# N
damp = 10000.0		# N/m/s
spring = 500000.0 	# N/m
lowestFreq = np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )
dt = 0.05/lowestFreq
dx = 0.1



def phi_m(m, x, L):
	phi = np.sin(m*np.pi*x / L)
	return(phi)

x = np.arange(0, span, dx)
time = np.arange(0, span/velocity, dt)

#intial conditions:
w = np.zeros((len(time), len(x)))
wDiff = np.zeros(len(time))
u = np.zeros(len(time))


for t in np.arange(2, len(time)):
	xQ = min( len(x)-1, int(round((velocity*time[t])/dx)) )
	wDiff[t] = (w[t-1, xQ] - w[t-2, xQ])/(dt)
	u[t] = -w[t-1, xQ] + 0.04

	for m in np.arange(1, excitationModes+1):

		eigenfreq = m**2*np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )

		duhamel = np.sin( eigenfreq*(time[t]-time[0:t+1]) )
		integral = np.trapz( (u[0:t+1]*spring - wDiff[0:t+1]*damp)*phi_m(m, velocity*time[0:t+1], span)*duhamel, time[0:t+1] )	
		w[t] += phi_m(m,x,span) * (2/(density*area*span*eigenfreq)) * integral


# Analytical w line
def analyticalW(F, x, L, EI, loc):
	W = np.zeros(len(x))
	bound = int(math.ceil(len(x) * loc/L))
	b = L-loc
	W[:bound+1] = F*b*x[:bound+1]*(L**2-b**2-(x[:bound+1])**2)/(6.0*L*EI)
	W[bound+1:] = F*b*x[bound+1:]*(L**2-b**2-(x[bound+1:])**2)/(6.0*L*EI) + F*(x[bound+1:]-loc)**3/(6.0*EI)
	return(W)

wAna = analyticalW(Q,x,span,EI, x[0])

fig, ax = plt.subplots()
ax.set_ylim([-0.05, 0.01])
wline, = ax.plot(x, -wAna, color='r')
beam, = ax.plot(x, w[0], color='k')
forceLine = ax.axvline(0, color='b', linestyle='-')


def animationFrames(i):
	beam.set_ydata( -w[i] )
	forceLine.set_xdata(velocity*time[i])
	wline.set_ydata( -analyticalW(Q, x, span, EI, velocity*time[i]) )
	return(beam, forceLine)


	
ani = animation.FuncAnimation(fig, animationFrames, np.arange(0, len(w)), interval=33)



plt.show()
	

	