import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.animation as animation


## dt INSTABILITY


EI = 90666667.0 	# [Nm^2]
density = 2400.0 	# kg/m^3
area = 0.4 			# m^2
span = 20.0 		# m
velocity = 2.0 		# m/s
excitationModes = 15
Q = 20000.0 		# N
spring = 500000.0 	# N/m
lowestFreq = np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )
dt = 0.01/lowestFreq
dx = 0.1
g = 9.81			# m/s^2
mass = 2000.0 		# kg
damp = 2*np.sqrt(spring*mass)-1		# N/m/s (just under critical damping)
print(damp)



def phi_m(m, x, L):
	phi = np.sin(m*np.pi*x / L)
	return(phi)

x = np.arange(0, span, dx)
time = np.arange(0, span/velocity, dt)

#intial conditions:
w = np.zeros((len(time), len(x)))
wDiff = np.zeros(len(time))
wQ = np.zeros(len(time))
u = np.zeros(len(time))
u[0] = 2.0*mass*g/spring
uDiff = np.zeros(len(time))

for t in np.arange(2, len(time)):
	xQ = min( len(x)-1, int(round((velocity*time[t])/dx)) )
	wDiff[t] = (w[t-1, xQ] - w[t-2, xQ])/dt
	wQ[t] = w[t-1, xQ]
	uDiff[t] = (u[t-1] - u[t-2])/dt

	for m in np.arange(1, excitationModes+1):

		eigenfreq = m**2*np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )

		Q_term = spring*( -wQ[1:t+1] + u[0:t] ) + damp*( -wDiff[1:t+1] + uDiff[1:t+1] )

		duhamel = np.sin( eigenfreq*(time[t]-time[1:t+1]) )
		integral_w = np.trapz( Q_term*phi_m(m, velocity*time[1:t+1], span)*duhamel, time[1:t+1] )	
		integral_u = np.trapz( Q_term*(time[t]-time[1:t+1]), time[1:t+1] )	
		w[t] += phi_m(m,x,span) * (2/(density*area*span*eigenfreq)) * integral_w
		u[t] = mass*g/spring + g*time[t]**2/2 - 1/mass * integral_u

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
ax.set_ylim([-0.05, 0.1])
wline, = ax.plot(x, -wAna, color='r')
beam, = ax.plot(x, w[0], color='k')
# forceLine = ax.axvline(0, color='b', linestyle='-')
massLoc, = ax.plot(0.0, -w[0,0]+u[0], 'o', color='r')

def animationFrames(i):
	beam.set_ydata( -w[i] )
	# forceLine.set_xdata(velocity*time[i])
	wline.set_ydata( -analyticalW(Q, x, span, EI, velocity*time[i]) )
	massLoc.set_xdata(velocity*time[i])
	xQ = min( len(x)-1, int(round((velocity*time[i])/dx)) )
	massLoc.set_ydata(-w[i, xQ]+u[i])
	return(beam, massLoc)


	
ani = animation.FuncAnimation(fig, animationFrames, np.arange(0, len(w)), interval=3)



plt.show()
	

	