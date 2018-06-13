import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.animation as animation


## TWO SPRING DAMPER MASSES AT CONSTANT SPEED.

def phi_m(m, x, L):
	phi = np.sin(m*np.pi*x / L)
	return(phi)




# input vehicle
mass = 80.0 		# kg
velocity = 2.0 		# m/s
spring = 5000.0 	# N/m
damp = 10000.0 		# N/m/s
wbase = 2.0 		# m
inertia = None
Q = 20000.0 	 	# N

# input beam
span = 20.0 		# m
area = 0.4 			# m^2
EI = 90666667.0 	# Nm^2
density = 2400.0 	# kg/m^3
excitationModes = 5

g = 9.81 	# gravitational constant
massOffset = mass*g/spring + 0.02 	# m




trajectory = span + 2*wbase 		# includes beam and land abutment on both ends

# time steps
lowestFreq = np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )
dt = 0.01/lowestFreq
time = np.arange(0, (span+wbase)/velocity+dt, dt)
t1 = int((wbase/velocity)/dt)	# Q1 enters the beam (as time index)
t2 = int((span/velocity)/dt)	# Q2 leaves the beam (as time index)

# x steps
dx = 0.1
x = np.arange(0, trajectory+dx, dx)	
x1 = int(wbase/dx)			# beginning of beam (as x index)
x2 = int((wbase+span)/dx)	# end of beam (as x index)

# initializing empty arrays for the numerical calculation
w = np.zeros((len(time), len(x)))		# w(t)
wDiff = np.zeros((len(time), len(x)))	# w'(t)			

u1 = np.zeros(len(time))			# u1(t)
u1[0] = mass*g/spring 				# initial condition for u1(0)
u1Diff = np.zeros(len(time)) 		# u1'(t)
u2 = np.zeros(len(time))			# u2(t)
u2[0] = mass*g/spring 				# initial condition for u2(0)
u2Diff = np.zeros(len(time)) 		# u2'(t)

Q1 = np.zeros(len(time))
Q2 = np.zeros(len(time))
Q1[:2] = mass*g
Q2[:2] = mass*g

xQ1 = np.zeros(len(time))
xQ2 = np.zeros(len(time))

for t in np.arange(1, len(time)):

	# determine x location (step) of forces at time t
	xQ1[t] = min( x2, int(round((velocity*time[t])/dx)) )			# x step location of force 1
	xQ2[t] = min( len(x)-1, int(round((velocity*time[t])/dx))+x1 )	# x step location of force 2


	# determine displacements u1(t) and u2(t)
	int_u1 = np.trapz( (time[t]-time[0:t+1])*Q1[0:t+1], time[0:t+1] )		
	int_u2 = np.trapz( (time[t]-time[0:t+1])*Q2[0:t+1], time[0:t+1] )

	u1[t] = mass*g/spring + g*time[t]**2/2 - 1/mass*int_u1		
	u2[t] = mass*g/spring + g*time[t]**2/2 - 1/mass*int_u2

	u1Diff[t] = (u1[t]-u1[t-1])/dt
	u2Diff[t] = (u2[t]-u2[t-1])/dt

	# determine beam displacement for m modes of excitation
	for m in np.arange(1, excitationModes+1):

		eigenfreq = m**2*np.pi**2 / (span**2) * np.sqrt( EI/(density*area) ) 	# CHECK IF THIS EQ IS CORRECT. YES IT IS.
	
		phi1 = phi_m(m, velocity*time[t]-wbase, span)
		phi2 = phi_m(m, velocity*time[t], span)

		# accounting for heaviside stepfunction for Q1
		if t>=t1:  
			duhamel1 = np.sin( eigenfreq*(time[t]-time[t1:t+1]) )
			int_w1 = np.trapz( phi1*Q1[t1:t+1]*duhamel1, time[t1:t+1] )
		else:
			int_w1 = 0

		# accounting for heaviside stepfunction for Q2
		if t<t2:   
			duhamel2 = np.sin( eigenfreq*(time[t]-time[0:t+1]) )
			int_w2 = np.trapz( phi2*Q2[0:t+1]*duhamel2, time[0:t+1] )	
		else:
			duhamel2 = np.sin( eigenfreq*(time[t]-time[0:t2]) )	
			int_w2 = np.trapz( phi2*Q2[0:t2]*duhamel2, time[0:t2] )
			

		phi = phi_m(m, np.arange(0, span, dx),span)
		w[t,x1:x2] += w[t,x1:x2] + phi*(2/(density*area*span*eigenfreq)*(int_w1+int_w2))

	wDiff[t] = ( w[t]-w[t-1] )/dt

	# determine the forces for next timestep
	if t<(len(time)-2):
		Q1[t+1] = spring*( u1[t] - w[t, xQ1[t]] ) + damp*( u1Diff[t] - wDiff[t, xQ1[t]] )
		Q2[t+1] = spring*( u2[t] - w[t, xQ2[t]] ) + damp*( u2Diff[t] - wDiff[t, xQ2[t]] )


# Analytical w line
def analyticalW(F, x, L, EI, loc):
	W = np.zeros(len(x))
	bound = int(math.ceil(len(x) * loc/L))
	b = L-loc
	W[:bound+1] = F*b*x[:bound+1]*(L**2-b**2-(x[:bound+1])**2)/(6.0*L*EI)
	W[bound+1:] = F*b*x[bound+1:]*(L**2-b**2-(x[bound+1:])**2)/(6.0*L*EI) + F*(x[bound+1:]-loc)**3/(6.0*EI)
	return(W)



fig, ax = plt.subplots()
ax.set_ylim([0.05, -0.1])
# wAna = analyticalW(Q, x, testBeam.span, testBeam.EI, x[0])
# wline, = ax.plot(x, -wAna, color='r')
beam, = ax.plot(x, w[0], color='k')
# forceLine = ax.axvline(0, color='b', linestyle='-')
Q1Loc, = ax.plot(0.0, u1[0]-massOffset, 'o', color='r')
Q2Loc, = ax.plot(0.0, u2[0]-massOffset, 'o', color='g')


# Generate frames for animation
def animationFrames(i):
	beam.set_ydata( w[i] )

	Q1Loc.set_xdata(velocity*time[i])
	Q1Loc.set_ydata(u1[i]-massOffset)

	Q2Loc.set_xdata(velocity*time[i]+wbase)
	Q2Loc.set_ydata(u2[i]-massOffset)

	return(beam, Q1Loc, Q2Loc)


ani = animation.FuncAnimation(fig, animationFrames, np.arange(0, len(w)), interval=3)
plt.show()
	

