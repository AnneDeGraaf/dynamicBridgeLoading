import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
import pickle


## HYPERLOOP POD RUNNING THROUGH HYPERLOOP TUBE


# input vehicle
mass = 8000.0		# kg	(including payload)
velocity = 277.8 	# m/s   (1000.0 km/h)   277.8
spring = 60000.0	# N/m
damp = 4100.0 		# N.s/m
wbase = 16.0		# m 	(simplification of actual hyperloop model)
inertia = 2890.0	# kg.m^2  

# input beam
span = 40.0		# m
area = 11.375		# m^2
EI = 1.1827e12		# N.m^2
density = 2400.0	# kg/m^3  (concrete)
excitationModes = 5

g = 9.81 		# gravitational constant

massOffset = 0.5*mass*g/spring + 0.02 	# m

trajectory = span + 2*wbase 		# includes beam and land abutment on both ends

# time steps
lowestFreq = np.pi**2 / (span**2) * np.sqrt( EI/(density*area) )
dt = 0.1/lowestFreq
time = np.arange(0, (span+wbase+100.0)/(velocity)+dt, dt)
t1 = int((wbase/velocity)/dt)	# Q1 enters the beam (as time index)
t2 = int((span/velocity)/dt)	# Q2 leaves the beam (as time index)
t3 = int(((span+wbase)/velocity)/dt) # Q1 leaves the beam (as time index)

# x steps
dx = 0.1
x = np.arange(0, trajectory+dx, dx)	
x1 = int(wbase/dx)			# beginning of beam (as x index)
x2 = int((wbase+span)/dx)		# end of beam (as x index)

# initializing empty arrays for the numerical calculation
w = np.zeros((len(time), len(x)))	# w(t,x)
wDiff = np.zeros((len(time), len(x)))	# w'(t,x)			

u = np.zeros(len(time))			# u(t)
u[0] = 0.5*mass*g/spring 		# initial condition for u(0)
uDiff = np.zeros(len(time)) 		# u'(t)

theta = np.zeros(len(time))		# theta(t)
theta[0] = 0  				# initial condition for theta(0)
thDiff = np.zeros(len(time))		# theta'(t)

Q1 = np.zeros(len(time))
Q2 = np.zeros(len(time))
Q1[:2] = 0.5*mass*g
Q2[:2] = 0.5*mass*g

xQ1 = np.zeros(len(time))
xQ2 = np.zeros(len(time))

# phi function
def phi_m(m, x, L):
	phi = np.sin(m*np.pi*x / L)
	return(phi)


for t in np.arange(1, len(time)):

	# determine x location (step) of forces at time t
	xQ1[t] = min( x2, int(round((velocity*time[t])/dx)) )			# x step location of force 1
	xQ2[t] = min( len(x)-1, int(round((velocity*time[t])/dx))+x1 )		# x step location of force 2


	# determine intergrals for u(t) and theta(t) (trapezoidal method)
	int_Q1 = np.trapz( (time[t]-time[0:t+1])*Q1[0:t+1], time[0:t+1] )		
	int_Q2 = np.trapz( (time[t]-time[0:t+1])*Q2[0:t+1], time[0:t+1] )

	# determine vehicle displacement u(t) and pitch theta(t) and their derivatives (backward Euler method)
	u[t] = 0.5*mass*g/spring + 0.5*g*time[t]**2 - 1/mass*(int_Q1 + int_Q2)	
	uDiff[t] = (u[t]-u[t-1])/dt

	theta[t] = 0.5*wbase/inertia*(int_Q2 - int_Q1)
	thDiff[t] = (theta[t]-theta[t-1])/dt


	# determine beam displacement for m modes of excitation
	for m in np.arange(1, excitationModes+1):

		eigenfreq = m**2*np.pi**2 / (span**2) * np.sqrt( EI/(density*area) ) 	
	
		phi1 = phi_m(m, velocity*time[t]-wbase, span)
		phi2 = phi_m(m, velocity*time[t], span)

		# accounting for heaviside stepfunction for Q1, integral using trapezoidal method
		if t>=t1:  
			if t<t3:
				duhamel1 = np.sin( eigenfreq*(time[t]-time[t1:t+1]) )
				int_w1 = np.trapz( phi1*Q1[t1:t+1]*duhamel1, time[t1:t+1] )
			else:
				duhamel1 = np.sin( eigenfreq*(time[t]-time[t1:t3]) )
				int_w1 = np.trapz( phi1*Q1[t1:t3]*duhamel1, time[t1:t3] )
		else:
			int_w1 = 0

		# accounting for heaviside stepfunction for Q2, integral using trapezoidal method
		if t<t2:   
			duhamel2 = np.sin( eigenfreq*(time[t]-time[0:t+1]) )
			int_w2 = np.trapz( phi2*Q2[0:t+1]*duhamel2, time[0:t+1] )	
		else:
			duhamel2 = np.sin( eigenfreq*(time[t]-time[0:t2]) )	
			int_w2 = np.trapz( phi2*Q2[0:t2]*duhamel2, time[0:t2] )
			
		# beam deflection at time step t
		phi = phi_m(m, np.arange(0, span, dx),span)
		w[t,x1:x2] += w[t,x1:x2] + phi*(2/(density*area*span*eigenfreq)*(int_w1+int_w2))

	# time derivative of beam deflection at step t (backward Euler method)
	wDiff[t] = ( w[t]-w[t-1] )/dt

	# determine the forces for next timestep
	if t<(len(time)-2):
		Q1[t+1] = spring*( u[t] + 0.5*wbase*theta[t] - w[t, xQ1[t]] ) + damp*( uDiff[t] + 0.5*wbase*thDiff[t] - wDiff[t, xQ1[t]] )
		Q2[t+1] = spring*( u[t] - 0.5*wbase*theta[t] - w[t, xQ2[t]] ) + damp*( uDiff[t] - 0.5*wbase*thDiff[t] - wDiff[t, xQ2[t]] )


# storing the calculated data:
with open('example_hyperloop.pkl', 'w') as f:  # Python 3: open(..., 'wb')
    pickle.dump([w, u, theta, Q1, Q2], f)



# ANIMATION:

# initializing beam
fig, ax = plt.subplots(1,1)
ax.set_ylim([6.0, -6.0])
ax.set_yticks(np.arange(6.0, -6.0, -2.0))
ax.set_xlim([wbase, wbase+span])
ax.set_xticks(np.arange(wbase, wbase+span+5.0, 5.0))

beam, = ax.plot(x, w[0]*1000.0, color='dodgerblue', label='hyperloop tube')

# initializing mass
Q1Loc, = ax.plot(0.0, 0.002, 'o', color='slategrey', label='location Q1')
Q2Loc, = ax.plot(wbase, 0.002, 'o', color='darkslategrey', label='location Q2')


# Generate frames for animation
def animationFrames(i):
	beam.set_ydata( w[i]*1000.0 )

	Q1Loc.set_xdata(velocity*time[i])
	# Q1Loc.set_ydata(u[i]*2.0-massOffset*4.0)

	Q2Loc.set_xdata(velocity*time[i]+wbase)
	# Q2Loc.set_ydata(u[i]*2.0-massOffset*4.0)

	return(beam, Q1Loc, Q2Loc)


# show animation
ani = animation.FuncAnimation(fig, animationFrames, np.arange(0, len(w)), interval=33)
plt.show()
