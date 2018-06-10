import numpy as np 
import math
import matplotlib.pyplot as plt 
import matplotlib.animation as animation


def main():
	# Perform beamDeflection

	# add some robustness
	return(0)


def beamDeflection(vehicleObject, beamObject, timeStepSize=0.1, spaceStepSize=0.1, excitationModes=4):
	## VERHAALTJE

	# extracting parameters from input vehicleObject 
	mass = vehicleObject.mass 
	velocity = vehicleObject.velocity
	spring = vehicleObject.spring
	damping = vehicleObject.damping
	wheelbase = vehicleObject.wheelbase
	inertia = vehicleObject.inertia 

	# extracting parameters from input beamObject
	span = beamObject.span
	area = beamObject.area
	EI = beamObject.EI
	density = beamObject.density

	g = 9.81 	# gravitational constant

	# setting up numerical time steps
	t0 = 0								# time of first wheel entering the beam
	t1 = wheelbase/velocity  			# time of second wheel entering the beam
	t2 = span/velocity					# time of first wheel leaving the beam
	t3 = (wheelbase+span)/velocity		# time of second wheel leaving the beam
	tMax = t3+5							# time of end of calculations

	t1Step = int(math.ceil(t1/timeStepSize)) + 1		# step number belonging to t1
	t2Step = int(math.ceil(t2/timeStepSize)) + 1		# step number belonging to t2
	t3Step = int(math.ceil(t3/timeStepSize)) + 1 		# step number belonging to t3
	ntSteps = int(math.ceil(tMax/timeStepSize)) + 1		# total number of time steps

	# setting up numerical space steps
	xL = span 												# span of beam
	xBeam = int(math.ceil(xL/spaceStepSize)) + 1 			# number of x steps along beam
	xPadding = int(math.ceil(wheelbase/spaceStepSize))  	# some padding before the beam for when vehicle is only partially on beam
	nxSteps = int(tMax*velocity/spaceStepSize) + xPadding	# total number of x steps during calculation

	# initializing empty arrays for the numerical calculation
	time = np.linspace(0, tMax, ntSteps)				# t
	vehicleDisplacement = np.zeros(ntSteps)				# u(t)
	vehiclePitch = np.zeros(ntSteps)					# theta(t)
	beamDisplacement = np.zeros((ntSteps, nxSteps))		# w(t): padding is added for the land abutment before the beam
	print(np.shape(beamDisplacement))
	force1 = np.zeros(ntSteps)							# Q_1(t)
	force2 = np.zeros(ntSteps)							# Q_2(t)

	# initial conditions
	force1[0] = mass*g/2
	force2[0] = mass*g/2


	for t in range(ntSteps):

		# determine x location (step) of forces at time t
		xF1 = int(round((velocity*time[t])/spaceStepSize))		
		xF2 = xF1 + xPadding


		# determine vehicle displacement u(t) and pitch theta(t)
		integral1 = np.trapz( time[t]*time[0:t+1]*force1[0:t+1], time[0:t+1] )		
		integral2 = np.trapz( time[t]*time[0:t+1]*force2[0:t+1], time[0:t+1] )

		vehicleDisplacement[t] = mass*g/spring + g*time[t]**2/2 - 1/mass*(integral1 + integral2)				
		vehiclePitch[t] = wheelbase/(2*inertia) * (integral2 - integral1)


		# determine beam displacement for m modes of excitation
		for m in range(excitationModes+1)[1:excitationModes+1]:

			eigenfreq = m**2*np.pi**2 / (span**2) * np.sqrt( EI/(density*area) ) 
		
			phi1 = np.sin(m*np.pi*(velocity*time[t]-wheelbase)/span)
			phi2 = np.sin(m*np.pi*(velocity*time[t])/span)

			# accounting for heaviside stepfunction for force1
			if t>=t1Step and t<t3Step:  
				duhamel1 = np.sin( eigenfreq*(time[t]-time[t1Step:t+1]) )
				integral3 = np.trapz( phi1*force1[t1Step:t+1]*duhamel1, time[t1Step:t+1] )
			elif t>=t1Step and t>t3Step:
				duhamel1 = np.sin( eigenfreq*(time[t]-time[t1Step:t3Step]) )	
				integral3 = np.trapz( phi1*force1[t1Step:t3Step]*duhamel1, time[t1Step:t3Step] )
			else:
				integral3 = 0

			# accounting for heaviside stepfunction for force2
			if t<t2Step:   
				duhamel2 = np.sin( eigenfreq*(time[t]-time[0:t+1]) )
				integral4 = np.trapz( phi2*force2[0:t+1]*duhamel2, time[0:t+1] )	
			else:
				duhamel2 = np.sin( eigenfreq*(time[t]-time[0:t2Step]) )	
				integral4 = np.trapz( phi2*force2[0:t2Step]*duhamel2, time[0:t2Step] )
				

			modeM = 2/(density*area*span*eigenfreq)*(integral3+integral4)

			for x in range(xBeam):														# MOET DIT IN N LOOP?
				phiX = np.sin(m*np.pi*x/span)
				beamDisplacement[t,x+xPadding] = beamDisplacement[t,x+xPadding] + phiX*modeM

		# determine the forces for next timestep
		if t<(ntSteps-1):
			force1[t+1] = spring*( vehicleDisplacement[t] - beamDisplacement[t, xF1] + wheelbase/2*vehiclePitch[t] )  
				# + damping*(  derivVehicleDisplacement - derivBeamDisplacement + wheelbase/2*derivVehiclePitch  ))  
			force2[t+1] = spring*( vehicleDisplacement[t] - beamDisplacement[t, xF2] + wheelbase/2*vehiclePitch[t] )  
				# + damping*(  derivVehicleDisplacement - derivBeamDisplacement - wheelbase/2*derivVehiclePitch  ))

	return(beamDisplacement)



class Vehicle:
	def __init__(self, mass, inertia, velocity, springStiffness, dampingCoefficient, wheelbase):
		self.mass = mass					## ADD SOMETHING ABOUT UNITS
		self.velocity = velocity
		self.spring = springStiffness
		self.damping = dampingCoefficient
		self.wheelbase = wheelbase
		self.inertia = inertia


class Beam:
	def __init__(self, span, crossSectionalArea, EI, materialDensity):
		self.span = span					## ADD SOMETHING ABOUT UNITS
		self.area = crossSectionalArea
		self.EI = EI
		self.density = materialDensity


testVehicle = Vehicle(300, 3.125, 2, 147, None, 2)
testBeam = Beam(20, 0.4, 90666666, 2400)

testResults = beamDeflection(testVehicle, testBeam)

fig, ax = plt.subplots()
line, = ax.plot(testResults[0])

def animationFrames(i):
	line.set_ydata( testResults[i] )
	return(line,)
	
ani = animation.FuncAnimation(fig, animationFrames, np.arange(0, len(testResults)), interval=100)

plt.show()