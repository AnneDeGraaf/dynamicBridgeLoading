
import scipy.integrate as integrate
import numpy as np 


def main():
	# Perform forceApproximation

	# Then use this force to calculate w of beam
	# plot the results.

	# add some robustness
	return(0)


def beamDeflection(vehicleObject, beamObject, timeStepSize, excitationModes):
	## VERHAALTJE

	# extracting parameters from input vehicleObject 
	mass = verhicleObject.mass 
	velocity = verhicleObject.velocity
	spring = verhicleObject.spring
	damping = verhicleObject.damping
	wheelbase = verhicleObject.wheelbase
	inertia = verhicleObject.inertia 

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

	t1Step = ceal(t1/timeStepSize) + 1		# step number belonging to t1
	t2Step = ceil(t2/timeStepSize) + 1		# step number belonging to t2
	t3Step = ceil(t3/timeStepSize) + 1 		# step number belonging to t3
	nSteps = ceil(tMax/timeStepSize) + 1	# total number of time steps

	# initializing empty arrays for the numerical calculation
	time = np.linspace(0, tMax, nSteps)		# t
	vehicleDisplacement = np.zeros(nSteps)	# u(t)
	vehiclePitch = np.zeros(nSteps)			# theta(t)
	beamDisplacement = np.zeros(nSteps)		# w(t)
	force1 = np.zeros(nSteps)				# Q_1(t)
	force2 = np.zeros(nSteps)				# Q_2(t)

	# initial conditions
	force1[0] = mass*g/2
	force2[0] = mass*g/2


	for step in range(len(nSteps)):

		integral1 = np.trapz( time[step]*time[0:step+1]*force1[0:step+1], time[0:step+1] )		# TIME STEPS NEED ADJUSTING BC HEAVISIDE
		integral2 = np.trapz( time[step]*time[0:step+1]*force2[0:step+1], time[0:step+1] )

		vehicleDisplacement[step] = mass*g/spring + g*time[step]^2/2 - 1/mass*(integral1 + integral2)				
		vehiclePitch[step] = wheelbase/(2*inertia) * (integral2 - integral1)

		# derivVehicleDisplacement = 
		# derivVehiclePitch = 

		for m in range(excitationModes+1)[1:excitationModes+1]:					## w MOET SOWIESO TWEE DIMENSIES HEBBEN 

			eigenfreq = m^2*np.pi*2 / (span^2) * np.sqrt( EI/(density*area) ) 
			duhamel = np.sin( eigenfreq*(time[step]-time[0:step+1]) )			

			phi1 = np.sin(m*np.pi*(velocity*time[step]-wheelbase)/span)			## CHECK DEZE PHIS FF WANT TAU OF t? WANT IS EIGENLIJK FUNCTIE VAN x
			phi2 = np.sin(m*np.pi*(velocity*time[step])/span)

			integral3 = np.trapz( phi1*force1[0:step+1]*duhamel, time[0:step+1] )		## HEAVISIDE
			integral4 = np.trapz( phi2*force2[0:step+1]*duhamel, time[0:step+1] )		## HEAVISIDE

			beamDisplacement[step] = beamDisplacement[step] + 2/(density*area*span*eigenfreq)*(integral3+integral4)

			# derivBeamDisplacement =  np.sin( m*np.pi*x/L ) * 1/(mass*eigenfreq) * int( eigenfreq * np.cos(eigenfreq(time[step]-tau)) * F, tau=0..time[step] )

		# ook hier iets met heaviside doen. w(vt-W) etc. 
		force1[step+1] = springStiffness*( vehicleDisplacement[step] - beamDisplacement[step] + wheelbase/2*vehiclePitch[step] )  # w at x=vt-W
			# + dampingCoefficient*(  derivVehicleDisplacement - derivBeamDisplacement + wheelbase/2*derivVehiclePitch  ))  
		force2[step+1] = springStiffness*( vehicleDisplacement[step] - beamDisplacement[step] + wheelbase/2*vehiclePitch[step] )  # w at x=vt
			# + dampingCoefficient*(  derivVehicleDisplacement - derivBeamDisplacement - wheelbase/2*derivVehiclePitch  ))

		



# some sort of validation needs to take place of the result. 

	# 	if abs((force1LHS - force1RHS)/force1RHS) < relativeError and abs((force2LHS - force2RHS)/force2RHS) < relativeError:
	# 		return (force1RHS, force2RHS)

	# return("Maximum number of iterations reached, no satisfactory result. Try changing the parameters.")




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


