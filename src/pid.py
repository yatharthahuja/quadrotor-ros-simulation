#-------------
import time
from quadrotor import linearQuadrotor
import scipy.linalg
#from control.matlab import *
import math
import numpy as np
from sympy.diffgeom import *
from sympy import sqrt,sin,cos
from sympy import *
#--------------------------------------------------------------------------------------------------------------------------------------
def c(ang_radian):
	return math.cos(ang_radian)

def s(ang_radian):
	return math.sin(ang_radian)

def t(ang_radian):
	return math.tan(ang_radian)

def lie_derivtive(B, f, g):  # LfV = delV*f(x)

	M = Manifold("M",3)
	P = Patch("P",M)

	coord          = CoordSystem("coord",P,["x1","x2","x3"])
	x1,x2,x3       = coord.coord_functions()
	e_x1,e_x2,e_x3 = coord.base_vectors()

	f  = x1**2*e_x1 + (sin(x1)+x3**2)*e_x2 + (cos(x3) + x1**2)*e_x3
	g1 = (cos(x1))*e_x1+(x1**2)*e_x2 + 0*e_x3
	g2 = 0*e_x1+ x2*e_x2 + 0*e_x3
	h = x1

	Lg1h = LieDerivative(g1,h)
	Lg2h = LieDerivative(g2,h)
	Lgh = [Lg1h, Lg2h]
	
	return LfB, LgB

def CBF_check(X, U, U_target):
	#QP with CBF
	# avoid locations: 
	U_err = U_target - U

	#constraint(CBF):
	lie_derivtive(B, f, g)[0] + lie_derivtive(B, f, g)[1] +gamma*B>=0

	return 0#U_updated

def dlqr(A,B,Q,R):
    """
    Solve the discrete time lqr controller.
    x[k+1] = A x[k] + B u[k]
    cost = sum x[k].T*Q*x[k] + u[k].T*R*u[k]
    """
    # first, solve the ricatti equation
    P = np.matrix(scipy.linalg.solve_discrete_are(A, B, Q, R))
    # compute the LQR gain
    K = np.matrix(scipy.linalg.inv(B.T*P*B+R)*(B.T*P*A))
    return K

def lqr(A,B,Q,R):
	"""Solve the continuous time lqr controller.
	dx/dt = A x + B u
	cost = integral x.T*Q*x + u.T*R*u
	"""
	#ref Bertsekas, p.151
	#first, try to solve the ricatti equation
	X = np.matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))
	#compute the LQR gain
	K = np.matrix(scipy.linalg.inv(R)*(B.T*X))
	eigVals, eigVecs = scipy.linalg.eig(A-B*K)	 
	return K, X, eigVals


def LQR(x, y, z, roll, pitch, yaw, f):
	#Define the global variables to prevent them from dying and resetting to zero, each time a function call occurs. Some of these variables 		may be redundant.
	global kp_roll, ki_roll, kd_roll, kp_pitch, ki_pitch, kd_pitch, kp_yaw, ki_yaw, kd_yaw, prevErr_roll, prevErr_pitch, prevErr_yaw, pMem_roll, pMem_yaw, pMem_pitch, iMem_roll, iMem_pitch, iMem_yaw, dMem_roll, dMem_pitch, dMem_yaw, flag, setpoint, sampleTime
	#-----------------------
	#Assign your PID values here. From symmetry, control for roll and pitch is the same.
	#kp_roll = 70
	#ki_roll = 0.0002
	#kd_roll = 89
	#kp_pitch = kp_roll
	#ki_pitch = ki_roll
	#kd_pitch = kd_roll
	#kp_yaw = 0.1
	#ki_yaw = 0
	#kd_yaw = 0
	flag = 0
	target_location = [5, 6, 6]
	#Define other variables here, and calculate the errors.
	g = 9.8#m/s^2
	m = 3#kg
	Q = np.eye(12)
	R = np.eye(4)
	k_f = 0.3#thrust coeff = b = bl in U input thrust eqns
	k_m = 0.07#moment coeff = d in U input thrust eqns
	quad = linearQuadrotor()
	A = quad.A
	B = quad.B
	Ix = 0.1
	Iy = 0.1
	Iz = 0.5
	#K, S, E = lqr(A, B, Q, R)
	K = dlqr(A, B, Q, R)
	#K = np.array([[ 0, 0, 1, 0, 0, 2.64575131, 0, 0, 0, 0, 0, 0], 
	#				[ 0, 0, 0, 0, 1.93827322, 0, 13.50882506, 0, 0, 5.29317014, 0, 0],
 	#				[-1, 0, 0, -1.93827322, 0, 0, 0, 13.50882506, 0, 0, 5.29317014, 0],
 	#				[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1.73205081]])

	sampleTime = 1
	setpoint = 0
	err_x = x - target_location[0]
	err_y = y - target_location[1]
	err_z = z - target_location[2]
	err_pitch = float(pitch)*(180 / 3.141592653) - setpoint 
	err_roll = float(roll)*(180 / 3.141592653) - setpoint
	err_yaw = float(yaw)*(180/3.14159263) - setpoint
	currTime = time.time()
	#-----------------------
	#Reset the following variables during the first run only.
	if flag == 0:
		prevTime = 0
		flag += 1
	#------------------------
	dTime = currTime - prevTime
	#-------------------------------------------------------------------------------------------------------------------------------
	#This is the Heart of the PID algorithm. PID behaves more accurately, if it is sampled at regular intervals. You can change the sampleTime to whatever value is suitable for your plant.
	if(dTime >= sampleTime):
		X = np.array([x, y, z, float(pitch)*(180 / 3.141592653), float(roll)*(180 / 3.141592653), float(yaw)*(180 / 3.141592653)])
		#print("state:")
		#print(X)
		#Define dt, dy(t) here for kd calculations.
	
		X_err = np.array([err_x, err_y, err_z, err_pitch, err_roll, err_yaw, 
			 				err_x/dTime, err_y/dTime, err_z/dTime, err_pitch/dTime, err_roll/dTime, err_yaw/dTime]).T
		X_dot = X_err/dTime

		#nonlinear dynamics(x_dot = f(x)+g(x)u):
		#f = X_dot
		
		p = X_err[9]
		q = X_err[10]
		r = X_err[11]
		psi=pitch#float(pitch)*(180 / 3.141592653)#(x) 
		theta=roll#float(roll)*(180 / 3.141592653)#(y) 
		phi=yaw#float(yaw)*(180 / 3.141592653)#(z)
		#f[0] = err_x/dTime
		#f[1] = err_y/dTime
		#f[2] = err_z/dTime
		#f[3] = q*(s(phi)/c(theta))+r(c(phi)/c(theta))
		#f[4] = q*(c(phi))+r(s(phi))
		#f[5] = p+q*(s(phi)*t(theta))+r(c(phi)*t(theta))
		#f[6] = 0
		#f[7] = 0
		#f[8] = g
		#f[9] = ((Iy-Iz)/Ix)*q*r
		#f[10] = ((Iz-Ix)/Iy)*p*r
		#f[11] = ((Ix-Iy)/Iz)*p*q

		#g = C = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		#		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		#		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		#		[0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Ix, 0, 0],
		#		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Iy, 0],
		#		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Iz]])

		#for reference- psi=pitch(x) theta=roll(y) yaw=phi(z)
		
		#g[0][7] = (-1/m)*( s(phi)*s(psi) + c(phi)*c(psi)*s(theta) )
		#g[1][8] = (-1/m)*( s(psi)*s(phi) + c(phi)*c(psi)*s(theta) )
		#g[2][9] = (-1/m)*( c(phi)*c(theta))

		U_target = np.array([m*g, 0, 0, 0]).T
		KX = np.dot(K, X_err)[0]
		print(K.shape)
		print(X_err.shape)
		print(KX.shape)
	
	#Store the current variables into previous variables for the next iteration.
	prevTime = currTime
	prevErr_x = err_x
	prevErr_y = err_y
	prevErr_z = err_z
	prevErr_roll = err_roll
	prevErr_pitch = err_pitch
	prevErr_yaw = err_yaw
	
	#output = Kp*e(t) + Ki*integral(e(t)) + Kd*derivative(e(t))
	U_err = (U_target - KX).T
	#print(KX)
	#print("U_target:")
	#print(U_target)
	#print("U_err:")
	print(U_err)

	#U = CBF_check(U, X) # check and update U values through CBF
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#Some Gazebo information for your reference.
	
	#Positive roll is right wing down
	#Positive pitch is front nose down
	#Positive yaw is rotate CCW about z-axis
	
	#Red is x-axis
	#Green is y-axis
	#Blue is z-axis
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#br: Back Right
	#bl: Back Left
	#fl: Front Left
	#fr: Front Right
	#Calculate the ESC pulses (1000us - 2000us PWM signal) for each of the motor.
	output_roll = U_err[1] 
	output_pitch = U_err[2]
	output_yaw = U_err[3]
	#br in my code is fr in gazebo's world
	esc_br = 100 + output_roll + output_pitch - output_yaw
	#bl in my code is br in gazebo's world
	esc_bl = 100 + output_roll - output_pitch + output_yaw
	#fl in my code is bl in gazebo's world
	esc_fl = 100 - output_roll - output_pitch - output_yaw
	#fr in my code is fl in gazebo's world
	esc_fr = 100 - output_roll + output_pitch + output_yaw

	print(esc_fr)
	print(esc_fl)
	print(esc_br)
	print(esc_bl)
	
	#Limit the ESC pulses to upper limit and lower limit, in case the PID algorithm goes crazy and high af.
	if(esc_br > 700): esc_br = 700
	if(esc_bl > 700): esc_bl = 700
	if(esc_fr > 700): esc_fr = 700
	if(esc_fl > 700): esc_fl = 700
	
	if(esc_br < 10): esc_br = 10
	if(esc_bl < 10): esc_bl = 10
	if(esc_fr < 10): esc_fr = 10
	if(esc_fl < 10): esc_fl = 10
	
	#Map the esc values to motor values
	br_motor_vel = ((esc_br - 100)/25) + 30
	bl_motor_vel = ((esc_bl - 100)/25) + 30
	fr_motor_vel = ((esc_fr - 100)/25) + 30
	fl_motor_vel = ((esc_fl - 100)/25) + 30
	
	#Provide the motor velocities to the object 'f' that will now exit out of this function, and gets published to gazebo, providing velocities to each motor. Note that the sign here is +,-,+,- i.e CW, CCW, CW, CCW in gazebo model. Change view of gazebo model (by scrolling) such that the green line comes to your left, red line goes forward, and blue line goes upward. This is the convention that i refer to as "Gazebo model" incase you get confused.
	f.data = [fr_motor_vel,-fl_motor_vel,bl_motor_vel, -br_motor_vel]
	#[float(U.T[0]), float(U.T[1]), float(U.T[2]), float(U.T[3])]
	#print(f)
	#Return these variables back to the control file.
	return f
#--------------------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------------------
def PID(roll, pitch, yaw, f):
	#Define the global variables to prevent them from dying and resetting to zero, each time a function call occurs. Some of these variables 		may be redundant.
	global kp_roll, ki_roll, kd_roll, kp_pitch, ki_pitch, kd_pitch, kp_yaw, ki_yaw, kd_yaw, prevErr_roll, prevErr_pitch, prevErr_yaw, pMem_roll, pMem_yaw, pMem_pitch, iMem_roll, iMem_pitch, iMem_yaw, dMem_roll, dMem_pitch, dMem_yaw, flag, setpoint, sampleTime
	#-----------------------
	#Assign your PID values here. From symmetry, control for roll and pitch is the same.
	kp_roll = 70
	ki_roll = 0.0002
	kd_roll = 89
	kp_pitch = kp_roll
	ki_pitch = ki_roll
	kd_pitch = kd_roll
	kp_yaw = 0.1
	ki_yaw = 0
	kd_yaw = 0
	flag = 0
	#Define other variables here, and calculate the errors.
	sampleTime = 0
	setpoint = 10
	err_pitch = float(pitch)*(180 / 3.141592653) - setpoint 
	err_roll = float(roll)*(180 / 3.141592653) - setpoint
	err_yaw = float(yaw)*(180/3.14159263) - 30
	currTime = time.time()
	#-----------------------
	#Reset the following variables during the first run only.
	if flag == 0:
		prevTime = 0
		prevErr_roll = 0
		prevErr_pitch = 0
		prevErr_yaw = 0
		pMem_roll = 0
		pMem_pitch = 0
		pMem_yaw = 0
		iMem_roll = 0
		iMem_pitch = 0
		iMem_yaw = 0
		dMem_roll = 0
		dMem_pitch = 0
		dMem_yaw = 0
		flag += 1
	#------------------------
	#Define dt, dy(t) here for kd calculations.
	dTime = currTime - prevTime
	dErr_pitch = err_pitch - prevErr_pitch
	dErr_roll = err_roll - prevErr_roll
	dErr_yaw = err_yaw - prevErr_yaw
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#This is the Heart of the PID algorithm. PID behaves more accurately, if it is sampled at regular intervals. You can change the sampleTime to whatever value is suitable for your plant.
	if(dTime >= sampleTime):
		#Kp*e(t)
		pMem_roll = kp_roll * err_roll
		pMem_pitch = kp_pitch * err_pitch
		pMem_yaw = kp_yaw * err_yaw
		
		#integral(e(t))
		iMem_roll += err_pitch * dTime
		iMem_pitch += err_roll * dTime
		iMem_yaw += err_yaw * dTime
		
		if(iMem_roll > 400): iMem_roll = 400
		if(iMem_roll < -400): iMem_roll = -400
		if(iMem_pitch > 400): iMem_pitch = 400
		if(iMem_pitch < -400): iMem_pitch = -400
		if(iMem_yaw > 400): iMem_yaw = 400
		if(iMem_yaw < -400): iMem_yaw = 400
		
		#derivative(e(t))
		dMem_roll = dErr_roll / dTime
		dMem_pitch = dErr_pitch / dTime
		dMem_yaw = dErr_yaw / dTime
	
	#Store the current variables into previous variables for the next iteration.
	prevTime = currTime
	prevErr_roll = err_roll
	prevErr_pitch = err_pitch
	prevErr_yaw = err_yaw
	
	#output = Kp*e(t) + Ki*integral(e(t)) + Kd*derivative(e(t))
	output_roll = pMem_roll + ki_roll * iMem_roll + kd_roll * dMem_roll
	output_pitch = pMem_pitch + ki_pitch * iMem_pitch + kd_pitch * dMem_pitch
	output_yaw = pMem_yaw + ki_yaw * iMem_yaw + kd_yaw * dMem_yaw 
	#-------------------------------------------------------------------------------------------------------------------------------
		#Ignore this.
	#br_motor_vel = 50.5 + output_pitch + output_roll + output_yaw
	#bl_motor_vel = 50.5 - output_pitch + output_roll - output_yaw
	#fl_motor_vel = 50.5 - output_pitch - output_roll + output_yaw
	#fr_motor_vel = 50.5 + output_pitch - output_roll - output_yaw
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#Some Gazebo information for your reference.
	
	#Positive roll is right wing down
	#Positive pitch is front nose down
	#Positive yaw is rotate CCW about z-axis
	
	#Red is x-axis
	#Green is y-axis
	#Blue is z-axis
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#br: Back Right
	#bl: Back Left
	#fl: Front Left
	#fr: Front Right
	#Calculate the ESC pulses (1000us - 2000us PWM signal) for each of the motor.
	
	#br in my code is fr in gazebo's world
	esc_br = 1500 + output_roll + output_pitch - output_yaw
	#bl in my code is br in gazebo's world
	esc_bl = 1500 + output_roll - output_pitch + output_yaw
	#fl in my code is bl in gazebo's world
	esc_fl = 1500 - output_roll - output_pitch - output_yaw
	#fr in my code is fl in gazebo's world
	esc_fr = 1500 - output_roll + output_pitch + output_yaw
	
	#Limit the ESC pulses to upper limit and lower limit, in case the PID algorithm goes crazy and high af.
	if(esc_br > 2000): esc_br = 2000
	if(esc_bl > 2000): esc_bl = 2000
	if(esc_fr > 2000): esc_fr = 2000
	if(esc_fl > 2000): esc_fl = 2000
	
	if(esc_br < 1100): esc_br = 1100
	if(esc_bl < 1100): esc_bl = 1100
	if(esc_fr < 1100): esc_fr = 1100
	if(esc_fl < 1100): esc_fl = 1100
	
	#Map the esc values to motor values
	br_motor_vel = ((esc_br - 1500)/25) + 50
	bl_motor_vel = ((esc_bl - 1500)/25) + 50
	fr_motor_vel = ((esc_fr - 1500)/25) + 50
	fl_motor_vel = ((esc_fl - 1500)/25) + 50
	#----------------------------------------------------------------------------------------------------------------------------------
	#Ignore this shit here.
	'''
	if(fl_motor_vel > 70): fl_motor_vel = 70
	if(fr_motor_vel > 70): fr_motor_vel = 70
	if(bl_motor_vel > 70): bl_motor_vel = 70
	if(br_motor_vel > 70): br_motor_vel = 70
	
	
	if(err_roll > 0 && err_pitch > 0):
		fl_motor_vel = 51
		fr_motor_vel = 45
		bl_motor_vel = 51
		br_motor_vel = 51
	elif(err_roll > 0 && err_pitch < 0):
		fl_motor_vel = 45
		fr_motor_vel = 51
		bl_motor_vel = 51
		br_motor_vel = 51
	elif(err_roll < 0 && err_pitch > 0):
	'''	
	#---------------------------------------------------------------------------------------------------------------------------------
	#Provide the motor velocities to the object 'f' that will now exit out of this function, and gets published to gazebo, providing velocities to each motor. Note that the sign here is +,-,+,- i.e CW, CCW, CW, CCW in gazebo model. Change view of gazebo model (by scrolling) such that the green line comes to your left, red line goes forward, and blue line goes upward. This is the convention that i refer to as "Gazebo model" incase you get confused.
	f.data = [fr_motor_vel,-fl_motor_vel,bl_motor_vel, -br_motor_vel]
	
	#Return these variables back to the control file.
	return f, err_roll, err_pitch, err_yaw
#--------------------------------------------------------------------------------------------------------------------------------------
