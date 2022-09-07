import rospy
from gazebo_msgs.msg import ModelStates
from std_msgs.msg import Float64MultiArray, Float32
from geometry_msgs.msg import Pose
from controller import lqr_controller
from quadrotor import linearQuadrotor
import control
import time
import numpy as np

def quad_control(A, B, X, X_d):
	g = 9.8#m/s^2
	m = 3#kg
	U_d = np.array([0, 0, 0, m*g]).T
	Q = np.eye(12)
	R = np.eye(4)
	k_f = 0.3#thrust coeff = b = bl in U input thrust eqns
	k_m = 0.07#moment coeff = d in U input thrust eqns

	#rotor velocities^2:
	sq_w1 = 0
	sq_w2 = 0
	sq_w3 = 0
	sq_w4 = 0
	cmd_vel = np.array([sq_w1, sq_w2, sq_w3, sq_w4]).T
	#K, _, _ = control.lqr(self.A, self.B, Q, R)
	#P=solve_continuous_are(self.A,self.B,Q,R)
	#eig_P=np.linalg.eig(P)
	#K=-np.linalg.inv(R)@np.transpose(self.B)@P
	K, S, E = control.lqr(A, B, Q, R)
	kx = np.dot(K, X_d-X)
	print(kx)
	print(kx.shape)
	print(U_d)
	print(U_d.shape)
	U = U_d - kx


	#now solving the the equation for cmd vels: M.cmd_vel=U => cmd_vel = Minv.U
	M = np.array([[k_f, k_f, k_f, k_f], [k_f, 0, -k_f, 0], [0, k_f, 0, -k_f], [k_m, -k_m, k_m, -k_m]])

	cmd_vel = np.sqrt(np.linalg.inv(M).dot(U))
	print(cmd_vel)

	#-------------------------------------------------------------------------------------------------------------------------------
	#Some Gazebo information for your reference.
	
	#Positive roll is right wing down
	#Positive pitch is front nose down
	#Positive yaw is rotate CCW about z-axis
	
	#Red is x-axis
	#Green is y-axis
	#Blue is z-axis

	#br: Back Right
	#bl: Back Left
	#fl: Front Left
	#fr: Front Right
	#Calculate the ESC pulses (1000us - 2000us PWM signal) for each of the motor.
	
	#br in my code is fr in gazebo's world
	#esc_br = 1500 + output_roll + output_pitch - output_yaw
	#bl in my code is br in gazebo's world
	#esc_bl = 1500 + output_roll - output_pitch + output_yaw
	#fl in my code is bl in gazebo's world
	#esc_fl = 1500 - output_roll - output_pitch - output_yaw
	#fr in my code is fl in gazebo's world
	#esc_fr = 1500 - output_roll + output_pitch + output_yaw
	print(cmd_vel.shape)
	cmd_vel += 1500
	
	#Limit the ESC pulses to upper limit and lower limit, in case the PID algorithm goes crazy and high af.
	if(cmd_vel[0]>2000): cmd_vel[0] = 2000
	if(cmd_vel[1]>2000): cmd_vel[1] = 2000
	if(cmd_vel[2]>2000): cmd_vel[2] = 2000
	if(cmd_vel[3]>2000): cmd_vel[3] = 2000
	#if(esc_br > 2000): esc_br = 2000
	#if(esc_bl > 2000): esc_bl = 2000
	#if(esc_fr > 2000): esc_fr = 2000
	#if(esc_fl > 2000): esc_fl = 2000
	
	if(cmd_vel[0]<1000): cmd_vel[0] = 1000
	if(cmd_vel[1]<1000): cmd_vel[1] = 1000
	if(cmd_vel[2]<1000): cmd_vel[2] = 1000
	if(cmd_vel[3]<1000): cmd_vel[3] = 1000
	#if(esc_br < 1100): esc_br = 1100
	#if(esc_bl < 1100): esc_bl = 1100
	#if(esc_fr < 1100): esc_fr = 1100
	#if(esc_fl < 1100): esc_fl = 1100
	
	cmd_motor_vel=np.zeros([4,1])
	#Map the esc values to motor values
	cmd_motor_vel[0] = ((cmd_vel[0] - 1500)/25) + 50
	cmd_motor_vel[1] = ((cmd_vel[1] - 1500)/25) + 50
	cmd_motor_vel[2] = ((cmd_vel[2] - 1500)/25) + 50
	cmd_motor_vel[3] = ((cmd_vel[3] - 1500)/25) + 50
	#----------------------------------------------------------------------------------------------------------------------------------
	print(cmd_vel)
	return cmd_motor_vel

if __name__ == "__main__":

	i = 0
	iterations=5

	rospy.init_node("quad_control")



	rospy.init_node("Control")

		#initiate publishers that publish errors (roll, pitch,yaw - setpoint) so that it can be plotted via rqt_plot /err_<name>  
		err_rollPub = rospy.Publisher('err_roll', Float32, queue_size=1)
		err_pitchPub = rospy.Publisher('err_pitch', Float32, queue_size=1)
		err_yawPub = rospy.Publisher('err_yaw', Float32, queue_size=1)

		#initialte publisher velPub that will publish the velocities of individual BLDC motors
		velPub = rospy.Publisher('/Kwad/joint_motor_controller/command', Float64MultiArray, queue_size=4)
		#print(velPub)
		#Subscribe to /gazebo/model_states to obtain the pose in quaternion form
		#Upon receiveing the messages, the objects msg, velPub, err_rollPub, err_pitchPub and err_yawPub are sent to "control_kwad" function.
		PoseSub = rospy.Subscriber('/gazebo/model_states',ModelStates,control_kwad,(velPub, err_rollPub, err_pitchPub, err_yawPub))

	while i<iterations:
		#declare a publisher:
		#pub = ....
		quad = linearQuadrotor()
		X_d = np.array([5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0]).T
		cmd = quad_control(quad.A, quad.B, quad.X, X_d)
		print(cmd)
		#read new states:
		#quad.X = sub()

		#calculate output:
		#Y = np.dot(quad.C, X)

		#declare a publisher:
		#pub = ....

		i+=1


		