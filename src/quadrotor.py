# Quadrotor Mathematical model for Control
# Following: https://arxiv.org/pdf/1908.07401.pdf 
# For Control: Apply PID to get a target Input Torque/rotation value as per every udated position error in the iterations
# Following: https://www.researchgate.net/profile/Anju-Pillai-2/publication/311950905_Modeling_and_simulation_of_quadcopter_using_PID_controller/links/5bf2e5fe299bf1124fde54a0/Modeling-and-simulation-of-quadcopter-using-PID-controller.pdf?origin=publication_detail

import numpy as np
from numpy import linalg as LA
#import control
import scipy
from scipy.linalg import solve_continuous_are
import math
import warnings

warnings.filterwarnings("ignore")

class linearQuadrotor():

	g = 9.8#m/s^2
	m = 3#kg
	Ix = Iy = Iz = 1	

	A = np.array([[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, -g, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, g, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
				[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
				[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

	B = np.array([[0, 0, 0, 0],
				[0, 0, 0, 0],
				[0, 0, 0, 0],
				[0, 0, 0, 0],
				[0, 0, 0, 0],
				[1/m, 0, 0, 0],
				[0, 0, 0, 0],
				[0, 0, 0, 0],
				[0, 0, 0, 0],
				[0, 1/Ix, 0, 0],
				[0, 0, 1/Iy, 0],
				[0, 0, 0, 1/Iz]])

	C = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
				[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
				[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
				[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]])
		
	#D = np.zeros([6, 4])	
	
	def __init__(self):
		self.x = 1
		self.y = 1
		self.z = 1
		self.x_dot = 1
		self.y_dot = 1
		self.z_dot = 1
		self.phi = 1
		self.theta = 1
		self.psi = 1
		self.phi_dot = 1
		self.theta_dot = 1
		self.psi_dot = 1

		self.X = np.array([self.x, self.y, self.z, self.x_dot, self.y_dot, self.z_dot, 
			 				self.phi, self.theta, self.psi, self.phi_dot, self.theta_dot, self.psi_dot]).T
		self.U = np.array([0, 0, 0, 0]).T	#X = np.arr	X = np.array([1,2,3,0,0,0,])ay([1,2,3,0,0,0,])
		#print(LA.eig(self.A))

	def get_X_dot_vector(self):
		self.X_dot = (np.dot(self.A, self.X) + np.dot(self.B, self.U)).T
		return self.X_dot

	def get_Y_output_vector(self):
		self.Y = (np.dot(self.C, self.X)).T #+ np.dot(self.D, self.U)).T
		return self.Y

	def update_states(self, X_new):
		self.X = X_new

	def lqr(self, A,B,Q,R):
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

	def get_K(self):
		Q = np.eye(12)
		R = np.eye(4)
		#K, _, _ = self.lqr(self.A, self.B, Q, R)
		#P=solve_continuous_are(self.A,self.B,Q,R)
		#eig_P=np.linalg.eig(P)
		#K=-np.linalg.inv(R)@np.transpose(self.B)@P
		#K, S, E = control.lqr(self.A, self.B, Q, R)
		#print(self.X.shape)
		#print(np.dot(K, self.X))
		#return K

	

if __name__ == "__main__":
	quad = linearQuadrotor()
	#print(quad.get_K())
	#a = quad.get_K()
	#print(LA.eig(quad.A))
	#print(a.shape)
	#print(quad.get_K())
	#print(quad.get_K())
	#print(quad.get_Y_output_vector().shape)
	#print(type((quad.get_Y_output_vector())))
	#print(quad.get_X_dot_vector().shape)