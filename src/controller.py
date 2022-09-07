# To add: CBF Control Certificate
# To explore: Model Predictive Control

import numpy as np
import control
import time
from quadrotor import linearQuadrotor as quadrotor

# we need to define higher level control as well
class pd_controller():
	def __init__():
		print("under construction!")

class lqr_controller():
	def __init__(self, A, b, Q, R, U_desired):
		K = control.lqr(A, B, Q, R)
		return K

class mpc_controller():
	def __init__():
		print("under construction!")

class cbf_controller():
	def __init__():
		print("under construction!")