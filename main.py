import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import tkinter as tk
import math as m
import numpy as np
from scipy.stats import norm

# https://www.matburo.ru/Examples/Files/umf_7.pdf
# http://scask.ru/i_book_sub.php?id=10


'''
dU	 	d^2U
-- = 13 ----
dt 		dx^x
function it's true function of U

'''

def function(x, t):
	denominator = 156*t + 1
	exp_numerator = -3*x*x+52*t+2*x
	return (1/m.sqrt(denominator)) * m.exp(exp_numerator/denominator)

class Calculating(tk.Frame):

	def __init__(self, parent, delay, time_start, time_delta, time_max, x_0, x_1, num_parts, show_implicit, show_explicit):
		tk.Frame.__init__(self, parent)
		self.first_start = True
		self.pack(expand=True, fill='both')
		self.delay = delay#time delay, to display graphics between time and timedelta

		self.time_now = time_start
		self.time_delta = time_delta
		self.time_max = time_max
		self.step_for_time = 0

		self.x_0 = x_0#coordinate of the beginning of the rod
		self.x_1 = x_1#coordinate of the ending of the rod
		
		self.num_parts = num_parts#[0,....,99]
		
		self.delta_x = (self.x_1-self.x_0)/(self.num_parts-1)

		self.show_implicit = show_implicit
		self.show_explicit = show_explicit
		
		self.init_2d_list()
		self.explicit_schema()
		self.implicit_schema()
		self.init_plotting()
		self.draw()
	
	def init_2d_list(self):
		
		self.time_list = []
		for t in np.arange(self.time_now, self.time_max+self.time_delta, self.time_delta):
			self.time_list.append(t)
		
		self.x_list = []
		x = self.x_0
		for i in range(self.num_parts):
			self.x_list.append(x)
			x = x + self.delta_x

	def init_plotting(self):
		fig = Figure(figsize=(10, 7), dpi=100)
		self.fig = fig.add_subplot(111)

		self.canvas = FigureCanvasTkAgg(fig, self)
		self.canvas.get_tk_widget().pack(expand=True, fill='both')

	def generate_real_temp(self, time_now):
		
		self.temp_list = []
		for i in range(self.num_parts):
			temp_for_xi_and_time = function(self.x_list[i], time_now)
			self.temp_list.append(temp_for_xi_and_time)

		if self.first_start:
			self.max_t = max(self.temp_list)
			self.min_t = min([0, min(self.temp_list)])
			self.first_start = False

	def init_nodes_with_initial_values(self):
		'''returns a 2d array with a known temperature value by initial conditions:
			* initial condition at time self.time_now
			* boundary condition, at the beginning and at the end of the rod at any time

			for example:
								t(j)
								^
				self.time_max - | 0.3 nan nan  nan nan nan  nan 0.3
								| 0.4 nan nan  nan nan nan  nan 0.4
								| 0.5 nan nan  nan nan nan  nan 0.5
								| 0.6 nan nan  nan nan nan  nan 0.6
				self.time_now -	| 0.7 0.7 0.75 0.8 0.9 0.75 0.7 0.7
								 ----------------------------------->x(i)
								   ^							 ^
								 self.x_0 					self.x_1
		'''
		temp_list = []
		for j in range(len(self.time_list)):
			temp_list_for_j = []
			for i in range(self.num_parts):
				temp = np.nan
				if j == 0:
					temp = function(self.x_list[i], self.time_list[j])
				elif i == 0 or i == self.num_parts - 1:
					temp = function(self.x_list[i], self.time_list[j])
				temp_list_for_j.append(temp)
			temp_list.append(temp_list_for_j)
		return temp_list
	
	def explicit_schema(self):
		self.explicit_temp_list = self.init_nodes_with_initial_values()
		const = 13*self.time_delta/(self.delta_x**2)

		#start from 1 index not 0, because we know in 0 index temperature value by initial conditions
		for j in range(1, len(self.time_list)):
			for i in range(1, self.num_parts - 1):
				if np.isnan(self.explicit_temp_list[j][i]):
					ui_left = self.explicit_temp_list[j-1][i-1]
					ui = self.explicit_temp_list[j-1][i]
					ui_right = self.explicit_temp_list[j-1][i+1]
					self.explicit_temp_list[j][i] = ui + const*(ui_left - 2*ui + ui_right)

		#was need for testing, but in explicit schema could be nan and +-inf, if system crashes
		# for i1 in self.explicit_temp_list:
		# 	for i2 in i1:
		# 		if np.isnan(i2):
		# 			exit('some nan in explicit_temp_list')

	def implicit_schema(self):
		self.implicit_temp_list = self.init_nodes_with_initial_values()
		K = 13*self.time_delta/self.delta_x
		
		Ai = K
		Ci = 2*K + 1
		Bi = K

		#start from 1 index not 0, because we know in 0 index temperature value by initial conditions
		for j in range(len(self.time_list)-1):
			alpha = [0]
			beta = [self.implicit_temp_list[j+1][0]]

			#forward
			for i in range(1, self.num_parts):
				new_alpha = Bi/(Ci-alpha[i-1]*Ai)
				alpha.append(new_alpha)
				new_beta = (Ai * beta[i-1] + self.implicit_temp_list[j][i-1]) / (Ci-alpha[i-1]*Ai)
				beta.append(new_beta)
			#back
			for i in range(self.num_parts-2, 0, -1):
				self.implicit_temp_list[j+1][i] = alpha[i+1]*self.implicit_temp_list[j+1][i+1] + beta[i+1]
		
		for i1 in self.implicit_temp_list:
			for i2 in i1:
				if np.isnan(i2):
					exit('some nan in implicit_temp_list')
	def draw(self):

		self.fig.clear()
		
		if not self.first_start:
			self.fig.axis((self.x_0-0.005, self.x_1+0.005, self.min_t-0.005, self.max_t+0.005))

		self.generate_real_temp(self.time_now)

		#plot real function of temperature
		self.fig.plot(self.x_list, self.temp_list, "r-", label=f'Real')
		
		
		#plot explicit
		if self.show_explicit:
			self.fig.plot(self.x_list, self.explicit_temp_list[self.step_for_time], "b-", label=f'explicit')

		#plot implicit
		if self.show_implicit:
			self.fig.plot(self.x_list, self.implicit_temp_list[self.step_for_time], "g-", label=f'implicit')

		self.fig.set_xlabel(f'time: {np.round(self.time_now, 2)}')
		self.fig.legend()
		self.canvas.draw()

		
		if self.time_now > self.time_max:
			if input('press q for exit:\n') == 'q':
				exit()

		self.time_now += self.time_delta
		self.step_for_time += 1
		#self.time_list[self.step_for_time] == self.time_now  ---->it's equal
		self.after(self.delay, self.draw)



if __name__ == '__main__':


	delay = 0#time delay in milisec, to display graphics between time and timedelta

	time_start = 0
	time_delta = 0.01#step for time
	time_max = 20#max time

	x_0 = -50#coordinate of the beginning of the rod
	x_1 = 50#coordinate of the ending of the rod
	
	num_parts = 300#[0,....,299] number of parts for breaking the rod
	
	show_implicit = True
	show_explicit = False

	root = tk.Tk()
	app = Calculating(root, delay, time_start, time_delta, time_max, x_0, x_1, num_parts, show_implicit, show_explicit)
	root.mainloop()







