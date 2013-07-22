
# This file is part of delaunay_linterp.

# delaunay_linterp is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import scipy, scipy.stats, scipy.integrate, scipy.optimize, scipy.special, operator, time, itertools, sympy
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

import _delaunay_2_python as _delaunay2
import pyublas
	
# Plot the triangulation in 3d, and allow the user to insert one point at a time	
class IncrInsert_2d:
	def __init__(self, fn, initial_points, n_steps=1):		
		self.plot3d = None		
		self.title = None
		self.n_steps = n_steps
		self.N = 2
		self.initial_points = initial_points
		self.fn = fn
		self.triang = _delaunay2.DelaunayInterp2(fn)				
		# insert the initial points
		for x in initial_points:
			f_val = fn(x)			
			self.triang.insert(x, f_val)			
		print("inserted %d initial points" % len(initial_points))
				
		self.init_draw()
					
	def init_draw(self):
		self.fig = plt.figure()
		self.ax = Axes3D(self.fig)
		self.update_plot()

		self.axbutton = plt.axes([0.81, 0.05, 0.1, 0.075])
		self.button = Button(self.axbutton, 'Insert')
		self.button.on_clicked(self.on_insert_click)		
		plt.draw()
		
	def update_plot(self):
		t0 = time.time()
		if (self.plot3d != None):
			for item in self.plot3d:
				[x.remove() for x in item]
			self.plot3d = None
		[item.remove() for item in self.ax.collections]
		del self.ax.collections[:]		
		
		t1 = time.time()
		(all_segments, min_coords, max_coords) = self.triang.get_line_segments()
		t2 = time.time()
		self.plot3d = []		

		for segment in all_segments:			
			(p0, p1, is_concave) = segment
			color = ('b' if (is_concave) else 'r')
			(xpair, ypair, zpair) = zip(p0, p1)
			self.plot3d.append(self.ax.plot3D(xpair, ypair, zpair, color=color))
		t3 = time.time()
		error_queue = self.triang.get_error_queue()
		if (len(error_queue) > 0):
			(err, f_val, x) = error_queue[0]	
			self.next_point = self.ax.scatter([x[0]], [x[1]], [f_val], color='r')

		self.ax.set_xlim(min_coords[0], max_coords[0])
		self.ax.set_ylim(min_coords[1], max_coords[1])
		self.ax.set_zlim(min_coords[2], max_coords[2])
					
	def on_insert_click(self, event):
		t1 = time.time()
		for i in range(self.n_steps):
			self.triang.insert_largest_error_point()
		t2 = time.time()
		#print("time of insertion: %f, %d points, avg %f" % (t2-t1, self.n_steps, (t2-t1)/self.n_steps))
		
		self.update_plot()
		plt.draw()	

def plot_triang_2d(triang, title=None, show_concave=False):
	fig = plt.figure()
	ax = Axes3D(fig)	
	(all_segments, min_coords, max_coords) = triang.get_line_segments()

	for segment in all_segments:			
		(p0, p1, is_concave) = segment
		if (show_concave):
			color = ('b' if (is_concave) else 'r')
		else:
			color = 'b'
		(xpair, ypair, zpair) = zip(p0, p1)
		ax.plot3D(xpair, ypair, zpair, color=color)
	ax.set_xlim(min_coords[0], max_coords[0])
	ax.set_ylim(min_coords[1], max_coords[1])
	ax.set_zlim(min_coords[2], max_coords[2])
	if (title != None):
		ax.set_title(title)
	return ax

# if x is a list, convert it to a ndarray
# if it is a 1d ndarray, convert it to a column vector
def shape_arg(x):
	if (not isinstance(x, scipy.ndarray)):
		x = scipy.vstack(x)
	if (x.ndim == 1):
		return x.reshape((len(x), 1))
	else:
		return x
		
# Cobb-Douglas utility (note: this is not strictly concave if alphas sum to >= 1.0)
# all methods should be vectorized
# bounds is an arg that can be passed to maximizer functions
class CobbDouglas:
	def __init__(self, alpha_list, bounds=None):
		self.alphas = scipy.array(alpha_list)
		self.n_dims = len(self.alphas)
		if (bounds == None):
			bounds = [(0.001, None)] * self.n_dims
		self.bounds = bounds
		# matrix for inverse gradient
		self.log_alphas = scipy.log(self.alphas)
		self.A = scipy.matrix([self.alphas] * self.n_dims) - scipy.eye(self.n_dims)
		try:
			self.A_is_singular = False
			self.A_inv = scipy.linalg.inv(self.A)
		except(scipy.linalg.LinAlgError):
			self.A_is_singular = True
			
	def __call__(self, x):		return self.eval(x)
	
	def eval(self, x0):
		x = shape_arg(x0)
		assert(len(x) == self.n_dims)
		powers = [xi**alpha for (xi,alpha) in zip(x, self.alphas)]
		return scipy.multiply.reduce(powers)	
	def gradient(self, x0):
		x = shape_arg(x0)
		assert(len(x) == self.n_dims)
		u_array = self(x)
		x2 = u_array * self.alphas.reshape((len(x), 1));		assert(x2.shape == (self.n_dims, len(u_array)))
		assert(x2.shape == x.shape)
		x3 = x2 / x
		assert(x3.shape == x.shape)
		return x3
	def inv_gradient(self, y0):
		y = shape_arg(y0)
		#BREAK()
		assert(len(y) == self.n_dims)
		assert(not self.A_is_singular)
		b = scipy.matrix(scipy.log(y) - self.log_alphas.reshape((self.n_dims, 1)));		assert(b.shape == (self.n_dims, len(y[0])))
		log_x = self.A_inv * b;															assert(log_x.shape == (self.n_dims, len(y[0])))
		x = scipy.exp(scipy.asarray(log_x))	
		assert(scipy.sum(x.imag == 0.0) == x.size)
		return x
		
def test_incr_2d(fn_obj=CobbDouglas([0.3, 0.6])):
	obj = IncrInsert_2d(lambda x: fn_obj(x)[0], [[0.01, 0.01], [2.0, 0.01], [0.01, 2.0], [2.0, 2.0]])
	return obj	
	
	