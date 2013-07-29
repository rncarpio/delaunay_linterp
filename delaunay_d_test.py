
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import _delaunay_d_python as _delaunayd

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
	
def fn(x):
#	return scipy.log(2*x[0] + x[1])
	return x[0]**0.3 * x[1]**0.6

def test1():	
	triang = _delaunayd.Delaunay_incremental_interp_d(2, fn)
	# insert boundary points
	for x in [ [0.01, 0.01], [0.01, 3.0], [3.0, 0.01], [3.0, 3.0] ]:
		triang.insert(x, fn(x))	
	# adaptively place points
	for i in range(50):
		print("inserted %d" % i)
		triang.insert_largest_error_point()	
	# plot triangulation
	plot_triang_2d(triang, title="f(x) = x[0]**0.3 * x[1]**0.6", show_concave=True)

	return triang