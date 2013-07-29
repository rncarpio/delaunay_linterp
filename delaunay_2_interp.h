
/*
	This file is part of delaunay_linterp.

    delaunay_linterp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>
#include <string>
#include <vector>
#include <tuple>

#include <boost/heap/skew_heap.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Cartesian_d.h>
#include <Eigen/Dense>

using std::vector;
using std::array;
typedef unsigned int uint;
typedef vector<int> iVec;
typedef vector<double> dVec;

template <class ContainerT, class BitSetT>
ContainerT subset_of_sequence(ContainerT const &set, BitSetT const &subset) {  
  ContainerT result(subset.count());
  auto iter1 = set.begin();
  auto iter2 = result.begin();
  for (int i=0; i<subset.size() && iter1 != set.end(); i++) {
    if (subset[i]) {
	  *iter2 = *iter1;
	  iter2++;
	}
	iter1++;	
  }
  return result;
}

////////////////////////////////////////////////////////////////////////
// CGAL typedefs

//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_euclidean_traits_2<K>  Gt;
  
// this is the class used in CGAL::Triangulation_face_base_with_info_2.
class Face_Info {
public:
  array<K::RT, 2> m_gradient;
};

typedef K::Point_2   Point;
typedef K::RT RT;
typedef K::Line_2 Line_2;
typedef CGAL::Cartesian_d<RT> Kernel_d;
typedef Kernel_d::LA LA;
typedef LA::Matrix Matrix;

typedef CGAL::Triangulation_vertex_base_with_info_2<RT, Gt> Vertex_base;
typedef CGAL::Triangulation_face_base_with_info_2<Face_Info, Gt> Face_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
typedef CGAL::Delaunay_triangulation_2<Gt, TDS> Delaunay_2;
typedef Delaunay_2::Vertex_handle Vertex_handle;
typedef Delaunay_2::Face_handle Face_handle;
typedef Delaunay_2::Edge Edge;
typedef Delaunay_2::Face_circulator Face_circulator;
typedef Delaunay_2::Edge_circulator Edge_circulator;
typedef CGAL::Oriented_side Oriented_side;

typedef array<RT,3> RT_3;
typedef array<RT,2> RT_2;
typedef vector<Face_handle> Face_handle_vec;

double RT_to_double(double const &x) {
  return x;
}

double RT_to_double(CGAL::Lazy_exact_nt<CGAL::Gmpq> const &x) {
  return CGAL::to_double(x);
}

/////////////////////////////////////////////////////////////////////////////////////////
// mathematical functions

template <class IterT, class IterT2>
void get_centroid(int n_dims, IterT points_begin, IterT points_end, IterT2 out_begin) {  
  int n_points = points_end - points_begin;
  std::fill(out_begin, out_begin+n_dims, 0.0);
  for (int i=0; i<n_points; i++) {
    for (int j=0; j<n_dims; j++) {
	  out_begin[j] += points_begin[i][j];
	}
  }
  for (int j=0; j<n_dims; j++) {
    out_begin[j] /= n_points;
  }
  return;
}

// use SVD to find the best-fit k-dimensional subspace, given N points.
// converts to doubles and uses Eigen library.
template <class IterT>
Eigen::MatrixXd best_fit_subspace(int N, int k_dim, IterT points_begin, IterT points_end) {
  int n_points = points_end - points_begin;
  // form the A array. each row is a vector (x_i - x_0)
  Eigen::MatrixXd A(n_points, N);
  for (int i_v=0; i_v<n_points; i_v++) {
    for (int i_dim=0; i_dim<N; i_dim++) {
	  A(i_v, i_dim) = RT_to_double(points_begin[i_v][i_dim]);
	}	
  }
  // compute SVD
  Eigen::JacobiSVD<decltype(A)> svd(A, Eigen::ComputeFullV);	
  // columns of V matrix will form an orthonormal basis for the best-fit subspace. return the first k cols
  return svd.matrixV().block(0, 0, N, k_dim);
}

// compute the squared distance from a point to a vector subspace. the subspace argument is a matrix whose columns are the basis vectors, assumed to be orthonormal.
// returns a tuple containing the squared distance and the closest point.
std::tuple<double, Eigen::VectorXd> sqdist_point_subspace(int N, Eigen::MatrixXd const &subspace, Eigen::VectorXd const &p) {
  // if X is a matrix of orthonormal basis vectors as columns, XX' is the projection matrix onto the subspace
  Eigen::MatrixXd projection = subspace * subspace.transpose();
  // project p onto the subspace
  Eigen::VectorXd p2 = projection * p;
  // the norm of p-p2 is the distance
  Eigen::VectorXd diff = p - p2;
  return std::make_tuple(diff.squaredNorm(), p2);
}

template<class IterT>
Eigen::VectorXd make_eigen_vec(IterT begin, IterT end) {
  int len = end - begin;
  Eigen::VectorXd result(len);
  for (int i=0; i<len; i++) {
    result[i] = RT_to_double(begin[i]);
  }
  return result;
}

// compute the squared distance from a point to the best-fit k-dimensional hyperplane passing through k+1 points
// if k=1, a line. k=2, a plane. etc
// returns the squared distance and the (approximate) closest point on the hyperplane.
template<class NumberT, class IterT, class IterT2, class IterT3>
double sqdist_point_hplane_of_points(int N, IterT p_begin, IterT2 points_begin, IterT2 points_end, IterT3 out_closest_begin) {
  int n_points = points_end - points_begin;
  assert(n_points >= 2);  
  // first, get the centroid of the points and translate everything so that the centroid becomes the origin.
  // this ensures that the best-fit subspace for the k+1 points passes through the centroid.
  vector<NumberT> centroid(N);
  get_centroid(N, points_begin, points_end, centroid.begin());
  Eigen::VectorXd centroid_v = make_eigen_vec(centroid.begin(), centroid.end());
  vector<Eigen::VectorXd> translated_points(n_points);
  for (int i=0; i<n_points; i++) {
    translated_points[i] = make_eigen_vec(points_begin[i].begin(), points_begin[i].begin()+N) - centroid_v;
  }
  Eigen::VectorXd translated_p = make_eigen_vec(p_begin, p_begin+N) - centroid_v;
  
  // find the best-fit subspace for the translated points
  Eigen::MatrixXd subspace(N, n_points-1);
  subspace = best_fit_subspace(N, n_points-1, translated_points.begin(), translated_points.end());
  // find the distance
  double sqdist;
  Eigen::VectorXd projected;
  std::tie(sqdist, projected) = sqdist_point_subspace(N, subspace, translated_p);
  for (int i=0; i<N; i++) {
    out_closest_begin[i] = projected(i);
  }
  return sqdist;
} 

// compute the signed volume of a N-d simplex with N+1 points in N-d space.
template <class IterT>
double compute_volume(IterT points_begin, IterT points_end) {
  int n_points = (points_end - points_begin);
  int n_dims = points_begin->dimension();
  assert(n_points == n_dims + 1);
  Eigen::MatrixXd m(n_dims, n_dims);
  vector<double> temp_m(n_dims*n_dims); //temp
  for (int j=0; j<n_dims; j++) {
	auto direction = points_begin[j+1] - points_begin[0];
	for (int i=0; i<n_dims; i++) {
      m(i, j) = RT_to_double(direction[i]);				// column j is point_{i+1} - point_0
	  temp_m[i*n_dims + j] = m(i,j);
	}
  }
  double det = m.determinant();
  int factorial = 1;
  for (int i=1; i<=n_dims; i++) {
    factorial *= i;
  }
  double result = (det / factorial); // temp
  if (result == 0.0) {
    //assert(false);
  }
  return result;
}

// compute barycentric coordinates for point x in simplex, based on the first N points in the simplex
// returns a vector of RTs
template <class IterT>
vector<RT> compute_barycentric_coords(IterT points_begin, IterT points_end, Point const &x) {  
  int N = x.dimension();
  vector<RT> result(N+1);
  assert(N+1 == points_end - points_begin);	
  auto whole_volume = compute_volume(points_begin, points_end);
  for (int i=0; i<N+1; i++) {
	vector<Point> sub_simplex(points_begin, points_end);
	// replace the ith point with x
	sub_simplex[i] = x;
	auto sub_volume = compute_volume(sub_simplex.begin(), sub_simplex.end());
	double vol1 = sub_volume; // temp
	double vol2 = whole_volume; // temp
	result[i] = sub_volume / whole_volume;
  }
  return result;
}

// uses CGAL's built-in matrix library, so we can use CGAL's number types.
template<int N, class NumberT, class IterT, class IterT2>
array<NumberT,N> simplex_gradient(IterT points_begin, IterT points_end, IterT2 f_begin, IterT2 f_end) {
  assert(points_end - points_begin == N+1);
  assert(f_end - f_begin == N+1);
  Matrix V(N,N);
  Matrix df(N, 1);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
	  V(j,i) = points_begin[i+1][j] - points_begin[0][j];
	}
	df(i,0) = f_begin[i+1] - f_begin[0];
  }
  RT D = 1.0;
  Matrix gradient = LA::transpose(LA::inverse(V, D)) * df;
  assert(gradient.row_dimension() == N);
  assert(gradient.column_dimension() == 1);
  array<NumberT,N> result;
  for (int i=0; i<N; i++) {
    result[i] = gradient(i,0);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////
// hash-related functions

template <class HandleT>
size_t Handle_hash_fn(HandleT const &h) {
  auto hash_fn = CGAL::Handle_hash_function();
  return hash_fn(h);
}

// hash a sequence, element by element
template <class IterT, class HashFnT>
size_t hash_sequence(IterT begin, IterT end, HashFnT const &hash_fn) {
  size_t result = 0;
  int n = end - begin;
  for (int i=0; i<n; i++) {
    result ^= hash_fn(begin[i]);
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////////
// conversion functions

std::string Point_to_string(Point const &p) {
  std::string result("(");
  char pcTemp[64];    
  for (int i=0; i<p.dimension(); i++) {
    sprintf(pcTemp, "%f", RT_to_double(p[i]));
    result += pcTemp;
	if (i < p.dimension()-1) {
	  result += ", ";
	}
  }
  result += ")";
  return result;
}

// a set of Vertex_handles.  Internally, sort by CGAL::Hash (i.e. address) and store in a vector.
// equality comparison is then comparison of the vector.
// hash is combining the hash of each Vertex_handle.
class Ordered_Vertex_List : public vector<Vertex_handle> {
public:
  typedef vector<Vertex_handle> super;
  
  Ordered_Vertex_List(size_t n)
    : super(n)
  {
  }  
  template <class IterT>
  Ordered_Vertex_List(IterT v_begin, IterT v_end)
    : super(v_begin, v_end)
  {
    std::sort(this->begin(), this->end(), [] (Vertex_handle const &a, Vertex_handle const &b) -> bool {
	  return (Handle_hash_fn(a) < Handle_hash_fn(b));
	});
  }
  
  bool operator<(Ordered_Vertex_List const &rhs) const {
    vector<size_t> hash1(this->size()), hash2(rhs.size());
	std::transform(this->begin(), this->end(), hash1.begin(), CGAL::Handle_hash_function());
	std::transform(rhs.begin(), rhs.end(), hash2.begin(), CGAL::Handle_hash_function());
	return std::lexicographical_compare(hash1.begin(), hash1.end(), hash2.begin(), hash2.end());
  }
  
  // return a vector of all ordered subsets of this list (including itself) with a specified minimum size.
  // we want line segments and anything with a higher dimension, so the default min size is 2.
  vector<Ordered_Vertex_List> get_subsets(int min_size=2) {
    vector<Ordered_Vertex_List> result;
	int n = this->size();	
	for (int i=0; i < (1<<n); i++) {
	  boost::dynamic_bitset<> subset(n, i);
	  if (subset.count() >= min_size) {
	    result.push_back(subset_of_sequence(*this, subset));
		//printf("size %d, adding subset of size %d\n", n, result.back().size());
	  }
	}
	return result;
  }
  
  std::string to_string() const {
    std::string result("(");
	for (int i=0; i<this->size(); i++) {
	  result += Point_to_string((*this)[i]->point());
	}
	result += ")";
	return result;
  }
};

namespace std {      
  template <> struct hash<Ordered_Vertex_List> {
    size_t operator()(Ordered_Vertex_List const &s) const {
	  CGAL::Handle_hash_function hash_fn;
	  return hash_sequence(s.begin(), s.end(), hash_fn);
	}
  };
};

bool cmp_lt_Point(Point const &a, Point const &b) {
  for (int i=0; i<a.dimension(); i++) {
    if (a[i] < b[i]) {
	  return true;
	}
  }
  return false;
}

// an ErrorTuple contains:
//  1) x, the point at which f() was evaluated
//  2) f(x)
//  3) the error, i.e. (f(x) - (interpolated value of f(x)))^2
class ErrorTuple {
public:
  double m_err, m_f;
  Point m_x;
  
public:
  ErrorTuple(double err, double f, Point const &x)
    : m_err(err), m_f(f), m_x(x)
  {
    //printf("ErrorTuple() %x\n", this);
  }
  ~ErrorTuple()
  {
    //printf("~ErrorTuple() %x\n", this);
  }
  
  bool operator<(ErrorTuple const &b) const {
    return (m_err < b.m_err);
  }
};

typedef boost::heap::skew_heap<ErrorTuple, boost::heap::mutable_<true> > Error_Priority_Queue;
// for each Face, there may be one or more sub-simplices (including the Face itself) with associated error information.
typedef boost::unordered_map<Face_handle, vector<Error_Priority_Queue::handle_type>, CGAL::Handle_hash_function > Face_ErrorList_Map;
typedef boost::unordered_set<Ordered_Vertex_List, std::hash<Ordered_Vertex_List> > Ordered_Vertex_List_set;

class Delaunay_incremental_interp_2 {
public:
  static const int N = 2;
  std::unique_ptr<Delaunay_2> m_pTriang;  
  Error_Priority_Queue m_err_pqueue;			// priority queue of (error, f, point)
  Face_ErrorList_Map m_face_to_pqueue_handles;	// map of (face -> [pqueue handles])
  std::function<double(int, double*)> m_fn;	// callback function
  
  Delaunay_incremental_interp_2() {
    m_pTriang.reset(new Delaunay_2);	
  }
  Delaunay_incremental_interp_2(std::function<double(int, double*)> fn) {
    m_pTriang.reset(new Delaunay_2);	
	m_fn = fn;
  }  
  virtual ~Delaunay_incremental_interp_2() {	
    m_pTriang.reset();	
  }

  template <class IterT>
  double eval_fn_at_coords(IterT coords_begin, IterT coords_end) const {
	int n_args = coords_end - coords_begin;
	assert(n_args == N);
	array<double,N> args;
	for (int i=0; i<N; i++) {
	  args[i] = RT_to_double(coords_begin[i]);
	}
	double result = m_fn(N, &args[0]);
    return result;
  }
    
  // returns a tuple: (list of modified faces, list of deleted faces)  
  std::tuple<Face_handle_vec, Face_handle_vec> insert_and_get_modifications(Point const &x, double f) {
    std::tuple<Face_handle_vec, Face_handle_vec> result; // first elt is modified list, second is deleted list
    Vertex_handle vh = m_pTriang->insert(x);	      // insert point

    // get all faces adjacent to the new point. these have been modified	
    Face_circulator fc = m_pTriang->incident_faces(vh), done2(fc);
	Face_handle_vec F1, F2;	
    if (fc != NULL) {           
      do {
          std::get<0>(result).push_back(fc);		
	  } while ( ++fc != done2 ); 
    } 
	vh->info() = f;                                 // set f value of vertex
    return result;
  }
 
  template <class IterT>
  void insert(IterT x_begin, IterT x_end, double f) {
    insert_point(Point(x_begin[0], x_begin[1]), f);
  }

  void insert_point(Point const &x, double f) {
    Face_handle_vec modified_faces, deleted_faces;
    std::tie(modified_faces, deleted_faces) = insert_and_get_modifications(x, f);
	
	if (m_fn) {
	  // remove error tuples associated with deleted and modified faces from the queue
	  Face_handle_vec faces_to_remove = modified_faces;
	  faces_to_remove.insert(faces_to_remove.end(), deleted_faces.begin(), deleted_faces.end());
      for (auto i_f=faces_to_remove.begin(); i_f != faces_to_remove.end(); i_f++) {
        auto map_iter = m_face_to_pqueue_handles.find(*i_f);
	    //printf ("Deleting entries for modified simplex %x\n", Handle_hash_fn(*i_f));
	    if (map_iter != m_face_to_pqueue_handles.end()) {	    
		  for (auto i_handle=map_iter->second.begin(); i_handle != map_iter->second.end(); i_handle++) {
		    //printf("erasing %f %f (%f, %f) from queue\n", (**i_handle).m_err, (**i_handle).m_f, (double) (**i_handle).m_x[0], (double) (**i_handle).m_x[1]);
		    m_err_pqueue.erase(*i_handle);
		  }
  	      m_face_to_pqueue_handles.erase(map_iter);	    
        }
	  }
    }	  
	// for modified faces, calculate associated error tuples and add to the queue	      
	for (auto i_s=modified_faces.begin(); i_s != modified_faces.end(); i_s++) {
	  if (!m_pTriang->is_infinite(*i_s)) {		// ignore infinite (i.e. outside the convex hull) faces			
		set_simplex_gradient(*i_s);				// set gradient
		if (m_fn) {
          vector<Vertex_handle> vertices = get_vertices_of_simplex(*i_s);		
	      Ordered_Vertex_List vertices2(vertices.begin(), vertices.end());	
          ErrorTuple err_tuple = compute_simplex_error(vertices2);	  
	      auto pqueue_handle = m_err_pqueue.push(err_tuple);
		  //printf("inserting simplex %x, entry %f %f (%f, %f) into queue\n", Handle_hash_fn(*i_s), err_tuple.m_err, err_tuple.m_f, (double) err_tuple.m_x[0], (double) err_tuple.m_x[1]);
	      m_face_to_pqueue_handles[*i_s].push_back(pqueue_handle);		  
	      //printf ("Inserting simplex %x\n", Handle_hash_fn(*i_s));
	      // handle the facets on the boundary	  
	      add_error_of_boundary_subsimplices(*i_s);	
		}
	  }
	}	  	
	return;
  }
    
  // take the tuple with the largest error and insert it into the triangulation.
  // will not evaluate fn() again, since the value is already stored.
  // returns (x, f(x))
  ErrorTuple insert_largest_error_point() {
    assert(m_err_pqueue.size() > 0);
	assert(m_fn);
	ErrorTuple result = m_err_pqueue.top();
	insert_point(result.m_x, result.m_f);
	return result;
  }

  ErrorTuple get_largest_error_tuple() const {    
    return m_err_pqueue.top();
  }

  template <class IterT>
  double interp(IterT x_begin, IterT x_end) {
    return interp_point(Point(x_begin[0], x_begin[1]));
  }
  
  double interp_point(Point const &x) const {  
    double result;    
    Delaunay_2::Locate_type lt;
    int li;
	Point p = x;
    Face_handle face;
	bool done = false;
	bool found = false;
	while (!done) {
	  lt = Delaunay_2::FACE;
	  face = m_pTriang->locate(p, lt, li);	// locate point p
	  if (face != NULL && !(lt == Delaunay_2::OUTSIDE_CONVEX_HULL || lt == Delaunay_2::OUTSIDE_AFFINE_HULL)) {		// found
	    done = true;
		found = true;
		break;
	  } else if (lt == Delaunay_2::OUTSIDE_CONVEX_HULL && face != NULL) { // outside convex hull
	    // If the point query lies outside the convex hull of the triangulation but in the affine hull, the returned face is an infinite face which is a proof of the point's location :
        // - for a two dimensional triangulation, it is a face (∞, P, Q) such that query lies to the left of the oriented line PQ (the rest of the triangulation lying to the right of this line).
	    // project p onto the line PQ, then try locating again.
		// find the non-infinite vertices.
		int finite_vertices[3];
		int j = 0;
		for (int i=0; i<N+1; i++) {
		  if (! (m_pTriang->is_infinite(face->vertex(i)))) {
		    finite_vertices[j] = i;
			j++;
		  }
		}
		assert(j == 2);
		Line_2 pq(face->vertex(finite_vertices[0])->point(), face->vertex(finite_vertices[1])->point());
		p = pq.projection(p);
	  } else { // not found 
        done = true;       		
	    char pcTemp[1024];
	    sprintf(pcTemp, "interp: no simplex containing the point could be found: face=0x%x, lt=%d", Handle_hash_fn(face), lt);
            throw std::invalid_argument(pcTemp);
        result = std::numeric_limits<double>::quiet_NaN();		
	  }	
    }
	if (found) {
      // get coordinates of vertices	
	  array<Point, N+1> points;	  
	  for (int i=0; i<N+1; i++) {
	    Vertex_handle vertex = face->vertex(i);		  
	    points[i] = vertex->point();
	  }
	  // compute barycentric coords of x
	  auto barycentric_coords = compute_barycentric_coords(points.begin(), points.end(), p);	    
	  RT sum = 0.0;	  
	  for (int i=0; i<barycentric_coords.size(); i++) {
	    Vertex_handle vertex = face->vertex(i);	
	    RT val = vertex->info();
        sum += barycentric_coords[i] * val;		
	  }
      result = RT_to_double(sum);
    }	  
    return result;
  }  

  array<RT,N> set_simplex_gradient(Face_handle f) {
    array<RT,N> result;
	std::fill(result.begin(), result.end(), std::numeric_limits<double>::quiet_NaN());
	//std::fill(result.begin(), result.end(), 0.0);
	if (!m_pTriang->is_infinite(f)) {
	  vector< array<RT, N+1> > vertices_coords = get_coordinates_of_vertices(f);
	  array<RT,N+1> f_vals;
	  for (int i=0; i<N+1; i++) {
	    f_vals[i] = vertices_coords[i][N];
	  }
      result = simplex_gradient<N,RT>(vertices_coords.begin(), vertices_coords.end(), f_vals.begin(), f_vals.end());	  
	}
	f->info().m_gradient = result;
	return result;
  }
 
  array<RT,N> get_simplex_gradient(Face_handle f)  const {
    return f->info().m_gradient;
  }
  
  template <class IterT>
  array<RT,N> get_avg_simplex_gradient(IterT faces_begin, IterT faces_end) const {
    int n_faces = faces_end - faces_begin;
    vector<array<RT,N>> gradients;	
	for (int i=0; i<n_faces; i++) {
	  if (!m_pTriang->is_infinite(faces_begin[i])) {
        gradients.push_back(get_simplex_gradient(faces_begin[i]));
	  }
	}
	array<RT,N> result;
	get_centroid(N, gradients.begin(), gradients.end(), result.begin());
	return result;
  }
  
  array<RT,N> gradient(Point const &x) const {
    array<RT,N> result;	
	std::fill(result.begin(), result.end(), std::numeric_limits<double>::quiet_NaN());      
	//std::fill(result.begin(), result.end(), 0.0);      
    Delaunay_2::Locate_type lt;
    int li;
    Face_handle face;
	lt = Delaunay_2::FACE;
	vector<Face_handle> adjacent_faces;
	face = m_pTriang->locate(x, lt, li);	// locate point p
	if (face != NULL && !(lt == Delaunay_2::OUTSIDE_CONVEX_HULL || lt == Delaunay_2::OUTSIDE_AFFINE_HULL)) {		// found
	  if (lt == Delaunay_2::VERTEX) {
	    assert(!m_pTriang->is_infinite(face));
	    Vertex_handle vh = face->vertex(li);
		// take average of gradients of adjacent faces
        Face_circulator fc = m_pTriang->incident_faces(vh), done2(fc);	
		array<double,N> sum;
		std::fill(sum.begin(), sum.end(), 0.0);
		int n_faces = 0;
        if (fc != NULL) {           
          do {
		    adjacent_faces.push_back(fc);
	      } while ( ++fc != done2 ); 
        }
        result = get_avg_simplex_gradient(adjacent_faces.begin(), adjacent_faces.end());
      } else if (lt == Delaunay_2::EDGE) {
        adjacent_faces.push_back(face);
	    adjacent_faces.push_back(face->neighbor(li));
		result = get_avg_simplex_gradient(adjacent_faces.begin(), adjacent_faces.end());
	  } else if (lt == Delaunay_2::FACE) {
	    result = get_simplex_gradient(face);
	  }
	}
	return result;
  }
  
  // return a vector of vertices of simplex s
  vector<Vertex_handle> get_vertices_of_simplex(Face_handle s, bool filter_empty=true) const {
    vector<Vertex_handle> result;	
	for (int i=0; i<N+1; i++) {
	  Vertex_handle v = s->vertex(i);
	  if (filter_empty == true && !m_pTriang->is_infinite(v)) {
	    result.push_back(v);
	  }
	}
	return result;
  }
  
  // given a d-dimensional simplex, return all simplices of dimension <d that are on the boundary of the triangulation
  vector<Ordered_Vertex_List> get_boundary_subsimplices(Face_handle s) const {
    vector<Ordered_Vertex_List> result;
    int n_dim = N;
	vector<Vertex_handle> vertices;
    for (int i=0; i<n_dim+1; i++) {
	  vertices.push_back(s->vertex(i));
	}
	Face_handle null_handle;
	for (int i=0; i<n_dim+1; i++) {
	  // get the simplex opposite vertex i
	  Face_handle opposite_s = s->neighbor(i);
	  if (m_pTriang->is_infinite(opposite_s)) {
	    // if it doesn't exist, then the facet opposite vertex i must be on the boundary, therefore all vertices of s except i are on the boundary.
		// we will compute all sub-simplices of this boundary simplex by finding all subsequences of this list of vertices on the boundary.
	    boost::dynamic_bitset<> subset(n_dim+1);
	    subset.flip();							// all bits set
	    subset.set(i, 0);							// unset the i-th bit
		vector<Vertex_handle> sub_vertices = subset_of_sequence(vertices, subset);
	    result.push_back(Ordered_Vertex_List(sub_vertices.begin(), sub_vertices.end()));
	    //printf("  adding subsimplex: %s\n", result.back().to_string().c_str());
	  }
	}
	return result;
  }  

  // given a Face s, find its sub-simplices that are on the boundary (if any), calculate their error tuples, and add to the queue
  void add_error_of_boundary_subsimplices(Face_handle const &s) {
	auto boundary_subsimplices_d = get_boundary_subsimplices(s);
	// for each size-d subset, get the d-1, d-2... 2 size subsets
	for (auto i_s2=boundary_subsimplices_d.begin(); i_s2 != boundary_subsimplices_d.end(); i_s2++) {
	  auto boundary_subsimplices_lessthan_d = i_s2->get_subsets();
	  // add each subset to the set of all boundary subsets
	  for (auto i_s3=boundary_subsimplices_lessthan_d.begin(); i_s3 != boundary_subsimplices_lessthan_d.end(); i_s3++) {
        // for each subset, calculate the error and insertion point	
        ErrorTuple err_tuple = compute_simplex_error(*i_s3);	  
	    auto pqueue_handle = m_err_pqueue.push(err_tuple);		
	    m_face_to_pqueue_handles[s].push_back(pqueue_handle);
        //printf("inserting simplex %x, entry %f %f (%f, %f) into queue\n", Handle_hash_fn(s), err_tuple.m_err, err_tuple.m_f, (double) err_tuple.m_x[0], (double) err_tuple.m_x[1]);		
	    //printf ("Inserting subset of size %d: %s\n", i_s3->size(), i_s3->to_string().c_str());	  
      }
	}
    return;
  }

  // given a sequence of vertices, return the coordinates of the points, plus the value of f(x) as an additional coordinate
  template <class IterT>
  vector< array<RT, N+1> > get_coordinates_of_vertices(IterT vertices_begin, IterT vertices_end) const {    
    int n_points = vertices_end - vertices_begin;
    vector< array<RT, N+1> > result(n_points);
	for (int i=0; i<n_points; i++) {
	  for (int j=0; j<N; j++) {
	    result[i][j] = vertices_begin[i]->point()[j];
	  }
	  result[i][N] = vertices_begin[i]->info();
	}
	return result;
  }
  vector< array<RT, N+1> > get_coordinates_of_vertices(Face_handle f) const {
    vector<Vertex_handle> vertices;
    for (int i=0; i<N+1; i++) {
	  if (!m_pTriang->is_infinite(f->vertex(i))) {
	    vertices.push_back(f->vertex(i));
	  }
	}
	return get_coordinates_of_vertices(vertices.begin(), vertices.end());
  }
  
  // compute the error of a simplex of any dimension (e.g. boundary facets)
  ErrorTuple compute_simplex_error(Ordered_Vertex_List const &vertices) const {      
	typedef array<RT, N+1> point_t;    // one additional dimension for the value of f(x)
	vector<point_t> vertex_points = get_coordinates_of_vertices(vertices.begin(), vertices.end());   // get the coordinates of the vertices	    	
	point_t centroid;
	get_centroid(N, vertex_points.begin(), vertex_points.end(), centroid.begin());	// get centroid
	double f_val = eval_fn_at_coords(centroid.begin(), centroid.begin() + N);	// evaluate fn at centroid
	point_t centroid_with_f_val = centroid;
	centroid_with_f_val[N] = RT(f_val);	
    // find distance from p0 to best-fit plane
    double sqdist;
    point_t closest_point;
    sqdist = sqdist_point_hplane_of_points<RT>(N+1, centroid_with_f_val.begin(), vertex_points.begin(), vertex_points.end(), closest_point.begin());			
	// remove the last coordinate from centroid_coords. this is where a new grid point will be added if this simplex were chosen
	Point centroid_point(centroid[0], centroid[1]);	
	ErrorTuple result(sqdist, f_val, centroid_point);
	return result;
  }  

  vector<ErrorTuple> get_error_queue() const {
    vector<ErrorTuple> result(m_err_pqueue.ordered_begin(), m_err_pqueue.ordered_end());
	return result;
  }  
  
  // test for concavity
  bool is_vertex_concave(Vertex_handle vh) const {
    bool all_edges_concave = true;
	Vertex_handle vP, vQ;
    Edge_circulator ec = m_pTriang->incident_edges(vh), done2(ec);	
    if (ec != NULL) {           
      do {
	    if (!m_pTriang->is_infinite(*ec)) {
		  if (!is_edge_concave(*ec, vP, vQ)) {
		    all_edges_concave = false;
			break;
		  }
		}
	  } while ( ++ec != done2 ); 	  
    } 
	return all_edges_concave;
  }
    
  bool is_edge_concave(Edge const &e1, Vertex_handle &vP, Vertex_handle &vQ) const {
    Face_handle face1 = e1.first;
	int i_e1 = e1.second;
	Edge e2 = m_pTriang->mirror_edge(e1);
	Face_handle face2 = e2.first;
	int i_e2 = e2.second;
    vP = face1->vertex(m_pTriang->ccw(i_e1));
	vQ = face1->vertex(m_pTriang->cw(i_e1));
	Point P = vP->point();
	Point Q = vQ->point();
	if (m_pTriang->is_infinite(face1) || m_pTriang->is_infinite(face2)) {
	  return true;
	}
	array<RT,N> gradient1 = this->get_simplex_gradient(face1);
	array<RT,N> gradient2 = this->get_simplex_gradient(face2);	
	RT dgrad_x = gradient2[0] - gradient1[0];
	RT dgrad_y = gradient2[1] - gradient1[1];
    Line_2 PQ(P, Q);	
	Point P_left = Point(P[0] - 1.0, P[1]);
	Point P_down = Point(P[0], P[1] - 1.0);		
    Oriented_side side1 = PQ.oriented_side(face1->vertex(i_e1)->point());
	Oriented_side side_P_left = PQ.oriented_side(P_left);
	Oriented_side side_P_down = PQ.oriented_side(P_down);
	if (side_P_left == CGAL::ON_ORIENTED_BOUNDARY) {  // PQ is horizontal
	  dgrad_x = 0.0;
	} else if (side_P_left != side1) {	     // face1 is to right of PQ
	  dgrad_x = -dgrad_x;
	}
	if (side_P_down == CGAL::ON_ORIENTED_BOUNDARY) {  // PQ is vertical
	  dgrad_y = 0.0;
	} else if (side_P_down != side1) {	     // face1 is above PQ
	  dgrad_y = -dgrad_y;
	}
	bool result = (dgrad_x <= 0.0 && dgrad_y <= 0.0);
	//printf("edge: (%f, %f)->(%f,%f) grad1: %f,%f grad2: %f,%f dgrad_x: %e, dgrad_y: %e, result: %d\n", RT_to_double(P[0]), RT_to_double(P[1]), RT_to_double(Q[0]), RT_to_double(Q[1]),
	  //RT_to_double(gradient1[0]), RT_to_double(gradient1[1]), RT_to_double(gradient2[0]), RT_to_double(gradient2[1]), RT_to_double(dgrad_x), RT_to_double(dgrad_y), result);
	return result;
  }
  //////////////////////////////////////////////////////////////////
  // Plotting-related functions
  vector<array<double,N+1> > get_all_vertices() const {
    vector<array<double,N+1> > result;
    Delaunay_2::Finite_vertices_iterator it;
	array<double,N+1> x;
	for (it = m_pTriang->finite_vertices_begin(); it != m_pTriang->finite_vertices_end(); it++) {
	  for (int i=0; i<N; i++) {
	    x[i] = RT_to_double(it->point()[i]);
	  }
	  x[N] = it->info();
	  result.push_back(x);
	}
	return result;
  }

  vector<std::tuple<Point, int>> get_all_vertices2() const {
    vector<std::tuple<Point, int>> result;
    Delaunay_2::Finite_vertices_iterator it;
	for (it = m_pTriang->finite_vertices_begin(); it != m_pTriang->finite_vertices_end(); it++) {
      bool is_concave = is_vertex_concave(it);
	  result.push_back(std::make_tuple(it->point(), is_concave));
	}
	return result;
  }
  
  vector<std::tuple<Vertex_handle, Vertex_handle> > get_line_segments_of_simplex(Face_handle s) const {
    int n_dims = N;  
    vector<std::tuple<Vertex_handle, Vertex_handle> > result(n_dims + 1);
	for (int i=0; i<n_dims+1; i++) {
	  Vertex_handle v1 = s->vertex(i);
	  Vertex_handle v2 = s->vertex((i+1)%(n_dims+1));
	  result[i] = std::make_tuple(v1, v2);
	}
	return result;
  }
  
  // for plotting. return a vector of all line segments, and the min/max at each coordinate
  std::tuple<vector<std::tuple<array<RT,N+1>, array<RT,N+1>, int> >, array<double,N+1>, array<double,N+1> > get_line_segments() const {
    vector<std::tuple<array<RT,N+1>, array<RT,N+1>, int> > result;
	array<double,N+1> min_coords, max_coords;	
	std::fill(min_coords.begin(), min_coords.end(), DBL_MAX);
	std::fill(max_coords.begin(), max_coords.end(), -DBL_MAX);	
	array<Vertex_handle, 2> vh;
    Delaunay_2::Finite_edges_iterator it;
	for (it = m_pTriang->finite_edges_begin(); it != m_pTriang->finite_edges_end(); it++) {
      bool is_concave = is_edge_concave(*it, vh[0], vh[1]);
	  vector<array<RT,N+1> > vertices_with_f = get_coordinates_of_vertices(vh.begin(), vh.end());
	  result.push_back(std::make_tuple(vertices_with_f[0], vertices_with_f[1], is_concave ? 1 : 0));
	  // find min/max
	  for (int i=0; i<vertices_with_f.size(); i++) {
	    for (int j=0; j<N+1; j++) {
		  min_coords[j] = std::min(min_coords[j], RT_to_double(vertices_with_f[i][j]));
		  max_coords[j] = std::max(max_coords[j], RT_to_double(vertices_with_f[i][j]));
		} 
      }			  
	}  
	return std::make_tuple(result, min_coords, max_coords);
  }  
    
};

		
