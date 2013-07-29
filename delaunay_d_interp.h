
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>
#include <string>
#include <vector>
#include <tuple>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/heap/skew_heap.hpp>
#include <boost/dynamic_bitset.hpp>

//#include <unordered_map>
//#include <unordered_set>
#include <map>
#include <set>

#include <functional>

#include <CGAL/Cartesian_d.h>
#include <CGAL/Quotient.h>
#include <CGAL/Delaunay_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Object.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Lazy_exact_nt.h>

#include <Eigen/Dense>

using std::vector;
using std::array;
typedef unsigned int uint;
typedef vector<int> iVec;
typedef vector<double> dVec;

//typedef CGAL::MP_Float NT;
//typedef CGAL::Quotient<NT> RT;	
typedef CGAL::Gmpz NT;
typedef CGAL::Gmpq RT;
//typedef double RT;
//typedef double NT;
//typedef CGAL::Lazy_exact_nt<CGAL::Gmpq> NT;
//typedef CGAL::Lazy_exact_nt<CGAL::Gmpq> RT;

typedef CGAL::Cartesian_d<RT> Kernel;
typedef Kernel::LA LA;                 
typedef LA::Matrix Matrix;
typedef LA::Vector LA_Vector;
typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
typedef CGAL::Delaunay_d<Kernel> Delaunay_d;
typedef CGAL::Hyperplane_d<Kernel> Hyperplane;
typedef CGAL::Direction_d<Kernel> Direction;
typedef CGAL::Ray_d<Kernel> Ray;
typedef Delaunay_d::Point_d Point;
typedef Delaunay_d::Vector_d Vector;
typedef Convex_hull_d::Facet_handle Facet_handle;
typedef Delaunay_d::Simplex_handle Simplex_handle;
typedef Delaunay_d::Vertex_handle Vertex_handle;
typedef CGAL::Object Object;

typedef vector<Simplex_handle> Simplex_handle_vec;
typedef vector<Vertex_handle> Vertex_handle_vec;

/////////////////////////////////////////////////////////////////////
// forward declarations

double RT_to_double(RT const &x);

double RT_to_double(double const &x) {
  return x;
}
double RT_to_double(CGAL::Gmpq const &x) {
  return to_double(x);
}
double RT_to_double(CGAL::Quotient<NT> const &x) {
  return to_double(x);
}

template <class HandleT>
size_t Handle_hash_fn(HandleT const &h) {
  auto hash_fn = CGAL::Handle_hash_function();
  return hash_fn(h);
}

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

/////////////////////////////////////////////////////////
// hash functions

// hash an arbitrary buffer of data.
// pad with zeros into a buffer that is a multiple of sizeof(size_t), xor, then hash
size_t hash_buffer(unsigned char const *pSrc, int nBytes) {
  int reps = nBytes / sizeof(size_t);
  size_t combined = 0;
  unsigned char *pCombined = (unsigned char *) &combined;
  int i;
  for (i=0; i<reps*sizeof(size_t); i += sizeof(size_t)) {
    combined ^= (* (size_t*) (pSrc + i));
  }
  for (int j=0; j<nBytes-i; j++) {
    pCombined[j] ^= pSrc[i+j];
  }  
  return std::hash<size_t>()(combined);
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

// put hash functions for custom types here
namespace std {
  // MP_Floats have a vector member v, and a scalar member exp
  template <> struct hash<CGAL::MP_Float> {  
    size_t operator()(CGAL::MP_Float const &x) const {
	  return hash_buffer( (unsigned char const*) &x.v[0], x.v.size()*sizeof(decltype(x.v[0])) ) ^ std::hash<CGAL::MP_Float::exponent_type>()(x.exp);
    }
  };

  template <class T> struct hash<CGAL::Lazy_exact_nt<T> > {  
    size_t operator()(CGAL::Lazy_exact_nt<T> const &x) const {
	  return std::hash<double>()(x.exact().to_double());
    }
  };

  template <> struct hash<CGAL::Gmpq> {  
    size_t operator()(CGAL::Gmpq const &x) const {
	  return std::hash<double>()(x.to_double());
    }
  };
  
  // Quotients have two members, num and den
  // template <> struct hash<CGAL::Quotient<CGAL::MP_Float> > {  
    // size_t operator()(CGAL::Quotient<CGAL::MP_Float> const &x) const {
	  // typedef CGAL::Quotient<CGAL::MP_Float>::NT NT;
	  // return std::hash<NT>()(x.num) ^ std::hash<NT>()(x.den);
    // }
  // };

  template <> struct hash<CGAL::Quotient<CGAL::MP_Float> > {  
    size_t operator()(CGAL::Quotient<CGAL::MP_Float> const &x) const {	  
	  return std::hash<CGAL::MP_Float>()(x.num) ^ std::hash<CGAL::MP_Float>()(x.den);
    }
  };

  // for a Point, get each cartesian coordinate
  template <> struct hash<Point> {  
    size_t operator()(Point const &p) const {
      return hash_sequence(p.cartesian_begin(), p.cartesian_end(), std::hash<RT>());	  
    }
  };
  
  // this version just combines hashes (which is the pointer to the underlying object) of each Vertex_handle.
  // should work as long as vertices are not deleted
  template <> struct hash<Simplex_handle> {
    Delaunay_d const *m_pDelaunay;
    hash(Delaunay_d const *pDelaunay)
      : m_pDelaunay(pDelaunay)
    {}
  
    size_t operator()(Simplex_handle const &s) const {
      int n = m_pDelaunay->dimension();
	  size_t result = 0;
	  for (int i=0; i<n+1; i++) {
	    Vertex_handle v = m_pDelaunay->vertex_of_simplex(s, i);
	    result ^= CGAL::Handle_hash_function()(v);
	  }
	  return result;
	}
  };
  
  template <> struct hash<Ordered_Vertex_List> {
    size_t operator()(Ordered_Vertex_List const &s) const {
	  CGAL::Handle_hash_function hash_fn;
	  return hash_sequence(s.begin(), s.end(), hash_fn);
	}
  };

  template <> struct hash<CGAL::Gmpzf> {  
    size_t operator()(CGAL::Gmpzf const &x) const {	  
	  return std::hash<double>()(x.to_double());
    }
  };
  
  // for a Simplex of a Delaunay triangulation, we're only interested in the first N coordinates of each vertex (the last one is from the lifting map).
  // lexicograpically sort the points, then hash
  // template <> struct hash2<Simplex_handle> {
    // Delaunay_d const *m_pDelaunay;
    // hash2(Delaunay_d const &triang)
      // : m_pDelaunay(&triang)
    // {}
  
    // size_t operator()(Simplex_handle const &s) const {
      // int n = m_pDelaunay->dimension();
	  // vector<Point> points;
	  // for (int i=0; i<n; i++) {
	    // Point const &p(m_pDelaunay->point_of_simplex(s, i));
	    // points.push_back( Point(n, p.cartesian_begin(), p.cartesian_end()-1) );
	  // }
	  // std::sort(points.begin(), points.end());		// Point has a built-in lexicographic comparison
	  // return hash_sequence(points.begin(), points.end(), std::hash<Point>);
	// }
  // };
};
  
// PointCd is a handle for Tuple_d<_FT,_LA>	  
// Tuple_d contains a Vector.. what's that?
// Vector is LinearAlgebra::Vector_ in Linear_algebraCd
//   which has members v_, d_

typedef boost::unordered_map<Point, double, std::hash<Point> > Point_d_map;
typedef boost::unordered_map<Simplex_handle, double, std::hash<Simplex_handle> > Simplex_d_map;		// this will work as long as the underlying vertices don't go away
typedef boost::unordered_map<Vertex_handle, double, CGAL::Handle_hash_function > Vertex_d_map;		

// typedef std::map<Point, double> Point_d_map;
// typedef std::map<Simplex_handle, double> Simplex_d_map;
// typedef std::map<Vertex_handle, double> Vertex_d_map;		



/*
// given at least N affinely independent points in N-dimensional space, 
// calculate the coefficents of the best-fit plane equation.  If N affinely independent points are given, the plane
// will pass through all points.  If >N points are given, use total least squares to find the best fit.
// returns a vector of N+1 coefficients, e.g. ax + by + cz + d = 0 for N=3. the first N coefs are the normal vector to the plane
template <class IterT>
Hyperplane best_fit_plane(IterT points_begin, IterT points_end) {
  int n_points = (points_end - points_begin);
  int n_dims = points_begin->dimension();
  dVec result(n_dims+1);
  assert(n_points >= n_dims);
  // A is (n_points x n_dims+1), last col is 1's
  // doubles are good enough for this purpose.
  Eigen::MatrixXd A(n_points, n_dims+1);
  for (int i_v=0; i_v<n_points; i_v++) {
    for (int i_dim=0; i_dim<n_dims; i_dim++) {
	  A(i_v, i_dim) = RT_to_double(points_begin[i_v].cartesian(i_dim));
	}
	A(i_v, n_dims) = 1.0;		// fill in last col with 1
  }
  //Eigen::JacobiSVD<decltype(A)> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);	
  Eigen::JacobiSVD<decltype(A)> svd(A, Eigen::ComputeFullV);	
  // the last column of V is the null vector
  for (int i=0; i<n_dims+1; i++) {
	result[i] = svd.matrixV()(i, svd.matrixV().cols()-1);
  }
  // if n_points == n_dims, use sign of determinant to find orientation
  if (n_points == n_dims) {
	Matrix m(n_dims, n_dims);
	int i;
	for (i=0; i<n_points-1; i++) {
	  auto v = points_begin[i+1] - points_begin[0];
	  for (int j=0; j<n_dims; j++) {
		m(i, j) = v[j];				// each row is point_{i+1} - point_0
      }
	}
	for (int j=0; j<n_dims; j++) {
      m(i, j) = result[j];			// last row is the normal vector
    }	    
    if (LA::sign_of_determinant(m) < 0) {
	  for (int i=0; i<result.size(); i++) {
	    result[i] = -result[i];
	  }
    }
  }
  vector<RT> coefs(result.begin(), result.end());
  Hyperplane hplane(n_dims, coefs.begin(), coefs.end());
  return hplane;
}
	
// squared distance and closest point from a hyperplane to a point
std::tuple<RT, Point> sqdist_point_hyperplane(Hyperplane const &hplane, Point const &p) {
  if (hplane.has_on_boundary(p)) {
    return std::make_tuple(RT(0.0), p);
  } else {
    Direction d = hplane.orthogonal_direction();
	if (hplane.has_on_positive_side(p)) {
	  d = -d;
	}
	Ray r(p, d);
	Object result = intersection(r, hplane);
    const Point *ipnt = CGAL::object_cast<Point>(&result);
	assert(ipnt != NULL);
	RT sqdist = CGAL::squared_distance(p, *ipnt);
	return std::make_tuple(sqdist, *ipnt);
  }
}

template <class IterT>
Point get_centroid(int n_dims, IterT points_begin, IterT points_end) {  
  int n_points = points_end - points_begin;
  vector<RT> centroid_coords(n_dims, 0.0);	
  for (int i=0; i<n_points; i++) {
    for (int j=0; j<n_dims; j++) {
	  centroid_coords[j] += points_begin[i][j];
	}
  }
  for (int j=0; j<n_dims; j++) {
    centroid_coords[j] /= n_points;
  }
  return Point(n_dims, centroid_coords.begin(), centroid_coords.end());
}

// use SVD to find the best-fit k-dimensional subspace, given N points.
// if N==k, returns a subspace that passes through the points.
template <class IterT>
std::tuple<Point, Eigen::MatrixXd> best_fit_subspace(int n_dims, int k_dim, IterT points_begin, IterT points_end) {  
  int n_points = points_end - points_begin;
  // get the centroid of the points and use it as the origin
  Point centroid = get_centroid(n_dims, points_begin, points_end);
  // subtract the centroid from each point to get a sequence of Vectors
  vector<Vector> demeaned_points(n_points);
  for (int i=0; i<n_points; i++) {
    demeaned_points[i] = points_begin[i] - centroid;
  }  
  // form the A array. each row is a Vector, then add a column of 1's for the intercept
  Eigen::MatrixXd A(n_points, n_dims);
  for (int i_v=0; i_v<n_points; i_v++) {
    for (int i_dim=0; i_dim<n_dims; i_dim++) {
	  A(i_v, i_dim) = RT_to_double(demeaned_points[i_v].cartesian(i_dim));
	}
	//A(i_v, n_dims) = 1.0;		// fill in last col with 1
  }
  // compute SVD
  Eigen::JacobiSVD<decltype(A)> svd(A, Eigen::ComputeFullV);	
  // columns of V matrix will form an orthonormal basis for the best-fit subspace. return the first k cols
  return std::make_tuple(centroid, svd.matrixV().block(0, 0, n_dims, k_dim));
}

// compute the squared distance from a point to a vector subspace. the subspace argument is a matrix whose columns are the basis vectors, assumed to be orthonormal.
std::tuple<double, Eigen::VectorXd> sqdist_point_subspace(Eigen::MatrixXd const &subspace, Vector const &point) {
  // convert p into a column vector of doubles
  Eigen::VectorXd p(point.dimension());
  for (int i=0; i<point.dimension(); i++) {
    p[i] = RT_to_double(point[i]);
  }    
  // if X is a matrix of orthonormal basis vectors as columns, XX' is the projection matrix onto the subspace
  Eigen::MatrixXd projection = subspace * subspace.transpose();
  // project p onto the subspace
  Eigen::VectorXd p2 = projection * p;
  // the norm of p-p2 is the distance
  Eigen::VectorXd diff = p - p2;
  return std::make_tuple(diff.squaredNorm(), p2);
}

// compute the squared distance from a point to the best-fit k-dimensional subspace of k+1 points
// if k=1, a line. k=2, a plane. etc
// returns the squared distance and the (approximate) closest point on the subspace
template<class IterT>
std::tuple<RT, Point> sqdist_point_subspace_of_points(Point const &p, IterT points_begin, IterT points_end) {
  int n_points = points_end - points_begin;
  assert(n_points >= 2);
  int n_dims = p.dimension();
  assert(p.dimension() == points_begin[0].dimension());
  Point centroid;
  Eigen::MatrixXd subspace(n_dims, n_points-1);
  // get the centroid of the points and the best-fit subspace
  std::tie(centroid, subspace) = best_fit_subspace(n_dims, n_points-1, points_begin, points_end);
  // use the centroid as our origin
  Vector v = p - centroid;
  // find the distance
  double sqdist;
  Eigen::VectorXd projected;
  std::tie(sqdist, projected) = sqdist_point_subspace(subspace, v);
  vector<double> coords(n_dims);
  for (int i=0; i<n_dims; i++) {
    coords[i] = projected(i);
  }
  std::tuple<RT, Point> result(sqdist, Point(n_dims, coords.begin(), coords.end()));
  return result;
}  

// compute the signed volume of a N-d simplex with N+1 points in N-d space.
template <class IterT>
RT compute_volume(IterT points_begin, IterT points_end) {
  int n_points = (points_end - points_begin);
  int n_dims = points_begin->dimension();
  assert(n_points == n_dims + 1);
  Matrix m(n_dims, n_dims);
  vector<double> temp_m(n_dims*n_dims); //temp
  for (int j=0; j<n_dims; j++) {
	auto direction = points_begin[j+1] - points_begin[0];
	for (int i=0; i<n_dims; i++) {
      m(i, j) = direction[i];				// column j is point_{i+1} - point_0
	  //temp_m[i*n_dims + j] = m(i,j).to_double();
	}
  }
  auto det = LA::determinant(m);
  int factorial = 1;
  for (int i=1; i<=n_dims; i++) {
    factorial *= i;
  }
  RT result = (det / factorial); // temp
  if (result == 0) {
//    assert(false);
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
//	double vol1 = sub_volume.to_double(); // temp
//	double vol2 = whole_volume.to_double(); // temp
	result[i] = sub_volume / whole_volume;
  }
  return result;
}
*/

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
double compute_volume(int n_dims, IterT points_begin, IterT points_end) {
  int n_points = (points_end - points_begin);  
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
  auto whole_volume = compute_volume(N, points_begin, points_end);
  if (whole_volume == 0.0) {
    for (int i=0; i<N+1; i++) {
	  std::cout << points_begin[i] << std::endl;
	}
	assert(false);
  }
  for (int i=0; i<N+1; i++) {
	vector<Point> sub_simplex(points_begin, points_end);
	// replace the ith point with x
	sub_simplex[i] = x;
	auto sub_volume = compute_volume(N, sub_simplex.begin(), sub_simplex.end());
	double vol1 = sub_volume; // temp
	double vol2 = whole_volume; // temp
	result[i] = sub_volume / whole_volume;
  }
  return result;
}

// uses CGAL's built-in matrix library, so we can use CGAL's number types.
template<class NumberT, class IterT, class IterT2>
vector<NumberT> simplex_gradient(int N, IterT points_begin, IterT points_end, IterT2 f_begin, IterT2 f_end) {
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
/*
  printf("points:\n");
  for (int i=0; i<N+1; i++) {
    for (int j=0; j<N; j++) {
	  printf("%f ", RT_to_double(points_begin[i][j]));
	}
	printf("\n");
  }
  printf("V:\n");
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
	  printf("%f ", RT_to_double(V(i,j)));
	}
	printf("\n");
  }
*/  

  RT D = 1.0;
  Matrix gradient = LA::transpose(LA::inverse(V, D)) * df;
  assert(gradient.row_dimension() == N);
  assert(gradient.column_dimension() == 1);
  vector<NumberT> result(N);
  for (int i=0; i<N; i++) {
    result[i] = gradient(i,0);
  }
  return result;
}

bool cmp_lt_Point(Point const &a, Point const &b) {
  for (int i=0; i<a.dimension(); i++) {
    if (a[i] < b[i]) {
	  return true;
	}
  }
  return false;
}

typedef boost::unordered_set<Simplex_handle, CGAL::Handle_hash_function> Simplex_handle_set;	// this one uses the pointer to the Simplex
typedef boost::unordered_set<Simplex_handle, std::hash<Simplex_handle> > Simplex_value_set;	// this one is by value (uses CGAL::Handle_hash_function on each Vertex_handle)
typedef boost::unordered_map<Simplex_handle, double, CGAL::Handle_hash_function> Simplex_handle_d_map;
typedef boost::unordered_map<Simplex_handle, double, std::hash<Simplex_handle> > Simplex_value_d_map;
typedef boost::unordered_set<Vertex_handle, CGAL::Handle_hash_function> Vertex_set;
typedef boost::unordered_map<Vertex_handle, Simplex_value_set, CGAL::Handle_hash_function> Vertex_Simplex_set_map;
typedef boost::unordered_set<Ordered_Vertex_List, std::hash<Ordered_Vertex_List> > Ordered_Vertex_List_set;

// typedef boost::unordered_set<Simplex_handle, CGAL::Handle_hash_function> Simplex_handle_set;	// this one uses the pointer to the Simplex
// typedef std::set<Simplex_handle> Simplex_value_set;	// this one is by value (uses CGAL::Handle_hash_function on each Vertex_handle)
// typedef boost::unordered_map<Simplex_handle, double, CGAL::Handle_hash_function> Simplex_handle_d_map;
// typedef std::map<Simplex_handle, double> Simplex_value_d_map;
// typedef std::set<Vertex_handle> Vertex_set;
// typedef std::map<Vertex_handle, Simplex_value_set> Vertex_Simplex_set_map;
// typedef std::set<Ordered_Vertex_List> Ordered_Vertex_List_set;

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
  std::string to_str() const {
    std::ostringstream result;
	result << "err: " << m_err << ", f: " << m_f << ", x: " << m_x;
	return result.str();
  }
};

typedef std::tuple<Vertex_handle, vector<Simplex_handle> > InsertResult;
//typedef boost::heap::skew_heap<ErrorTuple, boost::heap::compare<cmp_lt_ErrorTuple>, boost::heap::mutable_<true> > Error_Priority_Queue;
typedef boost::heap::skew_heap<ErrorTuple, boost::heap::mutable_<true> > Error_Priority_Queue;
typedef boost::unordered_map<Simplex_handle, vector<Error_Priority_Queue::handle_type>, CGAL::Handle_hash_function> Simplex_Error_Map;

class Delaunay_incremental_interp_d {  
public:  
  int N, m_N;
  boost::shared_ptr<Delaunay_d> m_pDelaunay;
  Vertex_d_map m_pointMap;
  std::function<double(int, double*)> m_fn;	// callback function
  
  // if we want to store additional info with a simplex, we have to do it ourselves.
  // we can assume that simplices are never destroyed, with the randomized incremetal algorithm.
  struct SimplexInfo {
    double m_volume;
    vector<RT> m_gradient;
	vector<Error_Priority_Queue::handle_type> m_pqueue_handles;
  };

  typedef boost::unordered_map<Simplex_handle, SimplexInfo, CGAL::Handle_hash_function> Simplex_Info_Map;  
  Simplex_Info_Map m_simplex_info_map;
  Error_Priority_Queue m_err_pqueue;

  Delaunay_incremental_interp_d(int n_dims) {
    N = m_N = n_dims;  
    m_pDelaunay.reset(new Delaunay_d(m_N));
  }
  
  Delaunay_incremental_interp_d(int n_dims, std::function<double(int, double*)> fn) {	
    N = m_N = n_dims;
    m_pDelaunay.reset(new Delaunay_d(m_N));
	m_fn = fn;
  }
  
  int get_dimension() const {
    return m_N;
  }

  template <class IterT>
  double eval_fn_at_coords(IterT coords_begin, IterT coords_end) const {
	int n_args = coords_end - coords_begin;
	assert(n_args == N);
	vector<double> args(n_args);
	for (int i=0; i<N; i++) {
	  args[i] = RT_to_double(coords_begin[i]);
	}
	double result = m_fn(N, &args[0]);
    return result;
  }
  
  RT get_val_of_vertex(Vertex_handle v) const {
    RT result;
    auto map_iter = m_pointMap.find(v);
	if (map_iter != m_pointMap.end()) {		// found
	  result = map_iter->second;	  
	} else {
	  assert(false);
	}
	return result;
  }

  vector<Ordered_Vertex_List> get_boundary_subsimplices(Simplex_handle s) const {
    vector<Ordered_Vertex_List> result;
    int n_dim = m_pDelaunay->current_dimension();
	vector<Vertex_handle> vertices = get_vertices_of_simplex(s);
	Ordered_Vertex_List vertices2(vertices.begin(), vertices.end());
	//printf("get_boundary_subsimplices: simplex %s\n", vertices2.to_string().c_str());
	for (int i=0; i<n_dim+1; i++) {
	  // get the simplex opposite vertex i
	  Simplex_handle opposite_s = m_pDelaunay->opposite_simplex(s, i);
	  if (opposite_s == Simplex_handle()) {
	    // if it doesn't exist, then the facet opposite vertex i must be on the boundary.  add all the vertices except i to a subset		
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
	  
  void add_error_of_boundary_subsimplices(Simplex_handle const &s) {
	auto boundary_subsimplices_d = get_boundary_subsimplices(s);
	// for each size-d subset, get the d-1, d-2... 2 size subsets
	for (auto i_s2=boundary_subsimplices_d.begin(); i_s2 != boundary_subsimplices_d.end(); i_s2++) {
	  auto boundary_subsimplices_lessthan_d = i_s2->get_subsets();
	  // add each subset to the set of all boundary subsets
	  for (auto i_s3=boundary_subsimplices_lessthan_d.begin(); i_s3 != boundary_subsimplices_lessthan_d.end(); i_s3++) {
        // for each subset, calculate the error and insertion point	
        ErrorTuple err_tuple = compute_simplex_error(*i_s3);	  
	    auto pqueue_handle = m_err_pqueue.push(err_tuple);
	    m_simplex_info_map[s].m_pqueue_handles.push_back(pqueue_handle);	
        //printf ("Inserting subsimplex %x -> %s\n", Handle_hash_fn(s), err_tuple.to_str());		
	    //printf ("Inserting subset of size %d: %s\n", i_s3->size(), i_s3->to_string().c_str());	  
      }
	}
    return;
  }
  
  // taken from the implementation of Delaunay_d::Simplex_iterator
  bool is_simplex_valid(Simplex_handle const &s) const {     
      Delaunay_d::type_of_facet tf = Delaunay_d::lower_hull; 
	  Delaunay_d *DT = m_pDelaunay.get();

	  if (m_pDelaunay->current_dimension() < m_pDelaunay->dimension()) {
			return false;
	  }	  
      bool cocirc = DT->is_S_cocircular();
	  if (cocirc) {
	    if (DT->is_bounded_simplex(s)) {
 		  return true;
	    } else {
		  return false;
		}
	  } else {
	    if (DT->is_unbounded_simplex(s) && DT->type_of(s) == tf) {
		  return true;
		} else {
		  return false;
		}
	  }
	  return false;
  }  
  void clear_simplex_info(Simplex_handle const &s_handle) {
    auto map_iter = m_simplex_info_map.find(s_handle);
    if (map_iter != m_simplex_info_map.end()) {      
	  vector<Error_Priority_Queue::handle_type> &pqueue_handles(map_iter->second.m_pqueue_handles);
	  //printf ("Removing %d subsimplices of simplex %x\n", pqueue_handles.size(), Handle_hash_fn(s_handle));
      for (auto i_handle=pqueue_handles.begin(); i_handle != pqueue_handles.end(); i_handle++) {
        m_err_pqueue.erase(*i_handle);
	  }
  	  m_simplex_info_map.erase(map_iter);	    
    }	
  }
  
  void set_simplex_info(Simplex_handle const &s_handle) {
    Delaunay_d *DT = m_pDelaunay.get();
    vector<Vertex_handle> vertices = get_vertices_of_simplex(s_handle);
	vector<Point> points;
	for (int i=0; i<N+1; i++) {
	  points.push_back(vertices[i]->point());
	}
    m_simplex_info_map[s_handle].m_volume = compute_volume(N, points.begin(), points.end());
	if (m_simplex_info_map[s_handle].m_volume != 0.0) {
      m_simplex_info_map[s_handle].m_gradient = compute_simplex_gradient(s_handle);	// set gradient	    	
	  if (m_fn) {      
	    Ordered_Vertex_List vertices2(vertices.begin(), vertices.end());	
        ErrorTuple err_tuple = compute_simplex_error(vertices2);	  
	    auto pqueue_handle = m_err_pqueue.push(err_tuple);
	  
	    m_simplex_info_map[s_handle].m_pqueue_handles.push_back(pqueue_handle);

        /*		
	    printf ("Simplex %x cocirc %d volume %f bounded: %d unbounded: %d type: %d: ", Handle_hash_fn(s_handle), DT->is_S_cocircular(), m_simplex_info_map[s_handle].m_volume, DT->is_bounded_simplex(s_handle), DT->is_unbounded_simplex(s_handle), DT->type_of(s_handle));
	    for (int i=0; i<N+1; i++) {
	      std::cout << points[i] << " , ";
	    }
	    std::cout << std::endl;
	    printf ("Inserting simplex %x -> %s\n", Handle_hash_fn(s_handle), err_tuple.to_str());
        */
		
	    // handle the facets on the boundary	  
	    add_error_of_boundary_subsimplices(s_handle);
	  }
	}
  }

  vector<RT> compute_simplex_gradient(Simplex_handle const &s_handle) const {
    vector<RT> result(N);
	//std::fill(result.begin(), result.end(), std::numeric_limits<double>::quiet_NaN());
	std::fill(result.begin(), result.end(), 0.0);
	if (is_simplex_valid(s_handle)) {
	  vector< vector<RT> > vertices_coords = get_coordinates_of_vertices(s_handle);
	  vector<RT> f_vals(N+1);
	  for (int i=0; i<N+1; i++) {
	    f_vals[i] = vertices_coords[i][N];
	  }
      result = simplex_gradient<RT>(N, vertices_coords.begin(), vertices_coords.end(), f_vals.begin(), f_vals.end());	  
	}	
	return result;
  }
  
  template <class IterT>
  void insert(IterT x_begin, IterT x_end, double f) {
    assert(x_end - x_begin == N);
    insert_point(Point(N, x_begin, x_end), f);
  }
  
  void insert_point(Point const &x, double f) {    
    //std::cout << "inserting " << x << " -> " << f << std::endl;
	
	std::tuple<Simplex_handle_vec, Simplex_handle_vec> modified_and_deleted = insert_and_get_modifications(x, f);
    Simplex_handle_set modified_simplices(std::get<0>(modified_and_deleted).begin(), std::get<0>(modified_and_deleted).end());

	// remove error tuples associated with deleted and modified faces from the queue
	for (auto i_s=modified_simplices.begin(); i_s != modified_simplices.end(); i_s++) {	  	  
      clear_simplex_info(*i_s);
	}

	for (auto i_s=modified_simplices.begin(); i_s != modified_simplices.end(); i_s++) {	  	  
	  if (m_pDelaunay->current_dimension() >= N && is_simplex_valid(*i_s) && is_simplex_finite(*i_s)) {
        set_simplex_info(*i_s);
	  }
	}
	return;
  }

  std::tuple<Simplex_handle_vec, Simplex_handle_vec> insert_and_get_modifications(Point const &x, double f) {   
    std::tuple<Simplex_handle_vec, Simplex_handle_vec> result;  // first elt is modified list, second is deleted list
	// new_simplices are newly created simplices in Convex_hull_d.  They are initially unbounded (vertex 0 = anti_origin_).  Add these to our map
	// bounded_simplices are existing simplices that have vertex 0 changed from anti_origin_ to a new vertex.  remove these from our map, since they are no longer in the Delaunay triangulation
    Vertex_handle new_vertex = m_pDelaunay->insert(x, &(std::get<0>(result)) );
	m_pointMap[new_vertex] = f;	
	return result;
  }
          
  vector<Vertex_handle> get_all_vertex_handles() const {
	std::list<Vertex_handle> vertices_list = m_pDelaunay->all_vertices();
	vector<Vertex_handle> result(vertices_list.begin(), vertices_list.end());
	return result;
  }
  
  vector<Simplex_handle> get_all_simplex_handles() const {
    vector<Simplex_handle> result;
	Delaunay_d::Simplex_iterator i_s;
	for (i_s = m_pDelaunay->simplices_begin(); i_s != m_pDelaunay->simplices_end(); i_s++) {
	  result.push_back(i_s);
	}
	return result;
  }

  template <class IterT>
  double interp(IterT x_begin, IterT x_end) {
    assert(x_end - x_begin == N);
    return interp_point(Point(N, x_begin, x_end));
  }  

  double interp_point(Point const &p) const {
    double result;    	
    Simplex_handle simplex = m_pDelaunay->locate(p);	
  	if (simplex != Simplex_handle()) {		// found
	  vector<Point> points(N+1);
	  for (int i=0; i<N+1; i++) {
	    Vertex_handle vertex = m_pDelaunay->vertex_of_simplex(simplex, i);		
	    points[i] = Point(N, vertex->point().cartesian_begin(), vertex->point().cartesian_begin() + N);  // get rid of the extra dimension
	  }
	  auto barycentric_coords = compute_barycentric_coords(points.begin(), points.end(), p);
	  RT sum = 0.0;
	  for (int i=0; i<barycentric_coords.size(); i++) {
	    Vertex_handle vertex = m_pDelaunay->vertex_of_simplex(simplex, i);
	    double val = m_pointMap.at(vertex);
		sum += barycentric_coords[i] * val;		  
	  }
      result = RT_to_double(sum);		
	} else {
	  //bpl_assert(false, "interp: no simplex containing the point could be found");
      result = std::numeric_limits<double>::quiet_NaN();	
	}
    return result;
  }
  
  // return a vector of vertices of simplex s
  vector<Vertex_handle> get_vertices_of_simplex(Simplex_handle s, bool filter_empty=true) const {
    vector<Vertex_handle> result;
	int n = m_pDelaunay->current_dimension();
	for (int i=0; i<n+1; i++) {
	  Vertex_handle v = m_pDelaunay->vertex_of_simplex(s, i);
	  if (filter_empty == true && v != Vertex_handle()) {
	    result.push_back(v);
	  }
	}
	return result;
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
  
  vector<ErrorTuple> get_error_queue() const {
    vector<ErrorTuple> result(m_err_pqueue.ordered_begin(), m_err_pqueue.ordered_end());
	return result;
  }
 

  // given a sequence of vertices, return the coordinates of the points, plus the value of f(x) as an additional coordinate
  template <class IterT>
  vector< vector<RT> > get_coordinates_of_vertices(IterT vertices_begin, IterT vertices_end) const {    
    int n_points = vertices_end - vertices_begin;
    vector< vector<RT> > result(n_points);
	for (int i=0; i<n_points; i++) {
	  Vertex_handle v = vertices_begin[i];
	  result[i].resize(N+1);
	  std::copy(v->point().cartesian_begin(), v->point().cartesian_begin()+N, result[i].begin());
	  result[i][N] = RT(m_pointMap.at(v));		// last dimension is the f value
	}
	return result;
  }
  vector< vector<RT> > get_coordinates_of_vertices(Simplex_handle const &s_handle) const {
    vector<Vertex_handle> vertices;
    for (int i=0; i<N+1; i++) {	  
	  vertices.push_back(m_pDelaunay->vertex_of_simplex(s_handle, i));
	}
	return get_coordinates_of_vertices(vertices.begin(), vertices.end());
  }

  bool is_simplex_finite(Simplex_handle const &s_handle) const {    
	Vertex_handle vh;
    for (int i=0; i<N+1; i++) {	  
	  if (m_pDelaunay->vertex_of_simplex(s_handle, i) == vh) {
	    return false;
	  }
	}
	return true;
  }
  
  // compute the error of a simplex of any dimension (e.g. boundary facets)
  ErrorTuple compute_simplex_error(Ordered_Vertex_List const &vertices) const {      
	typedef vector<RT> point_t;    // one additional dimension for the value of f(x)
	vector<point_t> vertex_points = get_coordinates_of_vertices(vertices.begin(), vertices.end());   // get the coordinates of the vertices	    	
	point_t centroid(N+1);
	get_centroid(N, vertex_points.begin(), vertex_points.end(), centroid.begin());	// get centroid
	double f_val = eval_fn_at_coords(centroid.begin(), centroid.begin() + N);	// evaluate fn at centroid	
	centroid[N] = RT(f_val);	
    // find distance from p0 to best-fit plane
    double sqdist;
    point_t closest_point(N+1);
    sqdist = sqdist_point_hplane_of_points<RT>(N+1, centroid.begin(), vertex_points.begin(), vertex_points.end(), closest_point.begin());			
	// remove the last coordinate from centroid_coords. this is where a new grid point will be added if this simplex were chosen
	Point centroid_point(N, centroid.begin(), centroid.begin()+N);	
	ErrorTuple result(sqdist, f_val, centroid_point);
	return result;
  }  
    
  //////////////////////////////////////////////////////////////////
  // Plotting-related functions
  vector< vector<RT> > get_all_vertices() const {
    std::list<Vertex_handle> vertex_list = m_pDelaunay->all_vertices();
	Vertex_handle_vec vec2(vertex_list.begin(), vertex_list.end());
	return get_coordinates_of_vertices(vec2.begin(), vec2.end());
  }
  
  vector<std::tuple<Vertex_handle, Vertex_handle> > get_line_segments_of_simplex(Simplex_handle s) const {
    int n_dims = m_pDelaunay->dimension();
    vector<std::tuple<Vertex_handle, Vertex_handle> > result(n_dims + 1);
	for (int i=0; i<n_dims+1; i++) {
	  Vertex_handle v1 = m_pDelaunay->vertex_of_simplex(s, i);
	  Vertex_handle v2 = m_pDelaunay->vertex_of_simplex(s, (i+1)%(n_dims+1));
	  result[i] = std::make_tuple(v1, v2);
	}
	return result;
  }
  
  // for plotting. return a vector of all line segments, and the min/max at each coordinate
  std::tuple<vector<std::tuple<vector<RT>, vector<RT>, int> >, vector<double>, vector<double> > get_line_segments() const {    
    Ordered_Vertex_List_set all_segments;
	std::list<Simplex_handle> all_simplices = m_pDelaunay->all_simplices();
	for (auto i_s=all_simplices.begin(); i_s != all_simplices.end(); i_s++) {
	  auto segments_of_s = get_line_segments_of_simplex(*i_s);
	  for (auto i_seg=segments_of_s.begin(); i_seg != segments_of_s.end(); i_seg++) {
	    array<Vertex_handle,2> seg = {std::get<0>(*i_seg), std::get<1>(*i_seg)};	    
	    all_segments.insert(Ordered_Vertex_List(seg.begin(), seg.end()));
	  }
	}
	vector<std::tuple< vector<RT>, vector<RT>, int> > result;
	vector<double> min_coords(N+1), max_coords(N+1);
	std::fill(min_coords.begin(), min_coords.end(), DBL_MAX);
	std::fill(max_coords.begin(), max_coords.end(), -DBL_MAX);
	for (auto i_seg=all_segments.begin(); i_seg != all_segments.end(); i_seg++) {	  
	  vector< vector<RT> > vertices_with_f = get_coordinates_of_vertices(i_seg->begin(), i_seg->end());
	  result.push_back(std::make_tuple(vertices_with_f[0], vertices_with_f[1], 1));
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

// return the N+1 dimensional lifted point
Point get_lifted_point(boost::shared_ptr<Delaunay_d> pDelaunay, Point const &x) { 
  Delaunay_d::Lifted_R::Lift_to_paraboloid_d lift = pDelaunay->lifted_kernel().lift_to_paraboloid_d_object();
  Delaunay_d::Lifted_point_d lp = lift(x);
  return lp;
}

                                         
