#include <boost/python.hpp>
#include <boost/python/operators.hpp>
#include "delaunay_d_interp.h"

#include <Python.h>
#include <numpy/arrayobject.h>
#include <pyublas/numpy.hpp>

namespace bpl = boost::python;

#define bpl_assert(exp, s)	\
	if (!(exp)) {				\
      PyErr_SetString(PyExc_ValueError, (s));	\
      boost::python::throw_error_already_set();	\
	}	

typedef pyublas::numpy_vector<double> dPyArr;
typedef vector<dPyArr> dPyArrVector;

/////////////////////////////////////////////////////////
// convert to/from Python types

// basic conversions:
// dPyArr, dVec -> vector<RT>
// dPyArr -> Point
// vector<dPyArr> -> vector<Point>
// double -> RT
// Point -> doubles -> dPyArr
template <class IterT>
vector<RT> doubles_to_RTvec(IterT doubles_begin, IterT doubles_end) {
  int n = doubles_end - doubles_begin;
  vector<RT> result(n);
  std::transform(doubles_begin, doubles_end, result.begin(), [] (double x) -> RT { return RT(x); } );
  return result;
}

Point coords_to_point(dPyArr const &x) {
  bpl_assert(x.ndim() == 1, "x should have dimension 1");
  int n_dims = x.dims()[0];
  vector<RT> coords(n_dims);
  for (int i=0; i<n_dims; i++) {
    coords[i] = RT(x.sub(i));
  }
  return Point(n_dims, coords.begin(), coords.end());
}

template <class IterT>
Point coords_to_point(IterT coords_begin, IterT coords_end) {
  int n_dims = (coords_end - coords_begin);  
  vector<RT> coords = doubles_to_RTvec(coords_begin, coords_end);
  return Point(n_dims, coords.begin(), coords.end());
}

template <class QuotientT>
struct Quotient_to_tuple {
  static PyObject* convert(QuotientT const &x) {    
	//bpl::tuple result = bpl::make_tuple(to_double(x.num), to_double(x.den));	
	bpl::object result(to_double(x));
    return boost::python::incref(result.ptr());
  }                                            
}; 

template <class IterT>
vector<Point> coordList_to_PointVec(IterT coords_begin, IterT coords_end) {
  int n_points = coords_end - coords_begin;
  vector<Point> result(n_points);
  std::transform(coords_begin, coords_end, result.begin(), [] (dPyArr const &x) -> Point { return coords_to_point(x); } );
  return result;
}
  
struct Point_to_list { 
  static PyObject* convert(Point const &x) {
    int n_dims = x.dimension();
	bpl::list result;
	for (int i=0; i<n_dims; i++) {
	  result.append(RT_to_double(x[i]));
	}
    return boost::python::incref(result.ptr());
  }                                            
}; 
  
dPyArr point_to_dPyArr(Point const &x) {
  int n_dims = x.dimension();
  dPyArr result(n_dims);
  std::transform(x.cartesian_begin(), x.cartesian_end(), result.begin(), [] (RT const &c) -> double { return RT_to_double(c); } );
  return result;
}

template <class IterableT>
struct forward_iterable_to_list { 
  static PyObject* convert(IterableT const &cont) {
    bpl::list result;
    for (auto iter=cont.begin(); iter != cont.end(); iter++) {
      result.append(*iter);
    }                      
    return boost::python::incref(result.ptr());
  }                                            
}; 

template <class IterableT>
struct forward_iterable_to_double_array { 
  static PyObject* convert(IterableT const &cont) {
    dPyArr result(cont.size());
	auto iter=cont.begin();
	auto iter2=result.begin();
    for (; iter != cont.end(); iter++, iter2++) {
	  (*iter2) = RT_to_double(*iter);
    }          
    bpl::object result2(result.to_python());
    return boost::python::incref(result2.ptr());
  }                                            
}; 

template <class MapT>
struct map_to_list {
  static PyObject* convert(MapT const &map) {
    bpl::list result;
	for (auto iter=map.begin(); iter != map.end(); iter++) {
	  result.append(bpl::make_tuple(iter->first, iter->second));
	}
    return boost::python::incref(result.ptr());
  }                                            
}; 

template <class T1, class T2>        
struct tuple_2_to_python { 
  static PyObject* convert(std::tuple<T1, T2> const &arg) {
    bpl::tuple result = bpl::make_tuple(std::get<0>(arg), std::get<1>(arg));
    return boost::python::incref(result.ptr());
  }                                            
}; 

template <class T1, class T2, class T3>        
struct tuple_3_to_python { 
  static PyObject* convert(std::tuple<T1, T2, T3> const &arg) {
    bpl::tuple result = bpl::make_tuple(std::get<0>(arg), std::get<1>(arg), std::get<2>(arg));
    return boost::python::incref(result.ptr());
  }                                            
}; 

///////////////////////////////////////////////////////////////////////////////////////
// convert from Python types

// Python list of doubles -> Point
struct Point_from_python_list {
  static void *check(PyObject *pObj) {    
    if (!PyList_Check(pObj)) {
	  return NULL;
	} else {
	  bpl::handle<> h(bpl::borrowed(pObj));
	  bpl::list l(h);
	  for (int i=0; i<bpl::len(l); i++) {
	    bpl::object obj(l[i]);		
	    if (!PyFloat_Check(obj.ptr())) {
		  return NULL;
		}
	  }
	}
	return pObj;
  }
  static void construct(PyObject *pObj, boost::python::converter::rvalue_from_python_stage1_data* data) {
	bpl::handle<> h(bpl::borrowed(pObj));
	bpl::list l(h);
	void *storage = ((boost::python::converter::rvalue_from_python_storage<Point>*)data)->storage.bytes;
	// copy elements
    vector<RT> coords(bpl::len(l));
    for (int i=0; i<bpl::len(l); i++) {
      coords[i] = RT(bpl::extract<double>(l[i]));
    }
    new (storage) Point(bpl::len(l), coords.begin(), coords.end());
	Point *p = static_cast<Point*>(storage);
	data->convertible = storage;
	//std::cout << Point_to_string(*p) << std::endl;
  }    
};

struct ErrorTuple_to_python { 
  static PyObject* convert(ErrorTuple const &arg) {
    bpl::tuple result = bpl::make_tuple(arg.m_err, arg.m_f, arg.m_x);
    return boost::python::incref(result.ptr());
  }                                            
}; 

///////////////////////////////////////////////////////////////////////////////////////////
// wrapper for math functions

/*
dPyArr best_fit_plane_wrap(vector<dPyArr> const &point_list) {
  vector<Point> points(point_list.size());
  int n_dims = point_list[0].dims()[0];
  //printf("best_fit_plane_wrap: %d points of %d-d points\n", point_list.size(), n_dims);
  for (int i=0; i<point_list.size(); i++) {
	bpl_assert(point_list[0].dims()[0] == n_dims, "points are not all the same size");
    points[i] = coords_to_point(point_list[i]);
  }
  Hyperplane hplane = best_fit_plane(points.begin(), points.end());
  dPyArr result(hplane.dimension() + 1);
  for (int i=0; i<result.size(); i++) {
    result[i] = RT_to_double(hplane.coefficient(i));
  }
  return result;
}
*/

/*
bpl::tuple sqdist_point_hyperplane_wrap(dPyArr const &hplane_coefs, dPyArr const &x) {
  Point p = coords_to_point(x);
  int n_dims = x.size();
  bpl_assert(hplane_coefs.size() == n_dims+1, "hyperplane must have size == n_dims+1");
  vector<RT> coefs(&hplane_coefs[0], &hplane_coefs[0] + hplane_coefs.size());
  //std::transform(&hplane_coefs[0], &hplane_coefs[0] + hplane_coefs.size(), coefs.begin(), [] (double d) -> RT { return RT(d) } );
  Hyperplane hplane(n_dims, coefs.begin(), coefs.end());
  RT sqdist;
  Point closest_point;
  std::tie(sqdist, closest_point) = sqdist_point_hyperplane(hplane, p);
  return bpl::make_tuple(RT_to_double(sqdist), point_to_dPyArr(p));
}
*/

class Delaunay_incremental_interp_d_wrap : public Delaunay_incremental_interp_d {
public:
  typedef Delaunay_incremental_interp_d super;
  bpl::object m_fn_obj;
  Delaunay_incremental_interp_d_wrap(int n_dims, bpl::object fn) 
  : Delaunay_incremental_interp_d(n_dims)
  {
    if (fn != bpl::object()) {
	  super::m_fn = [=](int n_args, double *args_begin)->double { return this->call_python_fn(n_args, args_begin); };  
      m_fn_obj = fn;	  
    }
  }
    
  double call_python_fn(int n_args, double *args_begin) const {
    bpl::list args;    
	for (int i=0; i<n_args; i++) {
	  args.append(args_begin[i]);
	}
    bpl::object result = m_fn_obj(args);
    double d_result = bpl::extract<double>(result);
    return d_result;
  }
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(insert_overloads, Delaunay_d::insert, 1, 2)

BOOST_PYTHON_MODULE(_delaunay_d_python)
{ 
  bpl::to_python_converter<CGAL::Quotient<NT>, Quotient_to_tuple<CGAL::Quotient<NT>> >();  
  bpl::to_python_converter<Point, Point_to_list>();
  bpl::to_python_converter<std::tuple<Point, Point>, tuple_2_to_python<Point, Point> >();
  bpl::to_python_converter<vector<RT>, forward_iterable_to_double_array<vector<RT> > >();
  bpl::to_python_converter<vector<double>, forward_iterable_to_double_array<vector<double> > >();
  bpl::to_python_converter<vector<Point>, forward_iterable_to_list<vector<Point> > >();
  bpl::to_python_converter<vector<ErrorTuple>, forward_iterable_to_list<vector<ErrorTuple> > >();
  bpl::to_python_converter<vector<Vertex_handle>, forward_iterable_to_list<vector<Vertex_handle> > >();  
  bpl::to_python_converter<vector<Simplex_handle>, forward_iterable_to_list<vector<Simplex_handle> > >();  
  //bpl::to_python_converter<vector<Facet_handle>, forward_iterable_to_list<vector<Facet_handle> > >();  
  
  bpl::to_python_converter<std::list<Vertex_handle>, forward_iterable_to_list<std::list<Vertex_handle> > >();  
  bpl::to_python_converter<std::list<Simplex_handle>, forward_iterable_to_list<std::list<Simplex_handle> > >();    
  //bpl::to_python_converter<std::list<Facet_handle>, forward_iterable_to_list<std::list<Facet_handle> > >();    
  //bpl::to_python_converter<ErrorTuple, tuple_3_to_python<RT, Point, double> >();  
  bpl::to_python_converter<ErrorTuple, ErrorTuple_to_python >();  
  bpl::to_python_converter<InsertResult, tuple_2_to_python<Vertex_handle, vector<Simplex_handle> > >();
  bpl::to_python_converter<Vertex_set, forward_iterable_to_list<Vertex_set> >();  
  bpl::to_python_converter<Simplex_handle_set, forward_iterable_to_list<Simplex_handle_set> >();  
  bpl::to_python_converter<std::tuple<Simplex_handle_set, Vertex_set>, tuple_2_to_python<Simplex_handle_set, Vertex_set> >();

  typedef std::tuple<vector<RT>, vector<RT>, int> point_point_int_tuple;
  bpl::to_python_converter<point_point_int_tuple, tuple_3_to_python<vector<RT>, vector<RT>, int>>();
  bpl::to_python_converter<vector<point_point_int_tuple>, forward_iterable_to_list<vector<point_point_int_tuple> >>();  
  bpl::to_python_converter<std::tuple<vector<point_point_int_tuple>, dVec, dVec>, tuple_3_to_python<vector<point_point_int_tuple>, dVec, dVec> >();
  
  bpl::converter::registry::push_back(&Point_from_python_list::check, &Point_from_python_list::construct, bpl::type_id<Point >());
  
  //boost::python::def("best_fit_plane", best_fit_plane_wrap);
  //boost::python::def("sqdist_point_hyperplane", sqdist_point_hyperplane_wrap);

  bpl::enum_<CGAL::Bounded_side>("Bounded_side")
        .value("ON_UNBOUNDED_SIDE", CGAL::ON_UNBOUNDED_SIDE)
        .value("ON_BOUNDARY", CGAL::ON_BOUNDARY)
		.value("ON_BOUNDED_SIDE", CGAL::ON_BOUNDED_SIDE)
        ;  
		
  bpl::enum_<Delaunay_d::Delaunay_voronoi_kind>("Delaunay_voronoi_kind")
        .value("NEAREST", Delaunay_d::NEAREST)
        .value("FURTHEST", Delaunay_d::FURTHEST)
        ;  

  bpl::class_<Vertex_handle>("Vertex_handle", bpl::init<>())      
		.def("__hash__", Handle_hash_fn<Vertex_handle>)
		.def(bpl::self == bpl::self)
    ;

  bpl::class_<Simplex_handle>("Simplex_handle", bpl::init<>())
		.def("__hash__", Handle_hash_fn<Simplex_handle>)
		.def(bpl::self == bpl::self)
    ;

  // bpl::class_<Facet_handle>("Facet_handle", bpl::init<>())
		// .def("__hash__", Handle_hash_fn<Facet_handle>)
		// .def(bpl::self == bpl::self)
    // ;
	
  bpl::class_<Convex_hull_d, boost::shared_ptr<Convex_hull_d>, boost::noncopyable>("Convex_hull_d", bpl::init<int>())
		.def("dimension", &Convex_hull_d::dimension)		
		.def("current_dimension", &Convex_hull_d::current_dimension)		
		.def("associated_point", static_cast< Point(Convex_hull_d::*)(Vertex_handle) const> ( &Convex_hull_d::associated_point ))		
		.def("vertex_of_simplex", static_cast< Vertex_handle(Convex_hull_d::*)(Simplex_handle, int) const> ( &Convex_hull_d::vertex_of_simplex ))
		.def("point_of_simplex", static_cast< Point(Convex_hull_d::*)(Simplex_handle, int) const> ( &Convex_hull_d::point_of_simplex ))
		.def("opposite_simplex", static_cast< Simplex_handle(Convex_hull_d::*)(Simplex_handle, int) const> ( &Convex_hull_d::opposite_simplex ))
		.def("simplex", static_cast< Simplex_handle(Convex_hull_d::*)(Vertex_handle) const> ( &Convex_hull_d::simplex ))
		.def("index", static_cast< int(Convex_hull_d::*)(Vertex_handle) const> ( &Convex_hull_d::index))
		.def("vertex_of_facet", static_cast< Vertex_handle(Convex_hull_d::*)(Facet_handle, int) const> ( &Convex_hull_d::vertex_of_facet ))
		.def("point_of_facet", static_cast< Point(Convex_hull_d::*)(Facet_handle, int) const> ( &Convex_hull_d::point_of_facet ))
		.def("opposite_facet", static_cast< Facet_handle(Convex_hull_d::*)(Facet_handle, int) const> ( &Convex_hull_d::opposite_facet ))
		.def("index_of_vertex_in_opposite_facet", static_cast< int(Convex_hull_d::*)(Facet_handle, int) const> ( &Convex_hull_d::index_of_vertex_in_opposite_facet ))
		.def("insert", &Convex_hull_d::insert)
		.def("is_dimension_jump", &Convex_hull_d::is_dimension_jump)
		.def("facets_visible_from", &Convex_hull_d::facets_visible_from)
		.def("bounded_side", &Convex_hull_d::bounded_side)
		.def("clear", &Convex_hull_d::clear)
		.def("number_of_vertices", &Convex_hull_d::number_of_vertices)
		.def("number_of_facets", &Convex_hull_d::number_of_facets)
		.def("print_statistics", &Convex_hull_d::print_statistics)
		
		.def("index_of_vertex_in_opposite_simplex", static_cast< int(Convex_hull_d::*)(Simplex_handle, int) const> ( &Convex_hull_d::index_of_vertex_in_opposite_simplex ))
		
		.def("is_unbounded_simplex", static_cast< bool(Convex_hull_d::*)(Simplex_handle) const> ( &Convex_hull_d::is_unbounded_simplex ))
		.def("is_bounded_simplex", static_cast< bool(Convex_hull_d::*)(Simplex_handle) const> ( &Convex_hull_d::is_bounded_simplex ))		
		
		//.def("all_points", static_cast< std::list<Point>& (Convex_hull_d::*)() const> (&Convex_hull_d::all_points ))		
		.def("all_vertices", static_cast< std::list<Vertex_handle>(Convex_hull_d::*)()> (&Convex_hull_d::all_vertices ))		
        .def("all_simplices", static_cast< std::list<Simplex_handle>(Convex_hull_d::*)() > (&Convex_hull_d::all_simplices ))		
        .def("all_facets", static_cast< std::list<Facet_handle>(Convex_hull_d::*)() > (&Convex_hull_d::all_facets ))		
		
	;

  bpl::class_<Delaunay_d, bpl::bases<Convex_hull_d>, boost::shared_ptr<Delaunay_d>, boost::noncopyable>("Delaunay_d", bpl::init<int>())
		.def("dimension", &Delaunay_d::dimension)		
		.def("current_dimension", &Delaunay_d::current_dimension)
		.def("is_simplex_of_nearest", static_cast< bool(Delaunay_d::*)(Simplex_handle) const> ( &Delaunay_d::is_simplex_of_nearest ))
		.def("is_simplex_of_furthest", static_cast< bool(Delaunay_d::*)(Simplex_handle) const> ( &Delaunay_d::is_simplex_of_furthest ))
		.def("vertex_of_simplex", static_cast< Vertex_handle(Delaunay_d::*)(Simplex_handle, int) const> ( &Delaunay_d::vertex_of_simplex ))
		.def("associated_point", static_cast< Point(Delaunay_d::*)(Vertex_handle) const> ( &Delaunay_d::associated_point ))
		.def("point_of_simplex", static_cast< Point(Delaunay_d::*)(Simplex_handle, int) const> ( &Delaunay_d::point_of_simplex ))
		.def("opposite_simplex", static_cast< Simplex_handle(Delaunay_d::*)(Simplex_handle, int) const> ( &Delaunay_d::opposite_simplex ))
		.def("index_of_vertex_in_opposite_simplex", static_cast< int(Delaunay_d::*)(Simplex_handle, int) const> ( &Delaunay_d::index_of_vertex_in_opposite_simplex ))
		.def("simplex", &Delaunay_d::simplex)		
		.def("index", &Delaunay_d::index)		
		.def("contains", &Delaunay_d::contains)		
		.def("empty", &Delaunay_d::empty)		
		.def("clear", &Delaunay_d::clear)		
		.def("insert", &Delaunay_d::insert, insert_overloads())		
		.def("locate", &Delaunay_d::locate)		
		.def("lookup", &Delaunay_d::lookup)		
		.def("nearest_neighbor", &Delaunay_d::nearest_neighbor)			
		.def("all_simplices", &Delaunay_d::all_simplices)		
		.def("all_vertices", &Delaunay_d::all_vertices)
    ;
	
  bpl::class_<Delaunay_incremental_interp_d_wrap, boost::noncopyable>("Delaunay_incremental_interp_d", bpl::init<int, bpl::object>())
		.def_readonly("triang", &Delaunay_incremental_interp_d::m_pDelaunay)
		.def_readonly("fn", &Delaunay_incremental_interp_d_wrap::m_fn_obj)		
		.def("insert", &Delaunay_incremental_interp_d_wrap::insert_point)		
		.def("interp", &Delaunay_incremental_interp_d_wrap::interp_point)
		.def("get_error_queue", &Delaunay_incremental_interp_d_wrap::get_error_queue)
		.def("get_largest_error_tuple", &Delaunay_incremental_interp_d_wrap::get_largest_error_tuple)
		.def("insert_largest_error_point", &Delaunay_incremental_interp_d_wrap::insert_largest_error_point)
		.def("get_line_segments", &Delaunay_incremental_interp_d_wrap::get_line_segments)
		.def("get_all_vertices", &Delaunay_incremental_interp_d_wrap::get_all_vertices)
		//.def("gradient", &Delaunay_incremental_interp_d_wrap::gradient)		

		//.def("insert_point", &Delaunay_incremental_interp_d::insert_point)
		// .def("insert_point_no_error_update", &Delaunay_incremental_interp_d::insert_point_no_error_update)
		// .def("locate", &Delaunay_incremental_interp_d::locate)
		// .def("all_vertices_with_f", &Delaunay_incremental_interp_d::all_vertices_with_f)
		// .def("all_vertices", &Delaunay_incremental_interp_d::all_vertices)
		// .def("all_simplices", &Delaunay_incremental_interp_d::all_simplices)
		// .def("interp", &Delaunay_incremental_interp_d::interp)
		// .def("interp_vec", &Delaunay_incremental_interp_d::interp_vec)
		// .def("get_vertices_of_simplex", &Delaunay_incremental_interp_d::get_vertices_of_simplex)
		// .def("compute_simplex_error2", &Delaunay_incremental_interp_d::compute_simplex_error2)
		// .def("get_largest_error_tuple", &Delaunay_incremental_interp_d::get_largest_error_tuple)
		// .def("get_error_queue", &Delaunay_incremental_interp_d::get_error_queue)	
		// .def("init_error_queue", &Delaunay_incremental_interp_d::init_error_queue)	
		// .def("get_val_of_vertex", &Delaunay_incremental_interp_d::get_val_of_vertex)		
		// .def("get_line_segments", &Delaunay_incremental_interp_d::get_line_segments)		
	;	  
	boost::python::def("get_lifted_point", get_lifted_point);
} 