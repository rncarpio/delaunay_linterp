
#include <boost/python.hpp>
#include <boost/python/operators.hpp>

#include "delaunay_2_interp.h"

////////////////////////////////////////////////////////////////////////////////////
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

struct RT_to_double_struct {
  static PyObject* convert(CGAL::Lazy_exact_nt<class CGAL::Gmpq> const &x) {    	
	bpl::object result(RT_to_double(x));
    return boost::python::incref(result.ptr());
  }                                            
}; 

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

Point coords_to_point(dPyArr const &x) {
  bpl_assert(x.ndim() == 1, "x should have dimension 1");
  int n_dims = x.dims()[0];
  vector<RT> coords(n_dims);
  for (int i=0; i<n_dims; i++) {
    coords[i] = RT(x.sub(i));
  }
  return Point(coords[0], coords[1]);
}

template <class QuotientT>
struct Quotient_to_tuple {
  static PyObject* convert(QuotientT const &x) {    
	//bpl::tuple result = bpl::make_tuple(to_double(x.num), to_double(x.den));	
	bpl::object result(to_double(x));
    return boost::python::incref(result.ptr());
  }                                            
}; 

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

// Python list of doubles -> Point
template <int N>
struct Point_from_python_list {
  static void *check(PyObject *pObj) {    
    if (!PyList_Check(pObj)) {
	  return NULL;
	} else {
	  bpl::handle<> h(bpl::borrowed(pObj));
	  bpl::list l(h);
	  if (bpl::len(l) != N) {
	    return NULL;
	  }
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
	bpl_assert(bpl::len(l) == N, "must be length N");
    vector<RT> coords(bpl::len(l));
    for (int i=0; i<bpl::len(l); i++) {
	  double x = bpl::extract<double>(l[i]);
      coords[i] = x;
    }
    new (storage) Point(coords[0], coords[1]);
	Point *p = static_cast<Point*>(storage);
	data->convertible = storage;
	//std::cout << Point_to_string(*p) << std::endl;
  }    
};

BOOST_PYTHON_MODULE(_delaunay_interp2)
{   
  bpl::converter::registry::push_back(&Point_from_python_list<2>::check, &Point_from_python_list<2>::construct, bpl::type_id<Point >());
  bpl::to_python_converter<CGAL::Lazy_exact_nt<class CGAL::Gmpq>, RT_to_double_struct >();
  bpl::to_python_converter<vector<ErrorTuple>, forward_iterable_to_list<vector<ErrorTuple> > >();
  bpl::to_python_converter<array<RT,2>, forward_iterable_to_list<array<RT,2> > >();
  bpl::to_python_converter<array<RT,3>, forward_iterable_to_list<array<RT,3> > >();
  bpl::to_python_converter<ErrorTuple, ErrorTuple_to_python >();
  bpl::to_python_converter<Point, Point_to_list>();
  bpl::to_python_converter<std::tuple<RT_3, RT_3, int>, tuple_3_to_python<RT_3, RT_3, int> >();
  bpl::to_python_converter<vector<std::tuple<RT_3, RT_3, int>>, forward_iterable_to_list<vector<std::tuple<RT_3, RT_3, int>> > >();  
  bpl::to_python_converter<std::tuple<Point, int>, tuple_2_to_python<Point, int> >();
  bpl::to_python_converter<vector<std::tuple<Point, int>>, forward_iterable_to_list<vector<std::tuple<Point, int>> > >();
  
  bpl::to_python_converter<std::tuple<vector<std::tuple<RT_3, RT_3, int> >, dPyArr, dPyArr>, tuple_3_to_python<vector<std::tuple<RT_3, RT_3, int> >, dPyArr, dPyArr> >();
  
  bpl::enum_<CGAL::Bounded_side>("Bounded_side")
        .value("ON_UNBOUNDED_SIDE", CGAL::ON_UNBOUNDED_SIDE)
        .value("ON_BOUNDARY", CGAL::ON_BOUNDARY)
		.value("ON_BOUNDED_SIDE", CGAL::ON_BOUNDED_SIDE)
        ;  

  bpl::class_<Delaunay_incremental_interp_2, boost::noncopyable>("DelaunayInterp2", bpl::init<bpl::object>())
		.def("insert_only", &Delaunay_incremental_interp_2::insert_only)
		.def("insert_with_update", &Delaunay_incremental_interp_2::insert_and_update_error_queue)		
		.def("interp", &Delaunay_incremental_interp_2::interp)
		.def("get_error_queue", &Delaunay_incremental_interp_2::get_error_queue)
		.def("get_largest_error_tuple", &Delaunay_incremental_interp_2::get_largest_error_tuple)
		.def("insert_largest_error_point", &Delaunay_incremental_interp_2::insert_largest_error_point)
		.def("get_line_segments", &Delaunay_incremental_interp_2::get_line_segments)
		.def("get_all_vertices", &Delaunay_incremental_interp_2::get_all_vertices)
		.def("gradient", &Delaunay_incremental_interp_2::gradient)
		;

}

