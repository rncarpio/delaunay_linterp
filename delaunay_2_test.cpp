
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

#include <stdio.h>
#include <ctime>
#include "delaunay_2_interp.h"

dVec linspace(double first, double last, int len) {
  dVec result(len);
  double step = (last-first) / (len - 1);
  for (int i=0; i<len; i++) {
    result[i] = first + i*step;
  }
  return result;
}

template <int N>
class CobbDouglas {
public:
  array<double,N> m_alphas;

  template <class IterT>
  CobbDouglas(IterT alpha_begin, IterT alpha_end) {
    assert(alpha_end - alpha_begin == N);
    std::copy(alpha_begin, alpha_end, m_alphas.begin());
  }
  template <class IterT>
  double operator() (IterT args_begin) const {
    double result = 1.0;
    for (int i=0; i<N; i++) {
      result *= (pow(args_begin[i], m_alphas[i]));
    }
    return result;
  }
};

int main(int argc, char **argv) {
  int n_points = 10;
  clock_t t1, t2;
  Delaunay_incremental_interp_2 triang;
  array<double,2> alphas = {0.3, 0.6};
  CobbDouglas<2> cd_fn(alphas.begin(), alphas.end());
  array<double,2> args;
  double f_val;

  // use rectangular grid
  dVec grid1 = linspace(0.01, 5, n_points);
  dVec grid2 = linspace(0.01, 5, n_points);
  t1 = clock();	
  for (int i=0; i<n_points; i++) {
    for (int j=0; j<n_points; j++) {
      args[0] = grid1[i];
      args[1] = grid2[j];
      f_val = cd_fn(args.begin());
      triang.insert(args.begin(), args.end(), f_val);
    }
  }
  t2 = clock();
  printf("%d insertions, %d clocks, %f sec\n", n_points*n_points, (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  // compare vs. actual function
  dVec interp_grid1 = linspace(0.01, 5, 2*n_points);
  dVec interp_grid2 = linspace(0.01, 5, 2*n_points);
  dVec true_f_vals;
  for (int i=0; i<interp_grid1.size(); i++) {
    for (int j=0; j<interp_grid2.size(); j++) {
	  args[0] = interp_grid1[i];
      args[1] = interp_grid2[j];
	  true_f_vals.push_back(cd_fn(args.begin()));
	}
  }
  dVec interp_vals;
  t1 = clock();	  
  for (int i=0; i<interp_grid1.size(); i++) {
    for (int j=0; j<interp_grid2.size(); j++) {
	  args[0] = interp_grid1[i];
      args[1] = interp_grid2[j];	
	  interp_vals.push_back(triang.interp(args.begin(), args.end()));
	}
  }
  t2 = clock();
  printf("%d interpolations, %d clocks, %f sec\n", interp_vals.size(), (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);
  double sse = 0.0, diff;
  for (int i=0; i<interp_vals.size(); i++) {
    diff = true_f_vals[i] - interp_vals[i];
    sse += diff*diff;
  }
  printf("sum of squared errors: %f\n", sse);
  
  // use adaptive grid
  printf("Adaptive grid:\n");
  std::function<double(int, double*)> fn = [&cd_fn] (int n, double *x_begin)->double { return cd_fn(x_begin); };
  Delaunay_incremental_interp_2 adaptive_triang(fn);
  // insert boundary points
  args[0] = grid1.front(); args[1] = grid2.front(); f_val = cd_fn(args.begin()); adaptive_triang.insert(args.begin(), args.end(), f_val);
  args[0] = grid1.front(); args[1] = grid2.back(); f_val = cd_fn(args.begin()); adaptive_triang.insert(args.begin(), args.end(), f_val);
  args[0] = grid1.back(); args[1] = grid2.front(); f_val = cd_fn(args.begin()); adaptive_triang.insert(args.begin(), args.end(), f_val);
  args[0] = grid1.back(); args[1] = grid2.back(); f_val = cd_fn(args.begin()); adaptive_triang.insert(args.begin(), args.end(), f_val);  
  for (int i=0; i<n_points*n_points; i++) {
    adaptive_triang.insert_largest_error_point();
  }	
  
  interp_vals.clear();
  t1 = clock();	  
  for (int i=0; i<interp_grid1.size(); i++) {
    for (int j=0; j<interp_grid2.size(); j++) {
	  args[0] = interp_grid1[i];
      args[1] = interp_grid2[j];	
	  interp_vals.push_back(adaptive_triang.interp(args.begin(), args.end()));
	}
  }
  t2 = clock();
  printf("%d interpolations, %d clocks, %f sec\n", interp_vals.size(), (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);
  sse = 0.0;
  for (int i=0; i<interp_vals.size(); i++) {
    diff = true_f_vals[i] - interp_vals[i];
    sse += diff*diff;
  }
  printf("sum of squared errors: %f\n", sse);
  
  return 0;
}

