
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

#include <ctime>
#include "delaunay_d_interp.h"

// return an evenly spaced 1-d grid of doubles.
std::vector<double> linspace(double first, double last, int len) {
  std::vector<double> result(len);
  double step = (last-first) / (len - 1);
  for (int i=0; i<len; i++) { result[i] = first + i*step; }
  return result;
}

// the function to interpolate.
double fn (int, double* x_begin) { return pow(x_begin[0],0.3) * pow(x_begin[1],0.6); }

int main(int argc, char **argv) {
  int n_points = 10;
  clock_t t1, t2;      

  char pcTemp[1024];
  gets(pcTemp);
  
  // first, try a rectangular grid
  Delaunay_incremental_interp_d triang(2);
  std::vector<double> grid = linspace(0.01, 5., n_points);
  t1 = clock();	
  for (int i=0; i<n_points; i++) {
    for (int j=0; j<n_points; j++) {
	  array<double,2> args = {grid[i], grid[j]};
      triang.insert(args.begin(), args.end(), fn(2, &args[0]));
    }
  }
  t2 = clock();
  printf("regular grid: %d insertions, %d clocks, %f sec\n", n_points*n_points, (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  // second, try adaptive point placement
  std::function<double(int, double*)> fn_obj(fn);
  Delaunay_incremental_interp_d adaptive_triang(2, fn_obj);
  // insert boundary points
  array<double,2> args;
  t1 = clock();
  args[0]=grid.front(); args[1]=grid.front();	adaptive_triang.insert(args.begin(), args.end(), fn(2, &args[0]));
  args[0]=grid.front(); args[1]=grid.back();	adaptive_triang.insert(args.begin(), args.end(), fn(2, &args[0]));
  args[0]=grid.back(); args[1]=grid.front();	adaptive_triang.insert(args.begin(), args.end(), fn(2, &args[0]));
  args[0]=grid.back(); args[1]=grid.back();		adaptive_triang.insert(args.begin(), args.end(), fn(2, &args[0]));
  for (int i=0; i<n_points*n_points-4; i++) {
    adaptive_triang.insert_largest_error_point();
  }	
  t2 = clock();
  printf("adaptive grid: %d insertions, %d clocks, %f sec\n", n_points*n_points, (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  // compare interpolated value vs. actual function
  std::vector<double> true_f_vals, interp_grid = linspace(0.01 + 0.005, 5. - + 0.005, 2*n_points);  
  for (int i=0; i<interp_grid.size(); i++) {
    for (int j=0; j<interp_grid.size(); j++) {
	  array<double,2> args = {interp_grid[i], interp_grid[j]};
	  true_f_vals.push_back(fn(2, &args[0]));
	}
  }

  // get the interpolated values  
  std::vector<double> regular_triang_vals, adaptive_triang_vals;
  t1 = clock();	  
  for (int i=0; i<interp_grid.size(); i++) {
    for (int j=0; j<interp_grid.size(); j++) {
	  array<double,2> args = {interp_grid[i], interp_grid[j]};
	  regular_triang_vals.push_back(triang.interp(args.begin(), args.end()));
	}
  }
  t2 = clock();
  printf("regular grid: %d interpolations, %d clocks, %f sec\n", regular_triang_vals.size(), (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);
  t1 = clock();	  
  for (int i=0; i<interp_grid.size(); i++) {
    for (int j=0; j<interp_grid.size(); j++) {
	  array<double,2> args = {interp_grid[i], interp_grid[j]};
	  adaptive_triang_vals.push_back(adaptive_triang.interp(args.begin(), args.end()));
	}
  }
  t2 = clock();
  printf("adaptive grid: %d interpolations, %d clocks, %f sec\n", adaptive_triang_vals.size(), (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);
  
  // compute sum of squared errors
  double sse1=0.0, sse2=0.0, diff;
  for (int i=0; i<true_f_vals.size(); i++) {
    diff = true_f_vals[i] - regular_triang_vals[i];
    sse1 += diff*diff;
	diff = true_f_vals[i] - adaptive_triang_vals[i];
	sse2 += diff*diff;
  }
  printf("regular grid: sum of squared errors: %f\n", sse1);
  printf("adaptive grid: sum of squared errors: %f\n", sse2);
  
  return 0;
}

