
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
#include "delaunay_2_interp.h"

// return an evenly spaced 1-d grid of doubles.
std::vector<double> linspace(double first, double last, int len) {
  std::vector<double> result(len);
  double step = (last-first) / (len - 1);
  for (int i=0; i<len; i++) { result[i] = first + i*step; }
  return result;
}

// the function to interpolate.
double fn (double x1, double x2) { return sin(x1 + x2); }

int main(int argc, char **argv) {
  int n_points = 10;
  clock_t t1, t2;  
  array<double,2> args;
  double f_val;

  // first, try a rectangular grid
  Delaunay_incremental_interp_2 triang;
  std::vector<double> grid = linspace(0., 5., n_points);
  t1 = clock();	
  for (int i=0; i<n_points; i++) {
    for (int j=0; j<n_points; j++) {
	  args = {grid[i], grid[j]};
      triang.insert(args.begin(), args.end(), fn(args[0], args[1]));
    }
  }
  t2 = clock();
  printf("%d insertions, %d clocks, %f sec\n", n_points*n_points, (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  // second, try adaptive point placement
  std::function<double(double, double)> fn_obj(fn);
  Delaunay_incremental_interp_2 adaptive_triang(fn_obj);
  // insert boundary points
  args = {grid.front(), grid.front()};	adaptive_triang.insert(args.begin(), args.end(), fn(args[0], args[1]));
  args = {grid.front(), grid.back()};	adaptive_triang.insert(args.begin(), args.end(), fn(args[0], args[1]));
  args = {grid.back(), grid.front()};	adaptive_triang.insert(args.begin(), args.end(), fn(args[0], args[1]));
  args = {grid.back(), grid.back()};	adaptive_triang.insert(args.begin(), args.end(), fn(args[0], args[1]));
  for (int i=0; i<n_points*n_points-4; i++) {
    adaptive_triang.insert_largest_error_point();
  }	
  
  return 0;
}

