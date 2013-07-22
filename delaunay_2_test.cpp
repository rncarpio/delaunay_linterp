
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
      args[2] = grid2[j];
      f_val = cd_fn(args.begin());
      triang.insert(args.begin(), f_val);
    }
  }
  t2 = clock();
  printf("%d insertions, %d clocks, %f sec\n", n_points*n_points, (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  return 0;
}

