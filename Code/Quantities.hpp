#ifndef QUANTITIES
#define QUANTITIES


double partitionfunction(const Eigen::Matrix<double, -1, 1, 0, -1, 1> eigvals, double beta)
{
  double Z = 0;
  double minval = eigvals[0];
  for (int i = 0; i < eigvals.size(); i++) Z += exp(-beta*(eigvals[i]-minval));
  return Z;
}

double SzSz(const Eigen::Matrix<double, -1, 1, 0, -1, 1> eigvals, Matrix<double,Dynamic,Dynamic> eigvecs, double beta, double t)
{
  return 1;
}


#endif
