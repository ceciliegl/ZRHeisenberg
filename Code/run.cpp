// Write a code that solves ZR+Heisenberg model using ED for triangular+honeycomb tube.

#include<vector>
#include<complex>
#include<iostream>
#include <iomanip>
#include<string>
#include <time.h>

using namespace std;

#include "ReadInputFiles.hpp"
#include "BasicFunctions.hpp"
#include "indexstate.hpp"
#include "Solver.hpp"
#include "Quantities.hpp"
#include "PrintFunctions.hpp"


// Symmetry: Conserved total Sz. So I can block diagonalise H into sectors of constant total Sz.

int main(int argc, char const *argv[])
{
  string project = argv[1];
  string run_number;

  double beta = 1;

  if(atoi(argv[2]) < 10)
  {
      run_number = std::string("00") + std::string(argv[2]);
  }
  else if(atoi(argv[2]) > 9 && atoi(argv[2]) < 100)
  {
      run_number = std::string("0") + std::string(argv[2]);
  }
  else if(atoi(argv[2]) > 99)
  {
      run_number = std::string(argv[2]);
  }

  string folder = std::string(getenv("HOME")) + std::string("/Documents/ZRHeisenberg/Data/") + std::string(project) + std::string("/Run") + std::string(run_number) + std::string("/");

  cout << folder << endl;

  ReadInputFiles params(folder + "parameters.txt");
  params.generate();

  Solver mysolver(folder, params.L, params.Nh, params.tl, params.tr, params.Jzl, params.Jzr, params.Jpml, params.Jpmr, params.EIGVECS, params.CORR);
  //PrintFunctions printer(folder, mysolver);

  if(params.RESETFILES){mysolver.resetdatafiles();}

  vector<double> Sz0Szj(2*params.L, 0.0);

  /*int maxnu = 2*params.L+1-params.Nh;

  vector<double> mineigvals(maxnu);
  vector<double> partfunc(maxnu);       //Partition function for each magnetisation sector, scaled by e^(-beta mineigvals)

  Eigen::Matrix<double, -1, 1, 0, -1, 1> eigvalsp;
  Matrix<double,Dynamic,Dynamic> eigvecsp;
  indexstate converttablep;*/

  //vector<Eigen::Matrix<double, -1, 1, 0, -1, 1>> allvals(maxnu);
  //vector<Matrix<double,Dynamic,Dynamic>> allvecs(maxnu);

  mysolver.solve();

  /*for (int mynu = 1; mynu < maxnu; mynu++)
  {
    //cout << "Making basis..." << endl;

    mysolver.nu = mynu;
    mysolver.solve();

    cout << fixed << setprecision(2) << setfill('0');

    for (int i = 0; i < mysolver.maxIndexValue; i++)
    {
      for (int j = 0; j < mysolver.maxIndexValue; j++)
      {
        cout << setw(5) << mysolver.H(i, j) << "   ";
      }
      cout << endl;
    }

    //cout << "Diagonalising..." << endl;

    //double start = clock(); //start clock

    //SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic>> es(solver.H,ComputeEigenvectors);

    mineigvals[mynu] = mysolver.eigenvals[0];

    //cout << "Etter diagonalisering" << endl;

    //allvals[mynu] = mysolver.eigenvals;
    //allvecs[mynu] = mysolver.eigenvecs;

    //double stop = clock(); //start clock

    //double time_used = (stop-start)/CLOCKS_PER_SEC; //print time used

    //cout << "L = " << L << "   " << "Matrix size = " << solver.converttable.index_to_state.size() << "x" << solver.converttable.index_to_state.size() << "   " << "Time used = " << time_used << "s" << endl;

    //cout << "Eigenvalues: ";
    //cout << es.eigenvalues().transpose() << endl;

    //partfunc += partitionfunction(mysolver.eigenvals, 100000); //Need to do something smart here to normalise w.r.t. the smallest eigenvalue?

    for (int j = 0; j < 2*params.L; j++)
    {
      Sz0Szj[j] += mysolver.SzCorr(0, j)/partfunc;
    }

    //partfunc[mynu] = partitionfunction(mysolver.eigenvals, beta);

    //mysolver.WriteEigvals();
    //mysolver.WriteSzStot();

    //Thinking a bout a way to only save some eigenvalues and eigenvectors at a time

    //eigvalspp = eigvalsp;
    //eigvecspp = eigvecsp;

    //Compute correlations here!

    eigvalsp = mysolver.eigenvals;
    eigvecsp = mysolver.eigenvecs;
    converttablep = mysolver.converttable;

  }

  double GS = findminimum(mineigvals);
  double partitionfunction = 0;
  for(int i = 0; i < maxnu; i++)
  {
    cout << partfunc[i]*exp(-beta*(mineigvals[i]-GS)) << endl;
    partitionfunction += partfunc[i]*exp(-beta*(mineigvals[i]-GS));//Sould double check this for a small system?
  }*/

  //cout << partitionfunction;

  //Having all the eigenvalues and eigenvectors, I can compute anything.
  //Compute Sx^2:



  //Compute Sy^2:

  //Copute Sz^2:



  /*for(int i = 0; i < solver.converttable.index_to_state.size(); i++)
  {
    cout << solver.converttable.index_to_state[i] << endl;
  }

  cout << fixed << setprecision(2) << setfill('0');

  for (int i = 0; i < solver.maxIndexValue; i++)
  {
    for (int j = 0; j < solver.maxIndexValue; j++)
    {
      cout << setw(5) << solver.H(i, j) << "   ";
    }
    cout << endl;
  }

  double partfuncsector = partitionfunction(es.eigenvalues(), 1);
  cout << partfuncsector << endl;

  cout << "Eigenvalues: ";
  cout << es.eigenvalues().transpose() << endl;
  cout << endl;
  cout << "Eigenvectors: " << endl;
  cout << es.eigenvectors() << endl;*/


  /*vector<short int> a = {1, -1, 0, 1, 0, 1, 1, -1, 1, 1, 1, 1, -1, -1};
  cout << solver.statevec_to_statenum(a) << endl;
  vector<short int> b = solver.statenum_to_statevec(solver.statevec_to_statenum(a));
  for (int i = 0; i < 2*L; i++) cout << b[i] << "   ";*/


  return 0;
}
