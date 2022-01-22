#ifndef SOLVER
#define SOLVER


// Use Eigen to diagonalize? Or LAPACK? Try Eigen.
#include <Eigen/Dense>
using namespace Eigen;

#include "Quantities.hpp"

#include <fstream>

class Solver
{
public:
  string dir;
  int L, TWOL, Nh, Ns, nu;
  int maxIndexValue;
  int maxStateValue;
  unsigned int twomax; //Are these used in  Solver?
  vector<int> TWOLpow; //Are these used in  Solver?

  double tl, tr, Jzl, Jzr, Jpml, Jpmr; //t, J_z and J_{\pm} for legs and rungs.

  indexstate converttable;

  Matrix<double,Dynamic,Dynamic> H;

  Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvals;
  Matrix<double,Dynamic,Dynamic> eigenvecs;

  bool EIGVEC, CORR;

  complex<double> zero;
  int Nb, Nt;


  Solver();
  Solver(string dir0, int L0, int Nh0, double tl0, double tr0, double Jzl0, double Jzr0, double Jpml0, double Jpmr0, bool EIGVEC0, bool CORR0);

  void solve();
  void makebasis();
  void fillH();
  void diagonalise();

  unsigned int statevec_to_statenum(vector<short int> statevec);
  vector<short int> statenum_to_statevec(unsigned int statenum);

  vector<double> Szmat(int siteind);
  //complex<double> SzCorr(int i, int j, double beta, double t);
  vector<vector<vector<complex<double>>>> SzCorr(vector<double> beta, vector<double> time);
  Matrix<double, Dynamic, Dynamic> makeSminus(indexstate converttablep);
  //complex<double> SpmCorr(int i, int j, Matrix<double, Dynamic, Dynamic> Sminusmat, Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, double beta, double t);
  vector<vector<vector<complex<double>>>> SpmCorr(Matrix<double, Dynamic, Dynamic> Sminusmat, Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time);

  double Sx2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);
  double Sy2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);
  double Sz2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);
  double S(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);

  void WriteEigvals();
  void WriteSzStot();
  void WritePartition(vector<double> beta, vector<double> partition);
  void WriteCorr(vector<double> beta, vector<double> time, vector<vector<vector<complex<double>>>> z, vector<vector<vector<complex<double>>>> pm);
  void resetdatafiles();

};

Solver::Solver(){}

Solver::Solver(string dir0, int L0, int Nh0, double tl0, double tr0, double Jzl0, double Jzr0, double Jpml0, double Jpmr0, bool EIGVEC0, bool CORR0)
{
  zero.real(0); zero.imag(0);

  Nb = 0; Nt =0;

  dir = dir0;

  L = L0;
  TWOL = 2*L;
  Nh = Nh0;
  nu = 0; //Default

  tl = tl0;
  tr = tr0;
  Jzl = Jzl0;
  Jzr = Jzr0;
  Jpml = Jpml0;
  Jpmr = Jpmr0;

  EIGVEC = EIGVEC0;
  CORR = CORR0;

  Ns = TWOL - Nh;
  twomax = 1<<Ns;
  TWOLpow = vector<int>(Nh);
  for (int i = 0; i < Nh; i++)
  {
    TWOLpow[i] = pow(TWOL, i);
  }
}


void Solver::solve()
{
  vector<double> mineigvals(Ns+1);

  Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp;
  Matrix<double,Dynamic,Dynamic> eigenvecsp;
  indexstate converttablep;

  vector<double> beta = {0, 0.5, 1, 2, 10, 100};
  vector<double> time = {0};

  Nb = beta.size();
  Nt = time.size();

  double start, stop;

  vector<vector<double>> partfunc(Ns+1, vector<double>(Nb, 0.0));       //Partition function for each magnetisation sector, scaled by e^(-beta mineigvals)
  vector<vector<vector<vector<complex<double>>>>> corrfunczz(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfuncpm(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.

  //Do nu=0 first:
  nu = 0;
  makebasis();
  fillH();
  diagonalise();
  mineigvals[0] = eigenvals[0];
  for(int b = 0; b < Nb; b++) partfunc[0][b] = partitionfunction(eigenvals, beta[b]);
  WriteEigvals();
  WriteSzStot();
  if(CORR)
  {
    cout << "IN CORR FOR nu=0" << endl;
    eigenvalsp = eigenvals;
    eigenvecsp = eigenvecs;
    converttablep = converttable;

    cout << "BEFORE SZCORR FOR nu=0" << endl;

    for(int j = 0; j < TWOL; j++) corrfunczz[0] = SzCorr(beta, time); //Only contribution from Sz?
    for(int j = 0; j < TWOL; j++) for(int b = 0; b < Nb; b++) for(int t = 0; t < Nt; t++) corrfuncpm[0][j][b][t] = zero;
  }

  Matrix<double, Dynamic, Dynamic> Sminus;

  for (int mynu = 1; mynu < Ns+1; mynu++)
  {
    nu = mynu;
    makebasis();
    fillH();
    diagonalise();
    mineigvals[mynu] = eigenvals[0];

    for(int b = 0; b < Nb; b++) partfunc[mynu][b] = partitionfunction(eigenvals, beta[b]);

    WriteEigvals();
    WriteSzStot();

    //Compute correlations here!

    if(CORR)
    {
      Sminus = makeSminus(converttablep);

      /*for(int j = 0; j < TWOL; j++)
        for(int b = 0; b < Nb; b++)
          for(int t = 0; t < Nt; t++)
          {
            start = clock(); //start clock
            corrfunczz[mynu][j][b][t] = SzCorr(0, j, beta[b], time[t]);
            stop = clock();

            cout << "Sz time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;

            start = clock(); //start clock
            corrfuncpm[mynu][j][b][t] = SpmCorr(0, j, Sminus, eigenvalsp, eigenvecsp, converttablep, beta[b], time[t]);
            stop = clock();

            cout << "Spm time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;
          }*/

      cout << "BEFORE SZCORR FOR nu=" << nu << endl;

      start = clock(); //start clock
      corrfunczz[mynu] = SzCorr(beta, time);
      stop = clock();

      cout << "Sz time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;

      start = clock(); //start clock
      corrfuncpm[mynu] = SpmCorr(Sminus, eigenvalsp, eigenvecsp, converttablep, beta, time);
      stop = clock();

      cout << "Spm time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;

      eigenvalsp = eigenvals;
      eigenvecsp = eigenvecs;
      converttablep = converttable;
    }
  }

  double GS = findminimum(mineigvals);
  cout << GS << endl;
  vector<double> partitionfunction(Nb, 0.0);
  vector<vector<vector<complex<double>>>> corrz(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));
  vector<vector<vector<complex<double>>>> corrpm(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  cout << endl;
  cout << endl;

  for(int i = 0; i < Ns+1; i++)
    for(int b = 0; b < Nb; b++)
    {
      //cout << partfunc[i]*exp(-beta*(mineigvals[i]-GS)) << endl;
      partitionfunction[b] += partfunc[i][b]*exp(-beta[b]*(mineigvals[i]-GS));//Sould double check this for a small system?
      for(int t = 0; t < Nt; t++)
        for(int j = 0; j < TWOL; j++)
        {
          //cout << "Here: " << corrfuncpm[i][j] << "   " << exp(-beta*(mineigvals[i]-GS)) << endl;
          corrz[j][b][t] += corrfunczz[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
          corrpm[j][b][t] += corrfuncpm[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
        }
    }

  for(int j = 0; j < TWOL; j++)
    for(int b = 0; b < Nb; b++)
      for(int t = 0; t < Nt; t++)
      {
        //cout << "Here: " << corrfuncpm[i][j] << "   " << exp(-beta*(mineigvals[i]-GS)) << endl;
        corrz[j][b][t] /= partitionfunction[b];
        corrpm[j][b][t] /= partitionfunction[b];
      }

  cout << "Feridg" << endl;

  if(CORR) WriteCorr(beta, time, corrz, corrpm);

  WritePartition(beta, partitionfunction);

  cout << "Ferdig" << endl;
}


void Solver::makebasis() //nu is number of up spins
{
  maxIndexValue = binomial(TWOL, Nh)*binomial(Ns, nu); //This is the number of states with Nh holes and nu up-spins.
  converttable.init_index_to_state(maxIndexValue);

  //Construct all states with nu up spins.

  //How do I want to represent a state? Possibly as a vector with 2L numbers -1, 0 or 1? Or as a vector of length Nh and a binary number with the rest?
  //Do it as a vector of length Nh and a binary number with the rest and translate it into spins on sites when computing quantities.
  //For 0 and 1 hole I can do s + 2^(2L-Nh)*holeposition. But what to do for two or more holes? Think about this for later.

  //First find all the possible binary numbers with nu ones. Loop over all numbers from 0 to 2^{2*L-Nh}

  vector<unsigned int> s = {}; //List of all spin configurations with nu up spins.

  //int nconf = 0;
  for (int i = 0; i < twomax; i++)
  {
    int nunow = 0;
    for (int j = 0; j < Ns; j++)
    {
      if (i >> j & 1) nunow += 1;
    }
    if (nunow == nu)
    {
      //cout << dectobin(i) << endl;
      //nconf += 1;
      s.push_back(i);
    }
  }

  int NsStates = s.size(); //Seems to be correct.

  short unsigned int index;
  unsigned int state;
  int maxStateValue = 0;

  //Then combine these with all possible hole positions to create a number representing a state. And make lists converting between state numbers and matrix indices.
  //Make an algorithm for general Nh:

  vector<int> holeind(Nh+1);
  vector<int> holemax(Nh+1);

  for (int i = 0; i < Nh; i++)
  {
    holeind[i] = Nh-1-i;
    holemax[i] = TWOL-i;
  }
  holeind[Nh] = 0;
  holemax[Nh] = 2;

  int p = 0;

  index = 0;
  while (holeind[Nh] == 0)
  {
    /*for (int kk = 0; kk < Nh; kk++)
    {
      cout << holeind[kk] << "   ";
    }
    cout << endl;*/
    for (int i = 0; i < NsStates; i++)
    {
      state = s[i];
      for (int j = 0; j < Nh; j++)
      {
        state += twomax*(TWOLpow[j]*holeind[Nh-1-j]); //NOT SURE WHETHER I AM COMBINING THE CORRECT TWOLpow AND ind?!
        //cout << TWOLpow[j] << "   " << holeind[Nh-1-j] << endl;
      }
      if (state > maxStateValue) maxStateValue = state;
      converttable.index_to_state[index] = state;
      //cout << state << endl;
      index++;
    }

    holeind[0]++;
    while (holeind[p] == holemax[p])
    {
      ++p;
      //cout << "p = " << p << endl;
      holeind[p]++;
      for (int q = p-1; q>=0; q--) holeind[q]=holeind[q+1]+1; //Feilen ligger her


      if(holeind[Nh] == 1){break;}
      if(holeind[p]!=holemax[p]){p=0;}
    }
  }

  converttable.init_state_to_index(maxStateValue+1);
  for (int i = 0; i < maxStateValue; i++) converttable.state_to_index[i] = maxIndexValue;
  for (int i = 0; i < maxIndexValue; i++)
  {
    converttable.state_to_index[converttable.index_to_state[i]] = i;
  }

  /*
  for (int i = 0; i < converttable.state_to_index.size(); i++)
  {
    cout << converttable.state_to_index[i] << endl;
  }

  for (int i = 0; i < converttable.index_to_state.size(); i++)
  {
    cout << converttable.index_to_state[i] << endl;
  }
  */

  //I think these work??? Look over it again tomorrow.


  //What to do about all the zeros which are not supposed to be zeros? I have set them to the maxIndexValue, i.e. out of range.
}

unsigned int Solver::statevec_to_statenum(vector<short int> statevec)
{
  int statenum = 0;

  int spin = 0;
  int holeexp = 0;

  for (int i = 0; i < statevec.size(); i++)
  {
    if (statevec[i] == 0)
    {
      statenum += i*twomax*pow(TWOL, holeexp);
      holeexp++;
      continue;
    }
    spin = spin << 1;
    if (statevec[i] == 1)
    {
      spin += 1;
    }
  }
  statenum += spin;
  return statenum;
}

vector<short int> Solver::statenum_to_statevec(unsigned int statenum)
{
  vector<short int> statevec(TWOL, -2);

  unsigned int spin = statenum%twomax;

  statenum /= twomax;
  //cout << statenum << endl;

  for (int i = 0; i < Nh; i++)
  {
    statevec[statenum/TWOLpow[Nh-1-i]] = 0;
    statenum %= TWOLpow[Nh-1-i];
  }

  for (int i = TWOL-1; i >= 0; i--)
  {
    if (statevec[i] != 0)
    {
      statevec[i] += 1 + 2*(spin%2);
      spin = spin >> 1;
    }
  }

  return statevec;
}

void Solver::fillH()
{
  //Compute all elements of H for a given nu?

  //maxIndexValue = binomial(2*L, Nh)*binomial(Ns, nu);

  H.resize(maxIndexValue,maxIndexValue);
  H.setZero();

  unsigned int innumber;
  vector<short int> instate;
  vector<short int> holevec(Nh);
  int holepos;
  int neighpos;
  int neigh2;
  int additionalsign;

  int j;

  for (int i = 0; i < maxIndexValue; i++)
  {

    innumber = converttable.index_to_state[i];
    instate = statenum_to_statevec(innumber);

    innumber /= twomax;

    for (int h = 0; h < Nh; h++)
    {
      holevec[Nh-1-h] = innumber/TWOLpow[Nh-1-h];
      innumber %= TWOLpow[Nh-1-h];
    }

    //Diagonal terms are easy: SziSzj for Jleg and Jrung add a Delta to tune z-component?:

    for (int site = 0; site < TWOL; site++)
    {
      H(i, i) += 0.25*2*Jzr*instate[site]*instate[(site+1)%TWOL]; //Should I include the 0.25 from s=1/2?
      H(i, i) += 0.25*Jzl*instate[site]*instate[(site+2)%TWOL]; //Should I include the 0.25 from s=1/2?
    }

    //Off-diagonal terms:

    //Heisenberg part:
    //Nearest neighbours may flip if they are one up and one down.

    for (int site = 0; site < TWOL; site++)
    {
      if (instate[site]*instate[(site+1)%TWOL] == -1)
      {
        instate[site] *= (-1);
        instate[(site+1)%TWOL] *= (-1);

        j = converttable.state_to_index[statevec_to_statenum(instate)];

        H(i,j) += Jpmr;

        instate[site] *= (-1);
        instate[(site+1)%TWOL] *= (-1);
      }
      if (instate[site]*instate[(site+2)%TWOL] == -1)
      {
        instate[site] *= (-1);
        instate[(site+2)%TWOL] *= (-1);

        j = converttable.state_to_index[statevec_to_statenum(instate)];

        H(i,j) += 0.5*Jpml;

        instate[site] *= (-1);
        instate[(site+2)%TWOL] *= (-1);
      }
    }

    //Hopping part:
    for (int neigh = -1; neigh < 2; neigh += 2)
    {
      neigh2 = 2*neigh;
      for (int hole = 0; hole < Nh; hole++)
      {
        holepos = holevec[hole];

        //Hop along rungs
        neighpos = (holepos+neigh+TWOL)%TWOL;

        if ((neigh == -1 && holepos == 0) || (neigh == +1 && holepos == TWOL-1)) {additionalsign = +1 + 2*(-1)*((TWOL-2-(Nh-1))%2);}
        else {additionalsign = +1;}

        instate[holepos] = instate[neighpos];
        instate[neighpos] = 0;
        j = converttable.state_to_index[statevec_to_statenum(instate)];

        H(i,j) += -2*tr*additionalsign*abs(instate[holepos]);

        instate[neighpos] = instate[holepos];
        instate[holepos] = 0;


        //Hop along legs.
        neighpos = (holepos+neigh2+TWOL)%TWOL;

        if ((neigh2 == -2 && holepos == 0) || (neigh2 == +2 && holepos == TWOL-2)){additionalsign = (+1) + 2*(-1)*((TWOL-3-(Nh-1+(int(abs(instate[TWOL-1]))-1)))%2);}
        else if ((neigh2 == -2 && holepos == 1) || (neigh2 == +2 && holepos == TWOL-1)) {additionalsign = (+1) + 2*(-1)*((TWOL-3-(Nh-1+(int(abs(instate[0]))-1)))%2);}
        else if (instate[(holepos+neigh+TWOL)%TWOL] == 0) {additionalsign = +1;}
        else {additionalsign = -1;}

        instate[holepos] = instate[neighpos];
        instate[neighpos] = 0;
        j = converttable.state_to_index[statevec_to_statenum(instate)];
        H(i,j) += -tl*additionalsign*abs(instate[holepos]);
        //cout << i << "   " << j << "   " << -tl*additionalsign*abs(instate[holepos]) << endl;

        instate[neighpos] = instate[holepos];
        instate[holepos] = 0;
      }
    }
  }


  return;
}

void Solver::diagonalise()
{
  SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic>> es;
  if (EIGVEC)
  {
    es = SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic>>(H,ComputeEigenvectors);
    eigenvecs = es.eigenvectors();
  }
  else {es = SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic>>(H,EigenvaluesOnly);}

  eigenvals = es.eigenvalues();

  return;
}

vector<double> Solver::Szmat(int siteind)
{
  //Sovle it as a sparse matrix maybe? It is diagonal?
  vector<double> ans(maxIndexValue);
  double statenumber;
  vector<short int> statevec;

  for (int i = 0; i < maxIndexValue; i++)
  {
    statenumber = converttable.index_to_state[i];
    statevec = statenum_to_statevec(statenumber);
    ans[i] = 0.5*statevec[siteind];
  }

  return ans;
}

vector<vector<vector<complex<double>>>> Solver::SzCorr(vector<double> beta, vector<double> time)
{
  vector<vector<vector<complex<double>>>> SzSz(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  double nSz0m = 0;
  double mSzjn = 0;

  double minval = eigenvals[0];

  //Compute SzSz correlations between site i and all other sites.
  vector<double> Sz0 = Szmat(0);
  vector<vector<double>> Szj(TWOL, vector<double>(maxIndexValue, 0.0));
  for(int j = 0; j < TWOL; j++) Szj[j] = Szmat(j);

  Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvecn, eigenvecm;
  double eigenvecnmalpha;
  vector<double> Szjnow(maxIndexValue);

  for (int n = 0; n < maxIndexValue; n++)
  {
    eigenvecn = eigenvecs.col(n);
    for (int m = 0; m < maxIndexValue; m++)
    {
      eigenvecm = eigenvecs.col(m);
      for(int j = 0; j < TWOL; j++)
      {
        Szjnow = Szj[j];
        for (int alpha = 0; alpha < maxIndexValue; alpha ++)
        {
          //This is probably an extremely stupid way to do this?
          eigenvecnmalpha = eigenvecn[alpha]*eigenvecm[alpha];
          nSz0m += eigenvecnmalpha*Sz0[alpha];
          mSzjn += eigenvecnmalpha*Szjnow[alpha];
        }
        //cout << eigenvals[n]-minval << endl;

        for(int b = 0; b < Nb; b++)
          for(int t = 0; t < Nt; t++)
            SzSz[j][b][t] += nSz0m*mSzjn*exponential(-(eigenvals[n]-eigenvals[m])*time[t])*exp(-beta[b]*(eigenvals[n]-minval));
        nSz0m = 0; mSzjn = 0;
      }
    }
  }
  return SzSz;
}

Matrix<double, Dynamic, Dynamic> Solver::makeSminus(indexstate converttablep)
{
  Matrix<double, Dynamic, Dynamic> Sminusmat(maxIndexValue, TWOL); //Matrix with one row for each basis state in the current mag sector and one column for each site i.

  vector<short int> statevec;
  unsigned int statenum;

  for(int i = 0; i < TWOL; i++)
  {
    for(int alpha = 0; alpha < maxIndexValue; alpha++)
    {
      statevec = statenum_to_statevec(converttable.index_to_state[alpha]);
      if(statevec[i] == +1)
      {
        statevec[i] = -1;
        statenum = statevec_to_statenum(statevec);
        Sminusmat(alpha, i) = converttablep.state_to_index[statenum];
      }
      else Sminusmat(alpha, i) = -1;
    }
  }

  return Sminusmat;
}

vector<vector<vector<complex<double>>>> Solver::SpmCorr(Matrix<double, Dynamic, Dynamic> Sminusmat, Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time)
{
  //Compute the correlations in the xy plane: SxiSxj + SyiSyj.

  //Only for one magnetisation sector. All mag sectors need to be added and divided by the partition function.
  //Will scale each term by the minimal energy in the corresponding sector.

  Eigen::Matrix<double, -1, 1, 0, -1, 1> Sminus0 = Sminusmat.col(0);
  //vector<Eigen::Matrix<double, -1, 1, 0, -1, 1>> Sminusj(TWOL);

  vector<vector<complex<double>>> W(Nb, vector<complex<double>>(Nt, zero));
  vector<vector<vector<double>>>partsum(TWOL, vector<vector<double>>(maxIndexValue, vector<double>(converttablep.MaxIndex, 0.0)));
  vector<vector<vector<complex<double>>>>totsum(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  double sum0, sumj;

  double minval = eigenvals[0];

  /*for(int n = 0; n < maxIndexValue; n++)
  {
    for(int m = 0; m < converttablep.MaxIndex; m++)
    {
      //cout << n << "   " << m << endl;

      W = 0.5*(exp(-beta*(eigenvals[n]-minval))*exponential(-(eigenvals[n]-eigenvalsp[m])*t)+exp(-beta*(eigenvalsp[m]-minval))*exponential(-(eigenvalsp[m]-eigenvals[n])*t)); //Exponential legger til en i i eksponenten.

      sumi = 0;
      sumj = 0;

      for(int alpha = 0; alpha < maxIndexValue; alpha++)
      {
        if(Sminusj[alpha] >= 0) sumj += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[Sminusj[alpha]];
        if(Sminusi[alpha] >= 0) sumi += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[Sminusi[alpha]];
      }
      totsum += W*sumi*sumj;
    }
  }*/

  for(int n = 0; n < maxIndexValue; n++)
  {
    for(int m = 0; m < converttablep.MaxIndex; m++)
    {
      sum0 = 0;
      for(int alpha = 0; alpha < maxIndexValue; alpha++)
      {
        int beta = Sminus0[alpha];
        if(beta >= 0) sum0 += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[beta];
      }
      for(int j = 0; j < TWOL; j++)
      {
        sumj = 0;

        for(int alpha = 0; alpha < maxIndexValue; alpha++)
        {
          int beta = Sminusmat.col(j)[alpha];
          if(beta >= 0) sumj += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[beta];
        }
        partsum[j][n][m] += sum0*sumj;
      }
    }
  }

  for(int n = 0; n < maxIndexValue; n++)
    for(int m = 0; m < converttablep.MaxIndex; m++)
      for(int b = 0; b < Nb; b++)
        for(int t = 0; t < Nt; t++)
        {
          W[b][t] = 0.5*(exp(-beta[b]*(eigenvals[n]-minval))*exponential(-(eigenvals[n]-eigenvalsp[m])*time[t])+exp(-beta[b]*(eigenvalsp[m]-minval))*exponential(-(eigenvalsp[m]-eigenvals[n])*time[t])); //Exponential legger til en i i eksponenten.
          for(int j = 0; j < TWOL; j++) totsum[j][b][t] += W[b][t]*partsum[j][n][m];
        }

  return totsum;
}

double Solver::Sx2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{

  //double Solver::Sx2(vector<double> statecoeffs, vector<double> statenumbers)

  //Compute Sx2 for a state with statecoeffs and corresponding statenumbers. I.e. state coeffs gives the coefficients for the different basis states found in statenumbers.

  //I guess I already have the statenumbers if i calculate Sx2 when I am computing the eigenvalues?
  double ans = 0;
  double maxstatenum = (Nh == 0) ? twomax : twomax*(1 + TWOL*TWOLpow[Nh-1]*Nh); //This is an overestimation for Nh>0!
  vector<double> intermediatestates(maxstatenum, 0.0);

  vector<short int> statevec;
  vector<short int> outvec;
  unsigned int outnum;

  for(int i = 0; i < maxIndexValue; i++)
  {
    statevec = statenum_to_statevec(converttable.index_to_state[i]);

    for(int j = 0; j < TWOL; j++)
    {
      outvec = statevec;
      if(statevec[j] == +1 || statevec[j] == -1)
      {
        outvec[j] *= (-1);
        outnum = statevec_to_statenum(outvec);
        intermediatestates[outnum] += statecoeffs[i];
      }
    }
  }

  for(int i = 0; i < intermediatestates.size(); i++)
  {
    ans += intermediatestates[i]*intermediatestates[i];
  }

  return ans*0.25;
}

double Solver::Sy2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{
  //Compute Sy2 for a state with statecoeffs and corresponding statenumbers. I.e. state coeffs gives the coefficients for the different basis states found in statenumbers.

  //I guess I already have the statenumbers if i calculate Sy2 when I am computing the eigenvalues?

  double ans = 0;
  double maxstatenum = (Nh == 0) ? twomax : twomax*(1 + TWOL*TWOLpow[Nh-1]*Nh); //This is an overestimation for Nh>0!
  vector<double> intermediatestates(maxstatenum, 0.0);

  vector<short int> statevec;
  vector<short int> outvec;
  unsigned int outnum;

  for(int i = 0; i < maxIndexValue; i++)
  {
    statevec = statenum_to_statevec(converttable.index_to_state[i]);

    for(int j = 0; j < TWOL; j++) //For each site
    {
      outvec = statevec;
      if(statevec[j] == +1 || statevec[j] == -1)
      {
        outvec[j] *= (-1);
        outnum = statevec_to_statenum(outvec);
        intermediatestates[outnum] += -statevec[j]*statecoeffs[i];
      }
    }
  }

  for(int i = 0; i < intermediatestates.size(); i++)
  {
    ans += intermediatestates[i]*intermediatestates[i];
  }

  return ans*0.25;
}

double Solver::Sz2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{
  //Compute Sz2 for a state with statecoeffs and corresponding statenumbers. I.e. state coeffs gives the coefficients for the different basis states found in statenumbers.

  //I think this just is summing the squares of the statecoeffs and mulitply by 0.25?
  //How does that make sense. Aren't the statecoeffs normalised?! Yes, it does not make sense, you need to include the total Sz?
  // So for a given nu, Sz2 is just 0.25*(2*nu - Ns)**2

  return 0.25*(2*nu - Ns)*(2*nu - Ns);
}

double Solver::S(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{
  double ans = Sx2(statecoeffs) + Sy2(statecoeffs) + Sz2(statecoeffs);
  return 0.5*(-1+sqrt(1+4*ans));
}


void Solver::WriteEigvals()
{
  ofstream Outfile(dir + "eigvals.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  Outfile.precision(17);
  Outfile << nu << "      " << eigenvals.transpose() << "\n";
}


void Solver::WriteSzStot()
{
  ofstream Outfile(dir + "SzStot.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  Outfile.precision(17);
  for(int i = 0; i < maxIndexValue; i++)
  {
    Outfile << 0.5*(2*nu - Ns) << "   " << S(eigenvecs.col(i)) << endl;
  }
  Outfile << "\n";
}

void Solver::WriteCorr(vector<double> beta, vector<double> time, vector<vector<vector<complex<double>>>> z, vector<vector<vector<complex<double>>>> pm)
{
  cout << "Starter her" << endl;

  ofstream Outfile(dir + "Corr.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  cout << "Inni her" << endl;

  Outfile.precision(17);
  for(int b = 0; b < beta.size(); b++)
    for(int t = 0; t < time.size(); t++)
    {
      Outfile << beta[b] << "   " << time[t] << "   ";
      for(int i = 0; i < z.size(); i++) Outfile << z[i][b][t] << "   ";
      for(int i = 0; i < pm.size(); i++) Outfile << pm[i][b][t] << "   ";
      Outfile << endl;
    }
}

void Solver::WritePartition(vector<double> beta, vector<double> partition)
{
  ofstream Outfile(dir + "Partition.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  for(int b = 0; b < beta.size(); b++)
  {
    Outfile << beta[b] << "   " << partition[b] << endl;
  }
}


void Solver::resetdatafiles()
{
  ofstream EigFile(dir + "eigvals.txt");
  if (!EigFile.is_open())
     cout<<"Could not open file" << endl;

  ofstream SFile(dir + "SzStot.txt");
  if (!SFile.is_open())
    cout<<"Could not open file" << endl;

  ofstream PartitionFile(dir + "Partition.txt");
  if (!PartitionFile.is_open())
    cout<<"Could not open file" << endl;

  ofstream CorrFile(dir + "Corr.txt");
  if (!CorrFile.is_open())
    cout<<"Could not open file" << endl;
}

#endif
