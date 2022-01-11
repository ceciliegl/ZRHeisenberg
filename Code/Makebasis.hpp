#ifndef MAKEBASIS
#define MAKEBASIS

class Makebasis
{
public:
  int nu;
  int TWOL, Nh, Ns;
  vector<unsigned int> s;

  indexstate converttable;

  unsigned int twomax;
  vector<int> TWOLpow;

  int NsStates;
  short unsigned int index;
  unsigned int state;
  int maxIndexValue, maxStateValue;

  vector<int> holeind;
  //vector<int> holemin;
  vector<int> holemax;

  Makebasis(int nu0, int TWOL0, int Nh0);

  void find_s();
  void RecursiveHoleLoops(int start, int stop, int holecount);
};

Makebasis::Makebasis(int nu0, int TWOL0, int Nh0)
{
  nu = nu0;
  TWOL = TWOL0;
  Nh = Nh0;

  Ns = TWOL - Nh;
  twomax = pow(2, Ns);
  TWOLpow = vector<int>(Nh);
  for (int i = 0; i < Nh; i++)
  {
    TWOLpow[i] = pow(TWOL, i);
  }

  maxIndexValue = binomial(TWOL, Nh)*binomial(Ns, nu);
  converttable.init_index_to_state(maxIndexValue);

  find_s();

  NsStates = s.size(); //Seems to be correct.
  maxStateValue = 0;

  holeind = vector<int>(Nh+1);
  //holemin = vector<int>(Nh+1);
  holemax = vector<int>(Nh+1);

  //RecursiveHoleLoops(Nh, TWOL, 0);

  for (int i = 0; i < Nh+1; i++)
  {
    holeind[i] = Nh-1-i;
    //holemin[i] = Nh-1-i;
    holemax[i] = TWOL-i;
  }
  holeind[Nh] = 0;

  int p = 0;

  index = 0;
  while (holeind[Nh] == 0)
  {
    for (int kk = 0; kk < Nh; kk++)
    {
      cout << holeind[kk] << "   ";
    }
    cout << endl;
    for (int i = 0; i < NsStates; i++)
    {
      state = s[i];
      for (int j = 0; j < Nh; j++)
      {
        state += twomax*(TWOLpow[j]*holeind[Nh-1-j]); //NOT SURE WHETHER I AM COMBINING THE CORRECT TWOLpow AND ind?!
      }
      if (state > maxStateValue) maxStateValue = state;
      converttable.index_to_state[index] = state;
      index++;
    }

    holeind[0]++;
    while (holeind[p] == holemax[p])
    {
      holeind[++p]++;
      holeind[p-1] = holeind[p]+1;

      if(holeind[p]!=holemax[p]){p=0;}
    }
  }


  //cout << NsStates << endl;

  converttable.init_state_to_index(maxStateValue+1);
  for (int i = 0; i < maxStateValue; i++) converttable.state_to_index[i] = maxIndexValue;
  for (int i = 0; i < maxIndexValue; i++)
  {
    converttable.state_to_index[converttable.index_to_state[i]] = i;
  }
}

void Makebasis::find_s()
{
  s = {}; //List of all spin configurations with nu up spins.

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
}

void Makebasis::RecursiveHoleLoops(int start, int stop, int holecount=0)
{
  //First start should be Nh?
  //First stop should be TWOL?

  cout << holecount << endl;

  if (holecount == Nh) return;

  for (int kk = 0; kk < Nh; kk++)
  {
    cout << holeind[kk] << "   ";
  }
  cout << endl;

  cout << "HI!" << endl;
  for (int i = 0; i < NsStates; i++)
  {
    state = s[i];
    for (int j = 0; j < Nh; j++)
    {
      state += twomax*(TWOLpow[j]*holeind[Nh-1-j]); //NOT SURE WHETHER I AM COMBINING THE CORRECT TWOLpow AND ind?!
    }
    if (state > maxStateValue) maxStateValue = state;
    converttable.index_to_state[index] = state;
    index++;
  }

  if (holecount == 0)
  {
    for (holeind[holecount] = Nh-1; holeind[holecount] < TWOL; holeind[holecount]++)
    {
      return RecursiveHoleLoops(start-1, stop-1, holecount+1);
    }
  }
  else
  {
    for (holeind[holecount] = Nh-1-holecount; holeind[holecount] < holeind[holecount-1]; holeind[holecount]++)
    {
      return RecursiveHoleLoops(start-1, stop-1, holecount+1);
    }
  }
}

#endif
