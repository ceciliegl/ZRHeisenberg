#ifndef INDEXSTATE
#define INDEXSTATE

using namespace std;

class indexstate
{
public:

  vector<short unsigned int> state_to_index;
  vector<unsigned int> index_to_state;

  indexstate();
  indexstate(int MaxIndex0, int MaxState0);
  void init_index_to_state(int MaxIndex0);
  void init_state_to_index(int MaxState0);
};

indexstate::indexstate(){}

indexstate::indexstate(int MaxIndex0, int MaxState0)
{
  index_to_state = vector<unsigned int>(MaxIndex0);
  state_to_index = vector<short unsigned int>(MaxState0);
}

void indexstate::init_index_to_state(int MaxIndex0)
{
  index_to_state = vector<unsigned int>(MaxIndex0);
}

void indexstate::init_state_to_index(int MaxState0)
{
  state_to_index = vector<short unsigned int>(MaxState0);
}

#endif
