#ifndef BASICFUNCTIONS
#define BASICFUNCTIONS

complex<double> exponential(double phi)
{
  //returns exp(i\phi)
  complex<double> a(cos(phi), sin(phi));
  return a;
}

string dectobin(int dec)
{
  string bin = "";
  while(dec >= 1)
  {
    bin += to_string(dec%2);
    dec /= 2;
  }

  int n = bin.length();
  for (int i = 0; i < n/2; i++)
        swap(bin[i], bin[n-i-1]);
  return bin;
}

unsigned int atotheb(int a, int b)
{
  unsigned int ans = 1;
  for (int i = 0; i < b; i++){ans*=a;} //Works
  return ans;
}

unsigned int factorial(int num)
{
  if (num > 0) {return num*factorial(num-1);}
  else if (num == 0) {return 1;}
  else {return 0;}
}

unsigned int binomial(const int n, const int k)
{
  int ans;
  if (k > 0)
  {
    ans = (n-k+1);
    for (int i = 1; i < k; i++)
    {
      ans *= (n-k+1 + i);
      ans /= (1 + i);
    }
  }
  else if (k == 0) ans = 1;
  else ans = 0;

  return ans;
}

unsigned int numstates(int TOT, int var1, int var2, int var3)
{
  //Find the number of possibilities to distribute var1+var2+var3=TOT.

  return factorial(TOT)/(factorial(var1)*factorial(var2)*factorial(var3)); //May be scary for large numbers, but I think we will work with small enough numbers?
}

double findminimum(vector<double> A)
{
  double ans = numeric_limits<double>::max();
  for(int i = 0; i < A.size(); i++)
  {
    if(A[i] < ans) ans = A[i];
  }

  return ans;
}

#endif
