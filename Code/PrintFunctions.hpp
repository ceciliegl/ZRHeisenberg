#ifndef PRINTFUNCTIONS
#define PRINTFUNCTIONS

#include <fstream>

class PrintFunctions
{
public:
  string dir;
  Solver mysystem;

  PrintFunctions(string dir0, Solver mysystem0);

  void WriteEigvals();
  void resetdatafiles();
};

PrintFunctions::PrintFunctions(string dir0, Solver mysystem0)
{
  dir = dir0;
  mysystem = mysystem0;
}


void PrintFunctions::WriteEigvals()
{
  ofstream Outfile(dir + "eigvals.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  cout << "Hei!" << mysystem.nu << endl;

  Outfile.precision(17);
  Outfile << mysystem.nu << "      " << mysystem.eigenvals.transpose() << "\n";
}


void PrintFunctions::resetdatafiles()
{
  ofstream EigFile(dir + "eigvals.txt");
  if (!EigFile.is_open())
     cout<<"Could not open file" << endl;
}

#endif
