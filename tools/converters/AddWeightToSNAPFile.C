
#include "../../core/common/parallel.h"
#include "../../core/common/parseCommandLine.h"
#include "../common/graphIO.h"
#include <stdio.h>
#include <iostream>
#include<string>
#include<fstream>
#include<cstdlib>
using namespace std;
int parallel_main(int argc, char *argv[]) {
  commandLine P(argc, argv, "<input SNAP file> <output Ligra file>");
  char *iFile = P.getArgument(1);
  char *oFile = P.getArgument(0);
  cout << "Rewrite File ... \n";
  ifstream iFileStream(iFile);
  ofstream oFileStream(oFile);
  string tmp;
  long  lines = 0;
  while(++lines){
    getline(iFileStream,tmp);
    if (iFileStream.eof()) break;
    oFileStream << tmp << "\t1" <<"\n";
  }
  cout << "Rewrite " <<  lines <<" lines. \n";
}