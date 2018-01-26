#ifndef FILE_OPERATION   /* Include guard */
#define FILE_OPERATION


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <armadillo>

using namespace std;
using namespace arma;

//--------------------------------------------------------------------------

mat getPos(char* myFile){

  ifstream myfilein(myFile);
  string line;
  double posX,posY,posZ;

  int j=0;

  for (int i = 0; getline(myfilein, line) ;++i) ++j;

  mat A(j - 2 ,3);

  ifstream myfilein2(myFile);

  j=0;

    while (getline(myfilein2, line))
    {
        istringstream ss(line);

        string name;
        ss >> name >> posX>> posY>> posZ;
        if (j >1) {
        A(j - 2,0) = posX;
        A(j - 2,1) = posY;
        A(j - 2,2) = posZ;
     }

        j++; 
    }
myfilein.close();
return A;

  
}
//--------------------------------------------------------------------------
string* getType(char* myFile){

  ifstream myfilein(myFile);
  string line;
  double posX,posY,posZ;

  int j=0;

  for (int i = 0; getline(myfilein, line) ;++i) ++j;

  string* mytype = new string[j - 2];

  ifstream myfilein2(myFile);

  j=0;

    while (getline(myfilein2, line))
    {
        istringstream ss(line);

        string name;
        ss >> name >> posX >> posY >> posZ;
        if (j >1) {

        mytype[j - 2] = name;
     }

        j++; 
    }

return mytype;
}
#endif // 
