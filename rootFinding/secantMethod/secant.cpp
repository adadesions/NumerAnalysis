#include <iostream>

using namespace std;

float f(float x)
{
  return x*x-2*x+1;
}

float secant(float pos1, float pos2)
{
  float numerator, denominator;
  numerator = pos2*f(pos1) - pos1*f(pos2);
  denominator = f(pos1) - f(pos2);
  return numerator/denominator;
}

float reSecant(float pos1, float pos2)
{
  float posN;
  posN = secant(pos1, pos2);
  cout << posN << "  |  " << pos1 << "|  " << pos2 << "|" << endl; 
  if( f(posN) == 0 )
  {
    return posN;
  }
  return reSecant(pos2, posN);
}

int main()
{
  float posN, pos1, pos2;
  cout << "posN |  pos1  | pos2" << endl;
  cout << "--------------------" << endl;
  posN = reSecant(-13,10);
  cout<< "posN = "  <<  posN << endl;
  cout<< "f( "<< posN << " ) = "<< f(posN) << endl;

}
