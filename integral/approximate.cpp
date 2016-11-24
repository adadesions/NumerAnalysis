/*
 * =====================================================================================
 *
 *       Filename:  approximate.cpp
 *
 *    Description:  Integral by numerical approximate method
 *
 *        Version:  1.0
 *        Created:  11/19/16 18:44:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ada Kaminkure (AdaCode), ada@adacode.io
 *        Company:  ADACODE.IO
 *
 * =====================================================================================
 */

#include <iostream>

using namespace std;

float f(float x)
{
	return x*x+2*x+1;
}

float integral(int a, int b, float (*f)(float x))
{
	int i, n;
	float cof, xi, result;
	n = 10000;
	cof = (float)(b-a)/n;
	result = 0;
	for(i = 0;i < n; i++)
	{
		xi = a+(i*cof);
		cout << "xi =" << xi << endl;
		cout << "cof =" << cof << endl;
		result = result + f(xi);
		cout << "f(" << xi << ")= " << f(xi) << endl;
	}
	return cof*result;
}

int main()
{
	cout <<" Result = "<< integral(0,1,&f) << endl;
	return 0;
}













