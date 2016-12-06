/*
 * =====================================================================================
 *
 *       Filename:  steepest.cpp
 *
 *    Description:  Root finding implement by steepest descant method
 *
 *        Version:  1.0
 *        Created:  11/25/16 00:13:35
 *       Revision:  none
 *       Compiler: g++ 
 *
 *         Author:  Ada Kaminkure (AdaCode), ada@adacode.io
 *        Company:  ADACODE.IO *
 * =====================================================================================
 */

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include "stdio.h"

using namespace std;

double f(double x[], int n);
double partial(double x[], int pos);
double grad(double x[], int n);
double norm(double x[], int n);
double * u_vector(double x[], int n);
double * cal_x1(double x[], double *u, int n, double t);
double steepest(double x[], double *u, int n, double t);

double f(double x[], int n)
{
	double a, b, c;
	a = 4*x[0]*x[0];
	b = 4*x[0]*x[1];
	c = 2*x[1]*x[1];
	return a-b+c;
}

double partial(double x[], int pos, int n)
{
	double temp[n], h;
	int i;
	h = 0.000001;
	for(i=0;i<n;i++)
	{
		if( i == pos)
			temp[i] = x[i]+h;
		else
			temp[i] = x[i];
	}

	return (f(temp,n)-f(x,n))/h;
}

double norm(double x[], int n)
{
	double p[n], sum;
	int i;
	sum = 0;

	for(i=0;i<n;i++)
	{
		p[i] = pow(partial(x, i, n), 2);
		sum = sum + p[i];
	}

	return sqrt(sum); 
}

double * u_vector(double x[], int n)
{
	static double temp[10];
	double myNorm;
	int i;
	
	myNorm = norm(x, n);

	for(i=0;i<n;i++)
	{
		temp[i] = (-1)*partial(x, i, n)/myNorm;
	}
	return temp;
}

double * cal_x1(double x[], double *u, int n, double t)
{
	static double temp[10];
	int i;

	for(i=0;i<n;i++)
	{
		temp[i] = x[i]+t*(*(u+i));
	}
	return temp;
}

double steepest(double x[], double *u, int n, double t)
{	
	double *x1;
	x1 = cal_x1(x, u, n, t);
	return f(x1, n);
}

double goldenSec(double min, double max, int n)
{
	const double esp = 0.001;
	const double ratio = 1.618;
	double x0[] = {2, 3};
	double differRatio, t[n], *u;
	double phi_0, phi_1;
	bool halt_check1, halt_check2;

	differRatio = (max - min)/ratio;
	t[0] = max - differRatio;
	t[1] = min + differRatio;
	
	u = u_vector(x0,n);
	phi_0 = steepest(x0, u, n, t[0]);
	phi_1 = steepest(x0, u, n, t[1]);
	halt_check1 = abs(max-min) < esp;
	halt_check2 = abs(phi_1-phi_0) < esp;

	printf("%.3f\t%.3f\t%.3f\t%.3f\t", max,min,t[0],t[1]);
	printf("%.3f\t%.3f\t%.3f\t%.3f\n", phi_0, phi_1,abs(max-min),abs(phi_1-phi_0));

	if( halt_check1 )
	{	
		cout << "HALT At t0: " << t[0] << endl;
		cout << "HALT At t1: " << t[1] << endl;
		return t[0];
	}
	else if( halt_check2 )
	{
		cout << "HALT At t0: " << t[0] << endl;
		cout << "HALT At t1: " << t[1] << endl;
		return t[0];
	}
	

	if( phi_0 > phi_1 )
		min = t[0];
	else if( phi_0 < phi_1 )
		max = t[1];
	goldenSec(min, max, n);
}

int main ()
{
	double *u, *x1, *t, result;
	int n, i;
	n = 2;
	printf("Max\tMin\tt[0]\tt[1]\tphi_0\tphi_1\tAbsMM\tAbsPhi\n");
	goldenSec(0,3,n);	
	return 0;
}


