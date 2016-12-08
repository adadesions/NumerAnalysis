/*
 * =====================================================================================
 *
 *       Filename:  RegularFalsy.cpp
 *
 *    Description: Line search by RegularFalsy Method 
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
#include <string>
#include "stdio.h"

using namespace std;

double f(double x[], int n);
double partial(double x[], int pos);
double grad(double x[], int n);
double norm(double x[], int n);
double * u_vector(double x[], int n);
double * cal_x1(double x[], double *u, int n, double t);
double steepest(double x[], double *u, int n, double t);
double phi_prime(double t, int n);
double regularFalsi(double t[], int n);

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

double grad(double x[], int n)
{
	double sum;
	int i;
	sum = 0;
	
	for(i=0;i<n;i++)
	{
		sum = sum + partial(x,i,n);
	}
	return sum;
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

double phi_prime(double t, int n)
{
	double x0[] = {2, 3};
	double *u, *x1, result, my_norm;
	int i;
	result = 0;
	u = u_vector(x0,n);
	x1 = cal_x1(x0, u, n, t);
	for(i=0;i<n;i++)
	{
		result = result+(*(u+i)*partial(x1,i,n));
	}
	return result;
}
double regularFalsi(double t[], int n)
{
	const double esp = 0.001;
	double firstTerm, secondTerm, t_new, phi_0, phi_1, phi_2;
	phi_0 = phi_prime(t[0],n);
	phi_1 = phi_prime(t[1],n);
	firstTerm = t[1];
	secondTerm = phi_1*((t[1]-t[0])/(phi_1-phi_0));	

	// Cal t(k+1) as t_new
	t_new = firstTerm-secondTerm;
	phi_2 = phi_prime(t_new,n);
	
	// Display Value
	printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t\n",t[0],t[1],t_new,phi_1,phi_2);
	//HALT CASE
	if(phi_2 < esp)
		return t_new;

	//Assign New Value
	if(phi_1*phi_2 > 0)
		t[1] = t_new;
	else if( phi_1*phi_2 < 0 )
		t[0] = t_new;

	regularFalsi(t, n);
}


int main ()
{
	int n;
	double t[] = {3,5}, t_star;
	n = 2;
	printf("t(k-1)\tt(k)\tt(k+1)\tphi(k)\tphi(k+1)\n");
	t_star = regularFalsi(t,n);
	cout << "Answer : " << t_star << endl;
	return 0;
}


