/*
 * =====================================================================================
 *
 *       Filename:  steepest.c
 *
 *    Description:  Implement Steepest Descant, Root finding Method
 *
 *        Version:  1.0
 *        Created:  11/24/16 02:24:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ada Kaminkure (AdaCode), ada@adacode.io
 *        Company:  ADACODE.IO
 *
 * =====================================================================================
 */

#include "stdio.h"

double f(double x[])
{
	return 4*x[0]*x[0]-4*x[0]*x[1]+2*x[1]*x[1];
}
double grad(double x[], int n)
{
	double xh[n], h, partial[n];
	int i, len;
	len = n;
	h = 0.001;
	for(i=0;i<len;i++)
	{
		double c1, c2;
		c1 = x[i] + h;
		c2 = x[i] - h;
		partial[i] = (f(c1) - f(c2))/2*h;
		printf("ParX = %.6f\n", partial[i]);
	}
	return partial[0] + partial[1];
}
int main()
{
	double x[] = {1,2};
	printf("Answer is %.6f\n",f(x));
	printf("Gradient is %.6f\n",grad(x, 2));
}

















