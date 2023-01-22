#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void make_spl(points_t*pts, spline_t *spl)
{
	matrix_t	*B = NULL;
	//Matrix     *B=NULL; //macierz rozszerzona do gausa
	//Matrix	   *z=NULL;
	matrix_t     *fprim = NULL;
	matrix_t     *fprim2 = NULL;
	matrix_t     *fprim3 = NULL;
	double      *x = pts->x;
	double      *y = pts->y;
	double      X[9];
	double      Y[5];
	int i,j,res;
	int nb = pts->n; //>10 ? 10 : pts->n;

	for (i=0;i<9;i++)
	{
	      X[i]=0;
	      for (j=0;j<nb;j++)
	          X[i] = X[i] +pow(x[j],i);
	}
	B = make_matrix( 5, 6);
	for (i=0;i<=4; i++)
	{
	      for (j=0;j<=4;j++)
	      	B->e[i*B->cn+j] = X[i+j];
	}
	for (i=0;i<=4;i++)
	{   Y[i]=0;
	    for (j=0;j<nb;j++) {
	        Y[i] = Y[i] + pow(x[j], i)*y[j];
	    }
	}
	for (i=0; i<=4;i++)
	    B->e[i*B->cn+5] = Y[i];
	#if DEBUG
	printf("Macierz rozszerzona:\n");   
	write_matrix(B,stdout);
	#endif

	if (piv_ge_solver(B)) {
		spl->n = 0;
		return;
	}            
	#if DEBUG
	printf("Rozwiązanie równania (współczynniki):\n");
	write_matrix(B,stdout);
	#endif
	fprim = make_matrix(4,1);
	fprim2 = make_matrix(3,1);
	fprim3 = make_matrix(2,1);
	int k=1;
	for (i=1;i<5;i++) {
		fprim->e[i-1] = B->e[((i+1)*5)+k]*i;
		k++; }
	for (i=1;i<4;i++)
		fprim2->e[i-1] = fprim->e[i]*i;
	for (i=1;i<3;i++)
		fprim3->e[i-1] = fprim2->e[i]*i;
	if (alloc_spl(spl,nb) == 0) {
		for(i=0;i<nb;i++) {
			spl->x[i] = pts->x[i];
			double sum=B->e[5];
			int m=1;
			for (j=1;j<5;j++) {
				sum+=B->e[((j+1)*5)+m]*pow(x[i], j);
				m++;
		}spl->f[i] = sum;
		}
		for (i=0;i<nb;i++) {
			double sum=0;
			for (j=0;j<4;j++)
				sum+=fprim->e[j]*pow(x[i], j);
			spl->f1[i] = sum;
		}
		for (i=0;i<nb;i++) {
			double sum=0;
			for (j=0;j<3;j++)
				sum += fprim2->e[j]*pow(x[i],j);
			spl->f2[i] = sum;
		}
		for (i=0; i<nb; i++) {
			double sum=0;
			for (j=0;j<2;j++)
				sum += fprim3->e[j]*pow(x[i],j);
			spl->f3[i] = sum;
		}
	}

	//free_matrix(B);
	//free_matrix(fprim);
	//free_matrix(fprim2);
	//free_matrix(fprim3);
	
	

}
