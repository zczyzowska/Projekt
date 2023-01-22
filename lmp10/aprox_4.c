#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>

typedef struct Matrix {
	int r,c;
	double **data;
} Matrix;

void printToScreen(Matrix *mat) {
	  int i,j;
	  printf("[ \n");
	  for (i = 0; i<mat->r; i++) {
	       printf("  ");
	       for (j = 0; j < mat->c; j++) {
	            printf("%f ", mat->data[i][j]);
	       }
	       printf("; \n");
	   }
	   printf("]\n");
}
Matrix * createMatrix(int r, int c) {
	    int i;
	    Matrix * mat = (Matrix*) malloc(sizeof(Matrix));
	    if (mat != NULL) {
	    	mat->r = r;
	    	mat->c = c;
	    	mat->data = (double**) malloc(sizeof(double*) * r);
		for (i=0; i < r; i++) {
		    mat->data[i] = (double*) malloc(sizeof(double) * c);
		}
	     }
	     return mat;
}

void freeMatrix(Matrix * mat) {
	   int i;
	   for (i=0;i < mat->r; i++)
	        free(mat->data[i]);
	   free(mat->data);
	   free(mat);
}
int eliminate( Matrix *AB)
{
	for(int k = 0; k <= AB->r - 2; k++){
		int index_max_element = k;
		double max_element = fabs(AB->data[k][k]);
		for(int m = k + 1; m < AB->r; m++){
			if(fabs(AB->data[m][k]) > max_element){
				max_element = fabs(AB->data[m][k]);
				index_max_element = m;
			}
		}
		if(max_element == 0)
			return 1; /* macierz osobliwa */
		if(index_max_element != k){
		/*Zamiana wierszy */
			double *temp = AB->data[k];
			AB->data[k] = AB->data[index_max_element];
			AB->data[index_max_element] = temp;
		}	
		for(int i = k + 1; i < AB->r; i++){
			double p = AB->data[i][k]/AB->data[k][k];
			for(int j = k + 1; j <= AB->r; j++)
				AB->data[i][j] = AB->data[i][j] - p*(AB->data[k][j]);
		}
	}
	/*Sprawdzenie ostatniego elementu na głównej przekątnej */
	if(AB->data[AB->r - 1][AB->r - 1] == 0) /* macierz osobliwa */
		return 1;
	/*zerowanie macierzy */
	for(int i = 1; i < AB->r; i++){
		for(int j = 0; j < i; j++)
			AB->data[i][j] = 0;
	}
	return 0;
}
int  backsubst(Matrix *x, Matrix *AB) {

    int n=AB->r;
    if (n != AB->c-1)
         return 2;

	for( int r = n-1; r >= 0; r-- ){
		double s = 0;
		for( int c = r; c < n; c++ )
			s += AB->data[r][c]*x->data[c][0];
		 if( AB->data[r][r] == 0 )
            return 1;
		 x->data[r][0] = (AB->data[r][n] - s) / AB->data[r][r];
	}
      return 0;
 }
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
	printf("Macierz rozszerzona:\n");   
	write_matrix(B,stdout);
	if (piv_ge_solver(B)) {
		spl->n = 0;
		return;
	}            
	//res = eliminate(B);
	//z = createMatrix(B->r, 1);
	//res = backsubst(z,B);
	printf("Rozwiązanie równania (współczynniki):\n");
	write_matrix(B,stdout);
	
	fprim = make_matrix(4,1);
	fprim2 = make_matrix(3,1);
	fprim3 = make_matrix(2,1);
	for (i=1;i<5;i++)
		fprim->e[i-1] = B->e[i]*i;
	for (i=1;i<4;i++)
		fprim2->e[i-1] = fprim->e[i]*i;
	for (i=1;i<3;i++)
		fprim3->e[i-1] = fprim2->e[i]*i;
	if (alloc_spl(spl,nb) == 0) {
		for(i=0;i<nb;i++) {
			spl->x[i] = pts->x[i];
			double sum=0;
			for (j=0;j<5;j++)
				sum+=B->e[j]*pow(x[i], j);
			spl->f[i] = sum;
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

	free_matrix(B);
	free_matrix(fprim);
	free_matrix(fprim2);
	free_matrix(fprim3);
	
	

}
