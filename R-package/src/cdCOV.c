//#include "stdafx.h"
#include "math.h"
#include "stdlib.h" 
#include "R.h"
#include "Rmath.h"
#include "stdio.h"


double **alloc_matrix(int r, int c)
{
    /* allocate a matrix with r rows and c columns */
    int i;
    double **matrix;
    matrix = (double **) calloc(r, sizeof(double *));
    for (i = 0; i < r; i++)
    matrix[i] = (double *) calloc(c, sizeof(double));
    return matrix;
}


void free_matrix(double **matrix, int r, int c)
{
    /* free a matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++){
		free(matrix[i]);	
	}
    free(matrix);
}


double brdet(double *a,int n)
{
	/*computes determinant of a matrix*/
	int i,j,k,is,js,l,u,v;
	double f,det,q,d;
	is=0;
	js=0;
	f=1.0;
	det=1.0;
	for(k=0; k<=n-2; k++){
	   q=0.0;
	   for (i=k; i<=n-1; i++)
		   for (j=k; j<=n-1; j++){
			   l=i*n+j;
			   d=fabs(a[l]);
			   if (d>q){
				   q=d;
				   is=i;
				   js=j;
			   }
		   }
		if (q+1.0==1.0){
			det=0.0;
			return(det);
		}
		if (is!=k){
			f=-f;
			for (j=k;j<=n-1;j++){
				u=k*n+j;
				v=is*n+j;
				d=a[u];
				a[u]=a[v]; 
				a[v]=d;
			}
		}
		if (js!=k)
		{ 
			f=-f;
			for (i=k; i<=n-1; i++){
				u=i*n+js; 
				v=i*n+k;
				d=a[u]; 
				a[u]=a[v];
				a[v]=d;
			}
		}
		l=k*n+k;
		det=det*a[l];
		for (i=k+1; i<=n-1; i++){
			d=a[i*n+k]/a[l];
			for (j=k+1; j<=n-1; j++){
				u=i*n+j;
				a[u]=a[u]-d*a[k*n+j];
			}
		}
	}
	det=f*det*a[n*n-1];
	return(det);
}


int brinv(double a[], int n)
{ 
	/*computes Inverse of a matrix*/
	int *is,*js,i,j,k,l,u,v;
    double d,p;
    is=(int *)malloc(n*sizeof(int));
    js=(int *)malloc(n*sizeof(int));
    for (k=0; k<=n-1; k++){ 
		d=0.0;
        for (i=k; i<=n-1; i++)
			for (j=k; j<=n-1; j++){ 
				l=i*n+j; p=fabs(a[l]);
				if (p>d) { 
					d=p; is[k]=i; js[k]=j;
			    }
            }
        if (d+1.0==1.0){ 
			free(is); free(js); //printf("err**not inv\n");
            return(0);
        }
        if (is[k]!=k)
			for (j=0; j<=n-1; j++){ 
				u=k*n+j; v=is[k]*n+j;
				p=a[u]; a[u]=a[v]; a[v]=p;
			}
        if (js[k]!=k)
			for (i=0; i<=n-1; i++){ 
				u=i*n+k; v=i*n+js[k];
				p=a[u]; a[u]=a[v]; a[v]=p;
			}
        l=k*n+k;
        a[l]=1.0/a[l];
        for (j=0; j<=n-1; j++)
			if (j!=k){ 
				u=k*n+j; 
				a[u]=a[u]*a[l];
			}
        for (i=0; i<=n-1; i++)
			if (i!=k)
			for (j=0; j<=n-1; j++)
				if (j!=k){ 
					u=i*n+j;
					a[u]=a[u]-a[i*n+k]*a[k*n+j];
				}
        for (i=0; i<=n-1; i++)
			if (i!=k){ 
				u=i*n+k; a[u]=-a[u]*a[l];
			}
    }
    for (k=n-1; k>=0; k--){ 
		if (js[k]!=k)
			for (j=0; j<=n-1; j++){ 
				u=k*n+j; v=js[k]*n+j;
				p=a[u]; a[u]=a[v]; a[v]=p;
			}
        if (is[k]!=k)
			for (i=0; i<=n-1; i++){ 
				u=i*n+k; v=i*n+is[k];
				p=a[u]; a[u]=a[v]; a[v]=p;
			}
    }
    free(is); 
	free(js);
    return(1);
}


void brmul(double *a, double *b,int m,int n,int k,double *c)
{ 
	/*Matrix multiplication, 
	  interpret a as m rows and n columns matrix
	  interpret b as n rows and r columns matrix*/
	int i,j,l,u;
    for (i=0; i<=m-1; i++) 
		for (j=0; j<=k-1; j++){ 
			u=i*k+j; c[u]=0.0;
			for (l=0; l<=n-1; l++)
			  c[u]=c[u]+a[i*n+l]*b[l*k+j];
		}
    return;
}


/*
 * 
 * Compute the kernel estimation for C_i - C_j
 * 
 * z: a vectorized matrix. z[i*d+k] is corrresponding to the i-th row, k-th column of original matrix C, 
 * C is a n \times d matrix.
 * Kernel: a vectorized kernel density estimator matrix. Kernel[i*n+j] = Kernel[j*n+i] are the density
 * estimator of C_i - C_j
 * 
 */
void dmvnorm(double *z, double *width, int d, int n, double *Kernel)
{
	/*Multivariate normal density*/
	int i, j, k;
	double density, det, *molecular1, *molecular, *expect, *sigma;
	molecular1 = (double *)calloc(d, sizeof(double));
	molecular = (double *)malloc(sizeof(double));
	expect = (double *)calloc(d, sizeof(double));
	sigma = (double *)calloc(d*d, sizeof(double));

	if(d==1){
		*sigma = (*width)*(*width);
		det = *sigma;
	}
	else{
		for(i=0; i<d*d; i++)
			sigma[i] = width[i];
		det = brdet(width, d);
	}
	
	brinv(sigma, d);
	density = 1.0/(pow(2*PI, d/2.0)*pow(det, 0.5));

	for(i=0; i<n; i++){
		Kernel[i*n+i] = density;
		for(j=0; j<i; j++){
			for(k=0; k<d; k++){
				expect[k] = z[i*d+k] - z[j*d+k];
			}			
			brmul(expect, sigma, 1, d, d, molecular1);
			brmul(molecular1, expect, 1, d, 1, molecular);
			Kernel[i*n+j] = Kernel[j*n+i] = exp(-(*molecular)/2.0)/(pow(2*PI, d/2.0)*pow(det, 0.5));
		}
	}
    
	free(molecular1);
	free(molecular);
	free(expect);
	free(sigma);
}


void index_distance(double *x, double **Dx, int n, int d, double index)
{
    /*
        interpret x as an n by d matrix, in row order (n vectors in R^d)
        compute the Euclidean distance matrix Dx
    */
    int i, j, k, p, q;
    double dsum, dif;
    for (i=1; i<n; i++) {
        Dx[i][i] = 0.0;
        p = i*d;
        for (j=0; j<i; j++) {
            dsum = 0.0;
            q = j*d;
            for (k=0; k<d; k++) {
                dif = *(x+p+k) - *(x+q+k);
                dsum += dif*dif;
            }
            Dx[i][j] = Dx[j][i] = pow(sqrt(dsum),index);
        }
    }
}


double Akl(double **akl, double **A, int n, double *Ker) 
{
    /* -computes the A_{kl} or B_{kl} distances from the
        distance matrix (a_{kl}) or (b_{kl}) for dCov, dCor, dVar
        dCov = mean(Akl*Bkl), dVar(X) = mean(Akl^2), etc.
    */
    int j, k;
	double Ker_sum;
    double *akbar;
    double abar; 

	akbar = (double *) calloc(n, sizeof(double));
	Ker_sum = 0.0;
    abar = 0.0;
	
	for (j=0; j<n; j++) {
		Ker_sum += Ker[j];
	}

    for (k=0; k<n; k++) {
        akbar[k] = 0.0;
        for (j=0; j<n; j++) {
            akbar[k] += akl[k][j] * Ker[j];
        }
        abar += akbar[k] * Ker[k];
        akbar[k] /= Ker_sum;
    }
    abar /= Ker_sum * Ker_sum;

    for (k=0; k<n; k++)
        for (j=k; j<n; j++) {
            A[k][j] = akl[k][j] - akbar[k] - akbar[j] + abar;
            A[j][k] = A[k][j];
        }
    free(akbar);
    return Ker_sum;
}


void cdCOV(double *x, double *y, double *z, int *n, int *p,
                  int *q, int *d, double *index, double *width, double *Kernel, double *CDCOV) 
{
    /*  computes cdCov(x,y)  */

    int    i, j, k;
    double **Dx, **Dy, **A, **B;
    double Ker_sum, n2;

	dmvnorm(z, width, *d, *n, Kernel);

	/* computes Euclidean distance */
    Dx = alloc_matrix(*n, *n);
    Dy = alloc_matrix(*n, *n);
   
	index_distance(x, Dx, *n, *p, *index);
	index_distance(y, Dy, *n, *q, *index);

    for (i=0; i<*n; i++)
        CDCOV[i] = 0.0;

	A = alloc_matrix(*n, *n);
	B = alloc_matrix(*n, *n);

	for(i=0; i<*n; i++){
		Akl(Dx, A, *n, Kernel+i*(*n));
		Ker_sum = Akl(Dy, B, *n, Kernel+i*(*n));
   
		n2 = Ker_sum * Ker_sum;

		/* compute dCov(x,y) */  
		for (k=0; k<*n; k++){
			for (j=0; j<*n; j++) {
				CDCOV[i] += A[k][j]*B[k][j]*(Kernel+i*(*n))[k]*(Kernel+i*(*n))[j];
			}
		}
		CDCOV[i] /= n2; 
	}

	free_matrix(A, *n, *n);
    free_matrix(B, *n, *n);
	free_matrix(Dx, *n, *n);
    free_matrix(Dy, *n, *n);
    return;
}


void cdCOVtest(double *x) 
{
	/*int i;
	for(i=0; i<1000000000; i++){}*/
	*x = runif(0.0, 1.0);
	//printf("%f ", x);
}
