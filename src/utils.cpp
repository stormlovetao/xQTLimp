/*
###############################################################################################
## This file provide basic functions for imputation process ,including matrix operations and operations on LD
## This file is adapted from work of Bogdan Pasaniuc et. al [1]. 
## Reference: [1]. Pasaniuc, Bogdan, et al. "Fast and accurate imputation of summary statistics enhances evidence of functional enrichment." Bioinformatics 30.20 (2014): 2906-2914.
################################################################################################
*/
#include <string.h>
#include "utils.h"
#include  <stdio.h>
#include  <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <cmath>
#include <iostream>

using namespace std;

#define MAX_ITS 500
static double gmean_c_;
#define GEOMETRIC_MEAN(a,b) (fabs(a) > fabs(b) \
		     ? (gmean_c_ = (b)/(a), \
			fabs(a)*sqrt(1 + gmean_c_*gmean_c_)) \
		     : (0 != (b) \
			? (gmean_c_ = (a)/(b), \
			   fabs(b)*sqrt(1 + gmean_c_*gmean_c_)) \
			: (double)0.0))
#define AMAX1(a,b) ((a) > (b) ? (a) : (b))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

FILE* safe_fopen(const char* filename, const char *op) {
	FILE* f = fopen(filename, op);
	if(!f) {
		fprintf(stderr,"Error: Could not open file %s!\n",filename);
		exit(1);
	}
	return f;
}

void* safe_calloc(size_t nelements, size_t sizeof_element) {
	void* ptr = calloc(nelements, sizeof_element);
	if(!ptr) {
		size_t nbytes = nelements*sizeof_element;
		fprintf(stderr, "Error: Could not allocate ");
		fprintf(stderr, "%u bytes of memory!\n", (unsigned int)nbytes);
		exit(1);
	}
	return ptr;
}



// search for the index of the snp whose pos is specified by pos
// assume the snps in typed_snps are sorted by snp positions
bool search_by_pos(const vector<typed_snp> &typed_snps,
			pos_t pos, size_t &idx) {
	size_t num_typed_snps = typed_snps.size();
	
	size_t left = 0;
	size_t right = num_typed_snps - 1;
	bool found = false;
	
	while(left <= right) {
		idx = left+(right-left)/2;
		if(idx > right) {
			break;
		}
		if(typed_snps[idx].snp_pos == pos) {
			found = true;
			break;
		}
		else if(typed_snps[idx].snp_pos < pos) {
			left = idx + 1;
		}
		else {
			right = idx - 1;
		}
	}
	
	return found;
}

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<ref_snp> &all_snps, pos_t pos, size_t &idx) {
	size_t num_total_snps = all_snps.size();
	
	size_t left = 0;
	size_t right = num_total_snps - 1;
	bool found = false;
	
	while(left <= right) {
		// to avoid overflow
		idx = left+(right-left)/2;
		
		if(idx > right) {
			break;
		}
		// found it
		if(all_snps[idx].snp_pos == pos) {
			found = true;
			break;
		}
		// too small
		else if(all_snps[idx].snp_pos < pos) {
			left = idx + 1;
		}
		// too large
		else {
			right = idx - 1;
		}
	}
	
	return found;
}

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<genotype> &all_genos, pos_t pos, size_t &idx) {
	size_t num_total_snps = all_genos.size();
	
	size_t left = 0;
	size_t right = num_total_snps - 1;
	bool found = false;
	
	while(left <= right) {
		// to avoid overflow
		idx = left+(right-left)/2;
		
		if(idx > right) {
			break;
		}
		// found it
		if(all_genos[idx].snp_pos == pos) {
			found = true;
			break;
		}
		// too small
		else if(all_genos[idx].snp_pos < pos) {
			left = idx + 1;
		}
		// too large
		else {
			right = idx - 1;
		}
	}
	
	return found;
}

// mark snps that are typed in snp_flags to 1
// mark snps that need to convert z-scores in convert_flags to 1
size_t mark_snps(vector<typed_snp> &typed_snps, 
			const vector<ref_snp> &snp_ref_alt,vector<char>& convert_flags,
			vector<char>& impute_flags) {
	size_t num_total_snps = impute_flags.size();
	size_t nsnps = typed_snps.size();
	size_t count = 0;
	
	for(size_t i = 0; i < nsnps; i++) {
		pos_t snp_pos = typed_snps[i].snp_pos;
		
		// search for the snp
		size_t j;
		bool found = search_by_pos(snp_ref_alt, snp_pos, j);
		
		if(!found) {
			fprintf(stderr, "Error: %s  %lldnot found in encoding file!\n",
				typed_snps[i].snp_name.c_str(),typed_snps[i].snp_pos);
			exit(1);
		}
		
		size_t idx = snp_ref_alt[j].idx;
		if(idx < num_total_snps) {
			impute_flags[idx] = 0;
		}
		else {
			fprintf(stderr, "Error! Flag index out of bound!\n");
			exit(1);
		}
		
		// remember the typed snp idx
		typed_snps[i].idx = idx;
		
		string ref = snp_ref_alt[j].ref_allele;
		string alt = snp_ref_alt[j].alt_allele;
		
		// matched
		if(ref == typed_snps[i].ref_allele 
			&& alt == typed_snps[i].alt_allele) {
			// do nothing
		}
		// switched
		else if(alt == typed_snps[i].ref_allele 
			&& ref == typed_snps[i].alt_allele) {
			convert_flags[i] = 1;
			count++;
		}
		// error
		else {
			fprintf(stderr, "Error: Undefined encoding for SNP %s!\n", 
				typed_snps[i].snp_name.c_str());
			exit(1);
		}
	}
	
	return count;
}

void pdinv(double *cinv, double *coeff, int n)   {
  double *tt;
  double *p ;
  double sum ;
  int i,j, k ;

  p = (double *)malloc(n*sizeof(double));
  tt = (double *)malloc(n*n*sizeof(double));

  for(i=0;i<n;i++)  {
    p[i] = 0.0;
    for(j=0;j<n;j++)  {
      tt[i*n+j] = 0.0;
    }
  }
   
  for(i=0;i<n*n;i++)  {
    tt[i] = coeff[i];
  }

  choldc (tt, n, p) ;
 
  for (i=0; i<n; i++) {
    tt[i*n+i] = 1.0/p[i] ;
    for (j=i+1; j<n; j++) {
      sum=0.0 ;
      for (k=i; k<j; k++) {
        sum -= tt[j*n+k]*tt[k*n+i] ;
      }
      tt[j*n+i] = sum/p[j] ;
    }
  }

  for (i=0; i<n; i++)   {
    for (j=i; j<n; j++) {
     sum=0.0 ;
     for (k=j; k<n; k++) {
       sum += tt[k*n+j]*tt[k*n+i] ;
     }
     cinv[i*n+j] = cinv[j*n+i] = sum ;
    }
  }

  free(tt) ;
  free(p) ;

}

void 
cholsl (double *a, int n, double p[], double b[], double x[])

/** 
 Numerical Recipes.  Must change 
*/

{
  /*AT: Changing the code*/

 int i, k; 
  double sum; 
     
 
  for (i = 0; i < n; i++) 
    { 
      sum = b[i]; 
      for (k = i - 1; k >= 0; k--) 
        sum -= a[i*n+k] * x[k]; 
      x[i] = sum / p[i]; 
    } 
 
  for (i = (n-1); i >= 0; i--) 
    { 
      sum = x[i]; 
      for (k = i + 1; k < n; k++) 
        sum -= a[k*n+i]* x[k]; 
      x[i] = sum / p[i]; 
    }        


}

int 
choldc (double *a, int n, double p[])
{
   int i, j,k;         
  double sum; 
     
  
  for (i = 0; i < n; i++) 
    {           
      for (j = i; j < n; j++) 
        {                 
          sum = a[i*n+j]; 
          for  (k = i - 1; k >= 0; k--) 
            sum -= a[i*n+k] * a[j*n+k]; 
          if (i == j) 
            {          
	      /**                                       
              printf("zzchol %d %20.10f %9.3f\n",i, sum, a[i][i]) ; 
	      */                                  
              if (sum <= 0.0) {                  
                return -1 ; // not pos def
              }  
              p[i] = sqrt (sum); 
            }   
          else  
            { 
              a[j*n+i] = sum / p[i]; 
               
            } 
        }                                    
    } 
  
    return 1 ;
                
}

void cholesky(double *cf, double *a, int n) 
{  
  int i, j, k ;
  double *tt ;
  double *p ;
  
  p = (double *)malloc(n*sizeof(double));
  tt = (double *)malloc(n*n*sizeof(double));

  for(i=0;i<n;i++)  {
    p[i] = 0.0;
    for(j=0;j<n;j++)  {
      tt[i*n+j] = 0.0;
    }
  }

  for(i=0;i<n*n;i++)  {
    tt[i] = a[i];
  }

  choldc(tt, n, p ) ;

  for(i=0;i<n*n;i++)  {
    cf[i] = a[i];
  }

  for (i = 0; i < n; i++) {
    tt[i*n+i] = p[i] ;
    for (j=0; j <= i ; j++) {  
      k = (i)*n+(j) ;
      cf[k] = tt[i*n+j] ;
    }
  }
  
  free(tt) ; 
  free(p) ;
}

int svdecomp(double **u, int m, int n, double *w, double **v)
{
  int flag, i, its, j, jj, k, l, l1, k1;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  double *rv1 = NULL;
  
  if (m < n) 
  { fprintf(stderr,"svdecomp: matrix is row deficient\n");
    return -1;
  }
  
  rv1 = (double *)calloc(n, sizeof(double));
  if (NULL == rv1) 
  { perror("svdecomp"); return -1; }
  
  for (i = 0; i < n; i++)
  { l = i + 1;
    rv1[i] = scale*g;

/*
 *	Householder reduction to bidiagonal form
 */
    g = s = scale = 0.0;
    if (i < m)
    { for (k = i; k < m; k++)
      { scale += fabs(u[k][i]);
      }
      if (0 != scale)
      { for (k = i; k < m; k++)
	{ u[k][i] /= scale;
	  s += u[k][i]*u[k][i];
	}
	f = u[i][i];
	g = -SIGN(sqrt(s), f);
	h = f*g - s;
	u[i][i] = f - g;
	if (i != n-1)
	{ for (j = l; j < n; j++)
	  { s = 0.0;
	    for (k = i; k < m; k++)
	    { s += u[k][i]*u[k][j];
	    }
	    f = s/h;
	    for (k = i; k < m; k++)
	    { u[k][j] += f*u[k][i];
	    }
	  }
	}
	for (k = i; k < m; k++)
	{ u[k][i] *= scale;
	}
      }
    }

    w[i] = scale * g;
    g = s = scale = 0.0;
    if (i < m && i != n-1)
    { for (k = l; k < n; k++)
      { scale += fabs(u[i][k]);
      }
      if (0 != scale)
      { for (k = l; k < n; k++) 
	{ u[i][k] /= scale;
	  s += u[i][k]*u[i][k];
	}
	f = u[i][l];
	g = -SIGN(sqrt(s),f);
	h = f*g - s;
	u[i][l] = f - g;
	for (k = l; k < n; k++)
	{ rv1[k] = u[i][k]/h;
	}
	if (i != m-1)
	{ for (j = l; j < m; j++)
	  { s = 0.0;
	    for (k = l; k < n; k++)
	    { s += u[j][k]*u[i][k];
	    }
	    for (k = l; k < n; k++)
	    { u[j][k] += s*rv1[k];
	    }
	  }
	}
	for (k = l; k < n; k++)
	{ u[i][k] *= scale;
	}
      }
    }

    anorm = AMAX1(anorm,(fabs(w[i]) + fabs(rv1[i])));
  }

/*
 *	accumulation of right-hand transformations
 */
  for (i = n-1; i >= 0; i--)
  { if (i < n-1)
    { if (0 != g) 
      { for (j = l; j < n; j++)
	{ v[j][i] = u[i][j] / u[i][l] / g;
	}
	for (j = l; j < n; j++) 
	{ s = 0.0;
	  for (k = l; k < n; k++)
	  { s += u[i][k]*v[k][j];
	  }
	  for (k = l; k < n; k++)
	  { v[k][j] += s*v[k][i];
	  }
	}
      }
      for (j = l; j < n; j++)
      { v[i][j] = v[j][i] = 0.0;
      }
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

/*
 *	accumulation of left-hand transformations
 */
  for (i = n-1; i >= 0; i--)
  { l = i + 1;
    g = w[i];
    if (i < n-1)
    { for (j = l; j < n; j++)
      { u[i][j] = 0.0;
      }
    }
    if (0 != g)
    { g = 1/g;
      if (i != n-1)
      { for (j = l; j < n; j++)
	{ s = 0.0;
	  for (k = l; k < m; k++)
	  { s += u[k][i]*u[k][j];
	  }
	  f = (s/u[i][i])*g;
	  for (k=i; k < m; k++)
	  { u[k][j] += f*u[k][i];
	  }
	}
      }
      for (j = i; j < m; j++)
      { u[j][i] *= g;
      }
    } 
    else
    { for (j = i; j < m; j++)
      { u[j][i] = 0.0;
      }
    }
    ++u[i][i];
  }

/*
 *	diagonalization of the bidiagonal form
 */
  for (k = n-1; k > 0; k--)
  { for (its = 1; its <= MAX_ITS; its++)
    { 
/*
 *	test for splitting
 */
      flag = 1;
      for (l = k; l >= 0; l--)
      { l1 = l - 1;
	if (fabs(rv1[l])+anorm == anorm)
	{ flag = 0;
	  break;
	}
	if (fabs(w[l1])+anorm == anorm) break;
      }
/*
 *	cancellation of rv1[l] for l greater than 0
 */
      if (flag)
      { c = 0;
	s = 1;
	for (i = l; i <= k; i++)
	{ f = s*rv1[i];
	  if (fabs(f)+anorm != anorm)
	  { g = w[i];
	    h = GEOMETRIC_MEAN(f, g);
	    w[i] = h;
	    c = g/h;
	    s = -f/h;
	    for (j = 0;j < m; j++)
	    { y = u[j][l1];
	      z = u[j][i];
	      u[j][l1] = y*c + z*s;
	      u[j][i] = z*c - y*s;
	    }
	  }
	}
      }

/*
 *	test for convergence
 */
      z = w[k];
      if (l == k) 
      { if (z < 0.0)
	{ w[k] = -z;
	  for (j = 0; j < n; j++)
	  { v[j][k] = -v[j][k];
	  }
	}
	break;
      }

      if (its == MAX_ITS) 
      { fprintf(stderr,"svdecomp: failed to converge in %d iterations\n",
		MAX_ITS);
	free(rv1);
	return k;
      }

/*
 *	shift from bottom 2 by 2 minor
 */
      x = w[l];
      k1 = k-1;
      y = w[k1];
      g = rv1[k1];
      h = rv1[k];
      f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
      g = GEOMETRIC_MEAN(f, 1.0);
      f = ((x - z)*(x + z) + h*((y / (f + SIGN(g, f))) - h)) / x;
      c = s = 1;
      for (j = l; j <= k1; j++)
      {	i = j + 1;
	g = rv1[i];
	y = w[i];
	h = s*g;
	g = c*g;
	z = GEOMETRIC_MEAN(f, h);
	rv1[j] = z;
	c = f/z;
	s = h/z;
	f = x*c + g*s;
	g = g*c - x*s;
	h = y*s;
	y = y*c;
	for (jj = 0; jj < n; jj++)
	{ x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x*c + z*s;
	  v[jj][i] = z*c - x*s;
	}
	z = GEOMETRIC_MEAN(f, h);
	w[j] = z;
	if (0 != z)
	{ c = f/z;
	  s = h/z;
	}
	f = c*g + s*y;
	x = c*y - s*g;
	for (jj = 0; jj < m; jj++)
	{ y = u[jj][j];
	  z = u[jj][i];
	  u[jj][j] = y*c + z*s;
	  u[jj][i] = z*c - y*s;
	}
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  free(rv1);
  return 0;
}

void pinv(double *A_inv, double *A, size_t n) {
	// allocate memory for matrix inverse
	double **a = (double**)malloc(n*sizeof(double*));
	if(!a) {
		fprintf(stderr, "Error: Failed to allocate memory in pinv\n");
		exit(1);
	}
	double **v = (double**)malloc(n*sizeof(double*));
	if(!v) {
		fprintf(stderr, "Error: Failed to allocate memory in pinv\n");
		exit(1);
	}
	for(size_t i = 0; i < n; i++) {
		a[i] = (double*)malloc(n*sizeof(double));
		if(!a[i]) {
			fprintf(stderr, "Error: Failed to allocate memory in pinv\n");
			exit(1);
		}
		v[i] = (double*)malloc(n*sizeof(double));
		if(!v[i]) {
			fprintf(stderr, "Error: Failed to allocate memory in pinv\n");
			exit(1);
		}
	}
	
	// copy data
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < n; j++) {
			a[i][j] = A[i*n+j];
		}
	}

	// compute svd of a
	double *w = (double*)malloc(n*sizeof(double));
	if(!w) {
		fprintf(stderr, "Error: Failed to allocate memory in pinv\n");
		for(size_t i = 0; i < n; i++) {
			for(size_t j = 0; j < n; j++) {
				printf("%lf ", A[i*n+j]);
			}
			printf("\n");
		}
		exit(1);
	}
	
	// apply svd
	int rc = svdecomp(a, n, n, w, v);
	if(rc) {
		fprintf(stderr, "Error: Failed to apply SVD in pinv\n");
		printf("mat = [");
                for(size_t i = 0; i < n; i++) {
                        for(size_t j = 0; j < n; j++) {
                                printf("%lf ", A[i*n+j]);
                        }
                        printf(";\n");
                }
		printf("];/");
                exit(1);

	}
	
	// transpose a and v
	for(size_t i = 0; i < n; i++) {
		for(size_t j = i; j < n; j++) {
			double tmp;
			tmp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = tmp;
		}
	}
	
	// compute w inverse
	for(size_t i = 0; i < n; i++) {
		if(w[i] != 0.0) {
			w[i] = 1.0/w[i];
		}
	}
	
	// compute v*diagonal matrix
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < n; j++) {
			v[i][j] *= w[j];
		}
	}
	
	// compute v*diag*u
	for(size_t i = 0; i < n; i++) {
		for(size_t j = i; j < n; j++) {
			double sum = 0.0;
			for(size_t k = 0; k < n; k++) {
				sum += v[i][k]*a[k][j];
			}
			A_inv[i*n+j] = sum;
			A_inv[j*n+i] = sum;
		}
	}
	
	// clean up
	free(w);
	for(size_t i = 0; i < n; i++) {
		free(a[i]);
		free(v[i]);
	}
	free(a);
	free(v);
}

void jacobi(double *a, int n, double d[], double *v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=(double*)malloc(n*sizeof(double));
	b = b - 1;
	z=(double*)malloc(n*sizeof(double));
	z = z - 1;
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip*n+iq]=0.0;
		v[ip*n+ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip*n+ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip*n+iq]);
		}
		if (sm == 0.0) {
			z = z + 1;
			free(z);
			b = b + 1;
			free(b);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip*n+iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip*n+iq]=0.0;
				else if (fabs(a[ip*n+iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip*n+iq])/h;
					else {
						theta=0.5*h/(a[ip*n+iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip*n+iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip*n+iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	fprintf(stderr, "Too many iterations in routine jacobi\n");
}

void pinv_jacobi(double *A_inv, double *A, size_t n) {
	double *tmp_A = A - n - 1;
    
	// allocate memory for matrix inverse
	double *v = (double*)malloc(n*n*sizeof(double));
	double *tmp_v = v - n - 1;

	// compute eigen decomposition of a
	double *w = (double*)malloc(n*sizeof(double));
	double *tmp_w = w - 1;
	int nrot;
	jacobi(tmp_A, n, tmp_w, tmp_v, &nrot);
	
	// compute a*diagonal matrix
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < n; j++) {
			A[i*n+j] = v[i*n+j]*(1.0/w[j]);
		}
	}
	
	// compute a*diag*u
	for(size_t i = 0; i < n; i++) {
		for(size_t j = i; j < n; j++) {
			double sum = 0.0;
			for(size_t k = 0; k < n; k++) {
				sum += A[i*n+k]*(v[j*n+k]);
			}
			A_inv[i*n+j] = sum;
			A_inv[j*n+i] = sum;
		}
	}
	
	// clear up
	free(w);
	free(v);
}

// load z-score typed snps file
size_t load_zscore_typed_snps(const vector<typed_snp>& typed_snps,
 vector<zscore_typed_snp> &ztyped_snps
,vector<long long int>* p_useful_typed_snps) {
	size_t idx = 0;
	int len = (*p_useful_typed_snps).size();
	int typed_len = typed_snps.size();
	for(int i = 0;i < len;i++)
	{
		zscore_typed_snp tmp;
		tmp.idx = idx;
		for(int j = 0;j < typed_len;j++)
		{
			if((*p_useful_typed_snps)[i] == typed_snps[j].snp_pos)
			{
				tmp.snp_name = typed_snps[j].snp_name;
				tmp.snp_pos = typed_snps[j].snp_pos;
				tmp.zscore = typed_snps[j].zscore;
				ztyped_snps.push_back(tmp);
				idx++;
				break; 
			}
		}
		
	}
	
	size_t num_ztyped_snps = ztyped_snps.size();
	
	return num_ztyped_snps;
}


// compute allele frequencies for all the snps
void get_all_freqs(const vector<string>& haps, vector<double> &freqs) {
	size_t num_total_snps = freqs.size();
	for(size_t i = 0; i < num_total_snps; i++) {
		freqs[i] = get_freq(haps, i);
	}
}

// compute the allele frequency of a snp, specified by snp_idx
double get_freq(const vector<string>& haps, size_t snp_idx) {
	double freq = 0.0;
	size_t nhaps = haps[0].size();
	size_t i;
	for(i = 0; i < nhaps; i++) {
		if(haps[snp_idx][i] == '1') {
			freq += 1.0;
		}
	}
	freq = freq/((double)nhaps);
	return freq;
}

// compute the h frequency of two snps, specified by snp_idx1, snp_idx2
double get_h_freq(const vector<string>& haps, 
			size_t snp_idx1, size_t snp_idx2) {
	double freq = 0.0;
	size_t nhaps = haps[0].size();
	size_t i;
	for(i = 0 ; i < nhaps; i++) {
		if(haps[snp_idx1][i] == '1' && haps[snp_idx2][i] == '1') {
			freq += 1.0;
		}
	}
	freq = freq/((double)nhaps);
	return freq;
}

double get_var(const vector<string>& haps,
		const vector<double>& freqs, int idx) {
	double mn = freqs[idx];
	double vr = 0.0;
	size_t nhaps = haps.size();
	for(size_t i = 0; i < nhaps; i++) {
		if(haps[i][idx] == '1') {
			vr += (1.0-mn)*(1.0-mn);
		}
		else {
			vr += (0.0-mn)*(0.0-mn);
		}
	}
	return (vr/((double)nhaps-1.0));
}


