#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemm */

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()


void permute(double *x, double *y, int *n, int *xlen, int *nperm, double *zstats, double *tmat, int *rarray)

{

int i, k, l, m;
double cumsum;
int temp;


/* Set random seed using Splus function */

RANDIN;

/* Calculate first z-statistic (unpermuted data). */

cumsum = 0;

for(k = 0; k < *xlen; k++) {
   cumsum += x[k] * y[k];
}

zstats[0] = cumsum / *xlen;

/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }

/* Convert x to a full matrix. */
   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         tmat[k * *n + l] = x[m];
         tmat[l * *n + k] = x[m];
         m++;
      }
   }


/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (int)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }


/* Reorder x and take lower triangle. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         x[m] = tmat[rarray[k] * *n + rarray[l]];
         m++;
      }
   }


/* Calculate new sum of products. */

   cumsum = 0;

   for(k = 0; k < *xlen; k++) {
         cumsum += x[k] * y[k];
   
   }

   zstats[i] = cumsum / *xlen;

}

/* Reset random seed using an Splus function. */

RANDOUT;

}







void permuteKendall(double *x, double *y, int *n, int *xlen, int *nperm, double *zstats, double *tmat, int *rarray)

{

int i, j, k, l, m;
int temp;
double paired;
double unpaired;

/* Set random seed using Splus function */

RANDIN;

/* Calculate first Kendall z-statistic (unpermuted data). */

paired = 0;
unpaired = 0;

for(k = 0; k < (*xlen-1); k++) {
	for(j = k; j < *xlen; j++) {
		if (x[k]<x[j]) 
		{
			if (y[k]<y[j]) 
			{
				paired++;
			}
			else
			{
				if (y[k]>y[j])
				{
					unpaired++;
				}
			}
		}
		if (x[k]>x[j]) 
		{
			if (y[k]>y[j]) 
			{
				paired++;
			}
			else
			{
				if (y[k]<y[j])
				{
					unpaired++;
				}
			}
		}
   }
}

zstats[0] = (paired-unpaired); // ((*xlen-(*xlen-1))/2);

/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }

/* Convert x to a full matrix. */
   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         tmat[k * *n + l] = x[m];
         tmat[l * *n + k] = x[m];
         m++;
      }
   }


/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (int)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }


/* Reorder x and take lower triangle. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         x[m] = tmat[rarray[k] * *n + rarray[l]];
         m++;
      }
   }


/* Calculate new sum of products. */

	paired = 0;
	unpaired = 0;
	
	for(k = 0; k < (*xlen-1); k++) {
		for(j = k; j < *xlen; j++) {
			if (x[k]<x[j]) 
			{
				if (y[k]<y[j]) 
				{
					paired++;
				}
				else
				{
					if (y[k]>y[j])
					{
						unpaired++;
					}
				}
			}
			if (x[k]>x[j]) 
			{
				if (y[k]>y[j]) 
				{
					paired++;
				}
				else
				{
					if (y[k]<y[j])
					{
						unpaired++;
					}
				}
			}
   		}
	}
	
	zstats[i] = (paired-unpaired); // ((*xlen-(*xlen-1))/2);

}

/* Reset random seed using an Splus function. */

RANDOUT;

}
