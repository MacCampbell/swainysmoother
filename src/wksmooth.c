
#define FOR_R 1
#ifdef FOR_R
	#include <R.h>
#endif

#ifndef FOR_R
	#include <stdio.h>
	#include <math.h>
#endif
/* Function to compute a weighted smooth. Parameters/variables:
	y  :  measurements at each point (vector)
	G  :  the x-position at each point (vector)
	L  :  the weight to be given to each (G,y) pair
	z  :  x-positions at which to calculate the smooth
	
	N  :  length of y and G
	M  :  length of z
	
	W  : the sliding window width
	
	res : vector where we want to return the result.  Should be of length M

*/

#define ECA_ABS(x) ((x) >= 0 ? (x) : -(x))

void wksmooth(
		double *G,
		double *y,
		double *L,
		double *z,
		int *N,
		int *M,
		double *W,
		double *res
	) 
{
	int i,j, start=0;
	double numer, denom, diff;
	int been_there, new_start;
	
	for(j=0;j<(*M);j++)  {
		numer=0;
		denom=0;
		been_there=0;
		for(i=start;i<(*N);i++) {
			new_start=start;
//	printf("j=%d  i=%d    z[j]=%f   G[i]=%f   y[i]=%f\n", j, i, z[j], G[i], y[i]);
			if(ECA_ABS(G[i] - z[j]) <= (*W)/2) {
				numer += y[i] * L[i];
				denom += L[i];
				been_there=1;
			}
			else if(been_there==1) {
				break; /* get it out of this i-loop */
			}
			else {
				new_start=i;
			}
		}
		
		/* having been through that all, if we did hit some points within the window, we update
		   the variable start.  If we didn't then we don't, because this will cause the thing
		   to run all the way off to the end! */
		if(been_there==1) start=new_start;
		
		if(denom==0.0) {
			res[j] = NAN;
		}
		else {
			res[j] = numer/denom;
//			printf("j=%d   res[j]=%f     start=%d  max-i=%d \n",j, res[j], start,i);
		}
	}
}

#ifndef FOR_R

int main(void) {
	double G[10000];
	double y[10000];
	double L[10000];
	double z[1000];
	double res[1000];
	
	int i,j, M=1000, N=10000;
	double W=100;
	
	for(i=0;i<10000;i++) {
		G[i]=1.0*i;
		y[i]=rand()/1036765915.4 + G[i];
		L[i]=3.0;
	}
	for(j=0;j<1000;j++) {
		z[j]=10.0*j;
		res[j]=0.0;
	}
	wksmooth(G, y, L, z, &N, &M, &W, res);
	printf("Hello World!\n");	
}	
#endif
