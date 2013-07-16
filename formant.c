/*	PORTING OF FORMAT.M
	this routine extracts formant frequencies
	from a speech sample file 

	coded by: Francesco D'Este (Universitiet Leiden)
*/

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <sndfile.h>
#include "formant.h"
#include "shorten.h"
#include "mconf.h"
#include "global.h"



# define E_BITS_PER_COEF (2 + LPCQUANT)
# define VERY_SMALL      (1e-36)
#define MAXN 20
#define BLOCK_SIZE 512
#define FEPS                     (float)2.220446049250313e-016
#define pi (double)3,14159265
#define MAX_WINDOWS 5000



int polrt( double *xcof, double *cof, int m, float *real, float *imag )
{
register double *p, *q;
int i, j, nsav, n, n1, n2, nroot, iter, retry;
int final;
double mag = 0; 
double cofj = 0;
cmplx x0, x, xsav, dx, t, t1, u, ud;

final = 0;
n = m;
if( n <= 0 )
	return(1);
if( n > 36 )
	return(2);
if( xcof[m] == 0.0 )
	return(4);

n1 = n;
n2 = n;
nroot = 0;
nsav = n;
q = &xcof[0];
p = &cof[n];
for( j=0; j<=nsav; j++ )
	*p-- = *q++;	/*	cof[ n-j ] = xcof[j];*/
xsav.r = 0.0;
xsav.i = 0.0;

nxtrut:
x0.r = 0.00500101;
x0.i = 0.01000101;
retry = 0;

tryagn:
retry += 1;
x.r = x0.r;

x0.r = -10.0 * x0.i;
x0.i = -10.0 * x.r;

x.r = x0.r;
x.i = x0.i;

finitr:
iter = 0;

while( iter < 1000 )
{
u.r = cof[n];
if( u.r == 0.0 )
	{		/* this root is zero */
	x.r = 0;
	n1 -= 1;
	n2 -= 1;
	goto zerrut;
	}
u.i = 0;
ud.r = 0;
ud.i = 0;
t.r = 1.0;
t.i = 0;
p = &cof[n-1];
for( i=0; i<n; i++ )
	{
	t1.r = x.r * t.r  -  x.i * t.i;
	t1.i = x.r * t.i  +  x.i * t.r;
	cofj = *p--;		/* evaluate polynomial */
	u.r += cofj * t1.r;
	u.i += cofj * t1.i;
	cofj = cofj * (i+1);	/* derivative */
	ud.r += cofj * t.r;
	ud.i -= cofj * t.i;
	t.r = t1.r;
	t.i = t1.i;
	}

mag = ud.r * ud.r  +  ud.i * ud.i;
if( mag == 0.0 )
	{
	if( !final )
		goto tryagn;
	x.r = xsav.r;
	x.i = xsav.i;
	goto findon;
	}
dx.r = (u.i * ud.i  -  u.r * ud.r)/mag;
x.r += dx.r;
dx.i = -(u.r * ud.i  +  u.i * ud.r)/mag;
x.i += dx.i;
if( (fabs(dx.i) + fabs(dx.r)) < 1.0e-6 )
	goto lupdon;
iter += 1;
}	/* while iter < 500 */

if( final )
	goto lupdon;
if( retry < 5 )
	goto tryagn;
return(3);

lupdon:
/* Swap original and reduced polynomials */
q = &xcof[nsav];
p = &cof[0];
for( j=0; j<=n2; j++ )
	{
	cofj = *q;
	*q-- = *p;
	*p++ = cofj;
	}
i = n;
n = n1;
n1 = i;

if( !final )
	{
	final = 1;
	if( fabs(x.i/x.r) < 1.0e-4 )
		x.i = 0.0;
	xsav.r = x.r;
	xsav.i = x.i;
	goto finitr;	/* do final iteration on original polynomial */
	}

findon:
final = 0;
if( fabs(x.i/x.r) >= 1.0e-5 )
	{
	cofj = x.r + x.r;
	mag = x.r * x.r  +  x.i * x.i;
	n -= 2;
	}
else
	{		/* root is real */
zerrut:
	x.i = 0;
	cofj = x.r;
	mag = 0;
	n -= 1;
	}
/* divide working polynomial cof(z) by z - x */
p = &cof[1];
*p += cofj * *(p-1);
for( j=1; j<n; j++ )
	{
	*(p+1) += cofj * *p  -  mag * *(p-1);
	p++;
	}

setrut:
real[nroot] = (float)x.r;
imag[nroot] = (float)x.i;
nroot += 1;
if( mag != 0.0 )
	{
	x.i = -x.i;
	mag = 0;
	goto setrut;	/* fill in the complex conjugate root */
	}
if( n > 0 )
	goto nxtrut;
return(0);
}

//	levdurb(lpc,reflex,r,(order));
void levdurb(
       float *a,       /* (o) lpc coefficient vector starting
                              with 1.0 */
      float *k,       /* (o) reflection coefficients */
      float *r,       /* (i) autocorrelation vector */
      int order       /* (i) order of lpc filter */
   ){
       float  sum, alpha;
       int     m, m_h, i;

       a[0] = 1.0;

       if (r[0] < FEPS) { /* if r[0] <= 0, set LPC coeff. to zero */
           for (i = 0; i < order; i++) {
               k[i] = 0;
               a[i+1] = 0;
           }
       } else {
           a[1] = k[0] = -r[1]/r[0];
           alpha = r[0] + r[1] * k[0];
           for (m = 1; m < order; m++){
               sum = r[m + 1];
               for (i = 0; i < m; i++){
                   sum += a[i+1] * r[m - i];
               }





               k[m] = -sum / alpha;
               alpha += k[m] * sum;
               m_h = (m + 1) >> 1;
               for (i = 0; i < m_h; i++){
                   sum = a[i+1] + k[m] * a[m - i];
                   a[m - i] += k[m] * a[i+1];
                   a[i+1] = sum;
               }
               a[m+1] = k[m];
           }
       }


   }



void autocorr(
	//check this!!!
       float *r,       /* (o) autocorrelation vector */
       float *x, /* (i) data vector */
       int N,          /* (i) length of data vector */
       int order       /* largest lag for calculated
                          autocorrelations */
  ){
       int     lag, n;
       float   sum;

       for (lag = 0; lag <= order; lag++) {
           sum = 0;
           for (n = 0; n < N - lag; n++) {
               sum += x[n] * x[n+lag];
           }
           r[lag] = sum;

       }

//	for (lag = 0; lag < N ; lag ++) 	printf("[D] Autocorr s[%d] -> %f :\n",lag,r[lag]);
	
   }


int count_nonzero(float *real, int N){
	int count = 0,i = 0;
	for(i = 0; i < N - 1; i++) {
		if (real[i] > (float)0.01) { count++;
		 	//printf("\n[*] %f is more then 0.001, total count: %d\n",real[i],count);
		}
	}
	int non_zero = count;
	return non_zero;
}


void filter_solutions(float *real, float *imag, float *real_f, float *imag_f, int N){
	//printf("[*] filtering Soltions...\n");

	int count_b = 0,i = 0;
	for(i = 0; i < N - 1; i++) {
		if (real[i] > (float)0.01){
			real_f[count_b] = real[i];
			imag_f[count_b] = imag[i];
			//printf("%f is more then 0.001, imag_f is %f. total count: %d\n",real_f[count_b],imag_f[count_b],count_b);

			count_b++;

		}
	}
}


/* swap:  interchange v[i] and v[j]. */
void swap(double *v, int i, int j)
{
     double temp;

     temp = v[i];
     v[i] = v[j];
     v[j] = temp;
}


/* quicksort: sort v[0]..v[n-1] into increasing order. */
void quicksort(double *v, int n)
{
     int i, last;
     if (n <= 1)                         /* nothing to do */
         return;
     swap(v,0,rand() % n);       /* move pivot element to v[0] */
     last = 0;
     for (i = 1; i < n ; i++)         /* partition */
           if (v[i] < v[0])
               swap(v,++last, i);
     swap(v, 0, last);                 /* restore pivot */
     quicksort(v,last);               /* recursively sort each part. */
     quicksort(v+last+1, n-last-1);
}


// remove exact matches duplicates. Returns the new size of ffreq
int remove_duplicates(double *ffreq, int N){
	int i,j;
	int count = 0;
	double dummy[N];
	
	//printf("[H] Looking for duplicates\n");
	for (i = 0; i < N ; i++){
		for (j = 0; j < N ; j++) {
			if (ffreq[i] == ffreq[j] && i != j && ffreq[i] != (double)0 ){
				ffreq[j] = 0;
				count++;
			}
		}
	}

	int new_size = N -count;
	//printf("[O] new size: %d, count %d, N %d\n",new_size,count,N);
	count = 0 ;
	
	//copy to dummy
	
	//printf("[H] copying to dummy\n");
	for (i = 0 ; i < N ; i++){
		if (ffreq[i] != 0){
			dummy[count] = ffreq[i];
			count++;
		}
	}

	
	//riallochiamo ffreq;
	//printf("[H] Looking for duplicates\n");
	
	double *dummer;
	dummer = malloc(new_size*sizeof(*dummer));
	
	//ffreq = (double*)realloc(ffreq,new_size*sizeof(*ffreq));

	//re copy
	for (i = 0 ; i < new_size ; i++){
		dummer[i] = dummy[i];
	}

	free(ffreq);	
	ffreq = dummer;	

	return new_size;
}



/*	FUNCTION FORMANTS
	(i) : buf -> speech signal buffer
	(i) : N   -> sample number
	(o) : out -> formant vector
*/

void formants(float *signal,int N, double* out)
{

	double *formants;
	double *scaled;
	float *real_f, *imag_f;
	double *ffreq;
	float *lpc, *reflex, *real, *imag;
	double *d_lpc, *s_lpc;
	int fs= FS;
	int ms10 = fs*.05;
	int ms30 = fs*.006;
	int ncoeff = 2+fs/1000;
	int pos = 0;
	float *y;
	float *r;
	float y_mean;
	y_mean = 0;
	int non_zero = 0;
	int i = 0, j = 0;
	int iter = 0;
	int order = ncoeff;
	int new_size = 0;

	//structure to keep the story
	double freqs[MAX_WINDOWS][FORMANTS];
	// clean it up!
	for (i = 0; i < MAX_WINDOWS ; i++){
			for (j = 0; j < FORMANTS ; j++){
				freqs[i][j] = (double)0;
			}
	}

      	y = malloc(ms10*sizeof(*y));
	r = malloc(ms10 * sizeof(*r));

	if ((d_lpc = malloc(order * sizeof(*d_lpc))) == NULL){
		printf("\n\nERROR ALLOCATING d_lpc\n");
		fflush(stdout);
		exit(1);
	}

	if ((s_lpc = malloc(order * sizeof(*s_lpc))) == NULL){
		printf("\n\nERROR ALLOCATING s_lpc\n");
		fflush(stdout);
		exit(1);
	}
	if ((lpc =  malloc(order * sizeof(*lpc))) == NULL ){
		printf("\n\nERROR ALLOCATING LPC buffer\n");
		fflush(stdout);
		exit(1);;
	}
	else {
		printf("\n\n LPC COEFFS ALLOCATED\n");
	}
	if ((reflex = malloc(order * sizeof(*reflex))) == NULL){
		printf("\n\nERROR ALLOCATING reflex\n");
		fflush(stdout);
		exit(1);;
	}
	else {
		printf("\n\n REFLEX COEFFS ALLOCATED\n");
	}
	if ((real =  malloc(order * sizeof(*real))) == NULL){
		printf("\n\nERROR ALLOCATING real\n");
		fflush(stdout);
		exit(1);
	}
	if ((imag =  malloc(order * sizeof(*imag))) == NULL){
		printf("\n\nERROR ALLOCATING imag\n");
		fflush(stdout);
		exit(1);
	}

	//initialize vectors
	for (i = 0; i < order ; i++){
		reflex[i] = 0;
		lpc[i] = 0;
		d_lpc[i] = 0;
		s_lpc[i] = 0;
	}
	
	
	while ((pos+ms30) <= N) {
		
		printf("[!] Processing window s[%d,%d]\n",pos,(pos+ms10));
		
		//get a window of speech
		for (i = 0; i < (ms10) ; i++){
			y[i] = signal[pos+i];
			//printf("\n [ ] Copyng subsignal : %f",y[i]);
		}

		//compute autocorrelation check windows! (lag)
		autocorr(r,y,ms10,ms10);
		
		//find LP filter (to be implemented) (datatype???)


		//was ncoeff
		levdurb(lpc,reflex,r,(order-1));
	
		//cast lpc values to double

		for (i = 0; i < order ; i++){
			//printf("\n writing double buffer lpc : i = %d, order = %d",i,order);
			//fflush(stdout);
			d_lpc[i] = (double)lpc[i];
			s_lpc[i] = (double)lpc[i];
		}
		polrt(d_lpc, s_lpc,(order-1),real,imag);

		//find polynomial roots 
		non_zero = count_nonzero(real,(order-1));

		printf("\n non - zero: %d",non_zero);
		fflush(stdout);
		real_f = malloc(non_zero*sizeof(*real_f));
		imag_f = malloc(non_zero*sizeof(*imag_f));
		ffreq = malloc(non_zero*sizeof(*ffreq));
		//printf("[+] Filtering Solutions\n");
		filter_solutions(real,imag,real_f, imag_f,(order-1));

				
		
		//printf("[-] Filtered Solutions\n");
		for (i = 0; i < non_zero; i++){
			ffreq[i] = atan2((double)real_f[i],(double)imag_f[i]);
			//scale
			//ffreq[i] = (ffreq[i]*(fs/2*pi));
			ffreq[i] = ffreq[i] * fs/(2*3.14);
			//printf("[!] [%d]-th Ffreq extracted: %f Real_f: %f Img : %f\n",i, ffreq[i],real_f[i],imag_f[i]);
		}
		//printf("[-] Inverted complex Solutions\n");		
		

		free(real_f);
		free(imag_f);
		new_size = 0;

		
		for (i = 0; i < non_zero; i++){
			//printf("[!] [%d]-th Formant Frequency extracted at %f\n",i, ffreq[i]);
		}
		

		//new_size = remove_duplicates(ffreq,non_zero);
		//edit
		new_size = non_zero;

		quicksort(ffreq,non_zero - 1);

		
		for (i = 0; i < new_size; i++){
			//printf("[DB] [%d]-th Formant Frequency extracted at %f\n",i, ffreq[i]);
		}
		


		//find duplicates

		//do some smart filtering...

		//save the results for later averaging


		if (new_size > 0){
			for (i = 0; i < FORMANTS; i++){
				freqs[iter][i] = ffreq[i]; 
			}
			iter++;
		}
		free(ffreq);
	
		pos=pos+ms10;
	
	}

// compute formant freq average


	formants = malloc(FORMANTS*sizeof(*formants));
	scaled = malloc(FORMANTS*sizeof(*scaled));

	//reinit formants
	for (i = 0; i < FORMANTS ; i++){
		formants[j] = 0;
	}

	for (j = 0; j < FORMANTS ; j++){
		for (i = 0; i < iter ; i++){
			formants[j] = formants[j] + freqs[i][j];
			//printf("[T] freqs: %f\n", freqs[i][j]/iter);
		}
		//printf("[T] formants[%d]: %f normaformants: %f\n",j,formants[j],formants[j]/iter);
		scaled[j] = formants[j] / (iter);

	}

	//reinit formants
	for (i = 0; i < FORMANTS ; i++){
		formants[i] = 0;
	}


	for (i = 0; i < FORMANTS; i++){
		out[i] = scaled[i];
		if (DISPLAY_FORMANT) {
			printf("\n[!] %d-th Peak Freq is %f hz",i,scaled[i]);
		}
	}

	//free!!!
	free(y);
	free(r);
	free(real);
	free(imag);		
	free(lpc);
	free(reflex);
	free(d_lpc);
	free(s_lpc);

}


/*
int
main (int argc, char * argv [])
{


	//open wav file variables
	char 		*progname, *infilename ;
	SNDFILE	 	*infile = NULL ;
	FILE		*outfile = NULL ;
	SF_INFO	 	sfinfo ;


	double *out;
	infilename = "test.wav";
	printf("\n[*] LPC Analisys started on sample: %s", infilename);
	int channels = 1;
	float buf [BLOCK_SIZE] ;
	int k, m, readcount ;
	int fs = 11050;
	float *signal;
	signal = malloc(11078 * sizeof(*signal));
	
	float mono = 0;
	int count = 0;
	int i = 0;
	
	infile = sf_open (infilename, SFM_READ, &sfinfo);

	printf("\n[*] START DECODE\n");

	
	while ((readcount = sf_readf_float (infile, buf, BLOCK_SIZE)) > 0)
	{	for (k = 0 ; k < readcount ; k++){
			signal[count] = buf[k];
			printf("\n signal: -  %.3f", signal[count]);
			count++; 
			} ;
		} ;
    	printf("\n[!] length: %d\n",count);
	//call formants
    	out = (double *) malloc((FORMANTS)*sizeof(*out));

	formants(signal,count,out);
	//free(signal);
	//free(out);
	//sf_close(infile);

	
	return 0;
	//wav2lpc(cbuffer, DEFAULT_BLOCK_SIZE, qlpc, maxnlpc,buffer1,&resn);
}	
*/
