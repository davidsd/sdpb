// Compare independently rounded MPFR reference caches.
// Usage: compare_references rows cache-a cache-b [required_log2_error]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>
static double lg(const mpfr_t x) {
  if(mpfr_zero_p(x)) return -1e9;
  long e; double d=mpfr_get_d_2exp(&e,x,MPFR_RNDN);
  return log2(fabs(d))+e;
}
int main(int argc,char **argv) {
  if(argc<4) return 2;
  int n=atoi(argv[1]),ca,cb,count=0,rc,*ids=NULL;
  if(n<=0) return 2;
  FILE *a=fopen(argv[2],"r"),*b=fopen(argv[3],"r");
  if(!a||!b) return perror("cache"),2;
  mpfr_t x,y,d,mx,md;mpfr_inits2(2048,x,y,d,mx,md,(mpfr_ptr)0);
  double worst=-1e9;
  while((rc=fscanf(a," %d",&ca))==1) {
    if(fscanf(b," %d",&cb)!=1 || ca!=cb) return fprintf(stderr,"column mismatch\n"),2;
    if(ca<0) return 2;
    for(int j=0;j<count;j++) if(ids[j]==ca) return fprintf(stderr,"duplicate column\n"),2;
    int *next=realloc(ids,(size_t)(count+1)*sizeof(int));
    if(!next) return 2;
    ids=next; ids[count]=ca;
    mpfr_set_zero(mx,1);mpfr_set_zero(md,1);
    for(int i=0;i<n;i++) {
      if(!mpfr_inp_str(x,a,16,MPFR_RNDN)||!mpfr_inp_str(y,b,16,MPFR_RNDN)
         || !mpfr_number_p(x) || !mpfr_number_p(y)) return 2;
      mpfr_sub(d,x,y,MPFR_RNDN);mpfr_abs(d,d,MPFR_RNDN);
      if(mpfr_cmp(d,md)>0) mpfr_set(md,d,MPFR_RNDN);
      mpfr_abs(y,y,MPFR_RNDN);if(mpfr_cmp(y,mx)>0) mpfr_set(mx,y,MPFR_RNDN);
    }
    double e=mpfr_zero_p(md)?-1e9:lg(md)-lg(mx);
    worst=fmax(worst,e);printf("REFERENCE_COL %d log2_error %.9f\n",ca,e);count++;
  }
  if(rc!=EOF || !count || fscanf(b," %d",&cb)!=EOF) return 2;
  printf("REFERENCE_CONVERGENCE columns %d worst_log2 %.9f\n",count,worst);
  int fail=argc>4 && worst>atof(argv[4]);
  free(ids);fclose(a);fclose(b);mpfr_clears(x,y,d,mx,md,(mpfr_ptr)0);
  return fail ? 1 : 0;
}
