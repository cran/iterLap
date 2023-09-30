#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>
#include<Rinternals.h>

void calcQuadform(double *beta, double *Q, int *dim, 
		  int *incr, double *out){
  // calculates quadratic form beta'Qbeta
  // Q - symmetric matrix (code uses only upper triangular part)
  // dim - dimension of beta and Q
  // incr - increment (to add to index)
  int i,j,incrMat, incrVec;
  incrVec = *incr * *dim;
  incrMat = incrVec * *dim;
  for(i=0;i<*dim;i++){
    for(j=i;j<*dim;j++){
      if(i==j){
        *out += Q[*dim*i+j+incrMat]*beta[i]*beta[i];
      } else {
        *out += 2*Q[*dim*i+j+incrMat]*beta[i]*beta[j];
      }
    }
  }
}


double tdens(double *x, int *dim, double *m, double *SigmaInv,
	   double *det, int *nu, double *work, double *lgamdif,
           double *dimlogpidf, int *incr){
  int i=0;
  double out=0.0;
  int incrVec;
  incrVec = *incr * *dim;
  for(i=0;i<*dim;i++){
    work[i] = x[i] - m[i+incrVec];
  }
  calcQuadform(work, SigmaInv, dim, incr, &out);
  out = -0.5*(*nu + *dim)*log(1 + out/ *nu);
  out = *lgamdif - 0.5*(log(*det) + *dimlogpidf) + out;
  return exp(out);
}


void calcmixtdens(double *x, int *dim, int *ncomp, double *m, double *SigmaInv,
		  double *dets, double *weights, int *nu, double *work,
		  double *out){
  int i=0;
  double lgamdif = 0.0;
  double dimlogpidf = 0.0;
  double temp = 0.0;
  lgamdif = lgammafn(((double)*nu + (double)*dim)/ 2) - lgammafn((double)*nu / 2);
  dimlogpidf = *dim * log((double)*nu * M_PI);
  for(i=0;i<*ncomp;i++){
    temp = tdens(x, dim, m, SigmaInv, &dets[i], nu,
		 work, &lgamdif, &dimlogpidf, &i);
    *out += weights[i]*temp;
  }
}

double mvdens(double *x, int *dim, double *m, double *SigmaInv,
	   double *det, double *work, double *log2pi, int *incr){
  int i=0;
  double out=0.0;
  int incrVec;
  incrVec = *incr * *dim;
  for(i=0;i<*dim;i++){
    work[i] = x[i] - m[i+incrVec];
  }
  calcQuadform(work, SigmaInv, dim, incr, &out);
  out = *log2pi - 0.5*log(*det) - 0.5*out;
  return exp(out);
}


void calcmixmvdens(double *x, int *dim, int *ncomp, double *m, double *SigmaInv,
		  double *dets, double *weights, double *work, double *out){
  int i=0;
  double log2pi = 0.0;
  double temp = 0.0;
  log2pi = -0.5 * (double) *dim * log(2.0 * M_PI);
  for(i=0;i<*ncomp;i++){
    temp = mvdens(x, dim, m, SigmaInv, &dets[i], 
		 work, &log2pi, &i);
    *out += weights[i]*temp;
  }
}
