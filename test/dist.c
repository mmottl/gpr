/* Compute matrix of pairwise distances between elements */
/* of two vectors */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) 
{
  int n0, n1;
  double *x0, *x1, *D;
  int i, j, count=0;

  if (nrhs != 2)
    mexErrMsgTxt("Two input vectors required.");
  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  if (mxGetN(prhs[0]) != 1 || mxGetN(prhs[1]) != 1)
    mexErrMsgTxt("Both inputs must be column vectors.");

  /* Get the length of each input vector. */
  n0 = mxGetM(prhs[0]);
  n1 = mxGetM(prhs[1]);

  x0 = mxGetPr(prhs[0]);
  x1 = mxGetPr(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(n0, n1, mxREAL);
  
  D = mxGetPr(plhs[0]);

  for (j = 0; j < n1; j++) {
    for (i = 0; i < n0; i++) {
      *(D+count) = x0[i] - x1[j];
      count++;
    }
  }
}
