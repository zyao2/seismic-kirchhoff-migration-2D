#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* begin of declaration */
    double *pVelocityModel, *pSource, *pData, *pSnapshot;
    double dz, dx, dt;
    int diffOrder, boundary;
    int l, i, j, t;
    mwSize nz, nx, nt;
    const mwSize *pDimsSource;
    mwSize pDimsSnapshot[3] = {0};

    pVelocityModel = mxGetPr(prhs[0]);
    pSource = mxGetPr(prhs[1]);
    diffOrder = *mxGetPr(prhs[2]);
    boundary = *mxGetPr(prhs[3]);
    dz = *mxGetPr(prhs[4]);
    dx = *mxGetPr(prhs[5]);
    dt = *mxGetPr(prhs[6]);
    
    pDimsSource = mxGetDimensions(prhs[1]);
    nz = pDimsSource[0];
    nx = pDimsSource[1];
    nt = pDimsSource[2];
  
    /* initialize storage */
    plhs[0] = mxCreateDoubleMatrix(nx, nt, mxREAL);
    pData = mxGetPr(plhs[0]);
    
    pDimsSnapshot[0] = nz;
    pDimsSnapshot[1] = nx;
    pDimsSnapshot[2] = nt;
    plhs[1] = mxCreateNumericArray(3, pDimsSnapshot, mxDOUBLE_CLASS, mxREAL);
    pSnapshot = mxGetPr(plhs[1]);
    
    FD_fowardMideling(dz, dx, dt, nz,nx, nt, diffOrder, boundary, pVelocityModel,pSource,
        pData, pSnapshot);
  
}