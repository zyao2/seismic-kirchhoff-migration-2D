#include "mex.h"
double *dCoef(int order, const char* type);
double *diffOperator2d(const double *pData, mwSize m, mwSize n, const double *pCoeff, int order, double dist, int dim);
double *dampPml(const double *pu, const double *pv, int m, int n, double L);
void FD_fowardMideling(double dz,double dx,double dt, int nz, int nx, int nt, int diffOrder, int boundary, double *pVelocityModel,
        double *pSource,double *pData, double *pSnapshot);

/* ====================================================================== */
double* dCoef(int order, const char* type)
{
    /* begin of declaration */
    double *pA, *pb, *plhs_mldivide, *pCoeff;
    mxArray *lhs_mldivide[1], *rhs_mldivide[2];
    
    int i, j;
    /* end of declaration */
    
    pA = (double*)mxCalloc(order * order, sizeof(double));
    pb = (double*)mxCalloc(order, sizeof(double));
    
    if (!strcmp(type, "r"))
    {
        pb[0] = 1.0/2.0;
        for (j = 0; j < order; j++)
            for (i = 0; i < order; i++)
                pA[j * order + i] = pow(j+1, 2*(i+1)-1);
    }
    else if (!strcmp(type, "s"))
    {
        pb[0] = 1;
        for (j = 0; j < order; j++)
            for (i = 0; i < order; i++)
                pA[j * order + i] = pow(2*(j+1)-1, 2*(i+1)-1);
    }
    else
    {
        mexErrMsgTxt("Type must be \'s\' or \'r\'!");
    }
    
    /* c = A \ b; */
    /* Point mxArray* to double* */
    rhs_mldivide[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetData(rhs_mldivide[0], pA);
    mxSetM(rhs_mldivide[0], order);
    mxSetN(rhs_mldivide[0], order);
    
    rhs_mldivide[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetData(rhs_mldivide[1], pb);
    mxSetM(rhs_mldivide[1], order);
    mxSetN(rhs_mldivide[1], 1);
    
    /* Call Matlab function to calculate c = A \ b, or equivalently, c = mldivide(A, b) */
    mexCallMATLAB(1, lhs_mldivide, 2, rhs_mldivide, "mldivide");
    
    /* get data for output from mxArray */
    plhs_mldivide = mxGetPr(lhs_mldivide[0]);
    pCoeff = (double*)mxCalloc(order, sizeof(double));
    for (i = 0; i < order; i++)
        pCoeff[i] = plhs_mldivide[i];
    
    /* ATTENTION: Don't forget to free dynamic memory allocated by the mxCreate* function(s) (except for output arrays), otherwise memory leak will occur */
    mxDestroyArray(rhs_mldivide[0]);
    mxDestroyArray(rhs_mldivide[1]);
    mxDestroyArray(lhs_mldivide[0]);
    
    return pCoeff;
}


double* diffOperator2d(const double *pData, mwSize m, mwSize n, const double *pCoeff, int order, double dist, int dim)
{
    /* begin of declaration */
    double *pDiffData;
    
    int l, iOrder, i, j, k, idx1, idx2;
    /* end of declaration */
    
    if (dim < 1)
        dim = 1;
    
    l = 2 * order - 1;
    
    if (dim == 1)
    {
        pDiffData = (double*)mxCalloc((m-l) * n, sizeof(double));
        for (iOrder = 0; iOrder < order; iOrder++)
        {
            for (j = 0; j < n; j++)
            {
                i = 0;
                idx2 = order - (iOrder + 1);
                for (idx1 = order+iOrder; idx1 <= m-order+iOrder; idx1++)
                {
                    pDiffData[j*(m-l) + i] = pDiffData[j*(m-l) + i] +
                            pCoeff[iOrder] * (pData[j*m + idx1] - pData[j*m + idx2]) / dist;
                    i++;
                    idx2++;
                }
            }
        }
    }
    else    /* dim == 2 */
    {
        pDiffData = (double*)mxCalloc(m * (n-l), sizeof(double));
        for (iOrder = 0; iOrder < order; iOrder++)
        {
            for (i = 0; i < m; i++)
            {
                j = 0;
                idx2 = order - (iOrder + 1);
                for (idx1 = order+iOrder; idx1 <= n-order+iOrder; idx1++)
                {
                    pDiffData[j*m + i] = pDiffData[j*m + i] +
                            pCoeff[iOrder] * (pData[idx1*m + i] - pData[idx2*m + i]) / dist;
                    j++;
                    idx2++;
                }
            }
        }
    }
    
    return pDiffData;
}

/* ====================================================================== */
double* dampPml(const double *pu, const double *pv, mwSize m, mwSize n, double L)
{
    /* begin of declaration */
    double *pd0, *pd;
    
    int i, j;
    /* end of declaration */
    
    const double R = 1e-6;
    const double logR = log(R);
    
    /* d0 = -(3 * v)/(2 * L) * log(R); */
    pd0 = (double*)mxCalloc(m * n, sizeof(double));
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            pd0[j * m + i] = -(3.0 * pv[j * m + i])/(2 * L) * logR;
    
    /* d = d0 .* (u / L).^2; */
    pd = (double*)mxCalloc(m * n, sizeof(double));
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            pd[j * m + i] = pd0[j * m + i] * (pu[j * m + i] / L) * (pu[j * m + i] / L);
    mxFree(pd0);
    
    return pd;
}


void FD_fowardMideling(double dz,double dx,double dt, int nz, int nx, int nt, int diffOrder, int boundary, double *pVelocityModel,
        double *pSource,double *pData, double *pSnapshot)
{
    int l, i, j, t;
    double *pCoeff, *pOldFdm, *pCurFdm, *pNewFdm;
    double *puDampLeft, *pvDampLeft, *puDampRight, *pvDampRight, *puDampDown, *pvDampDown;
    double *pxDampLeft, *pxDampRight, *pxDamp, *pxb, *pzDampDown, *pzDamp, *pzb;
    double *pVdtSq;
    double *pzPhi, *pxPhi, *pzA, *pxA, *pzPsi, *pxPsi, *pzP, *pxP;
    double *pCurFdm_diffIn_zPhi, *pCurFdm_diffOut_zPhi, *pCurFdm_diffIn_xPhi, *pCurFdm_diffOut_xPhi;
    double *pCurFdm_diffIn_zA, *pCurFdm_diffOut_zA, *pCurFdm_diffIn_xA, *pCurFdm_diffOut_xA;
    double *pzA_diffIn, *pzA_diffOut, *pxA_diffIn, *pxA_diffOut;
    pCoeff = dCoef(diffOrder, "s");
    l = 2 * diffOrder - 1;
   
   
    puDampLeft = (double*)mxCalloc(nz * boundary, sizeof(double));
    for (j = 0; j < boundary; j++)
        for (i = 0; i < nz; i++)
            puDampLeft[j * nz + i] = (boundary - j) * dx;
    pvDampLeft = (double*)mxCalloc(nz * boundary, sizeof(double));
    memcpy(pvDampLeft, pVelocityModel, sizeof(double) * nz * boundary);
    pxDampLeft = dampPml(puDampLeft, pvDampLeft, nz, boundary, boundary * dx);
   
    puDampRight = (double*)mxCalloc(nz * boundary, sizeof(double));
    for (j = 0; j < boundary; j++)
        for (i = 0; i < nz; i++)
            puDampRight[j * nz + i] = (j + 1) * dx;
    pvDampRight = (double*)mxCalloc(nz * boundary, sizeof(double));
    memcpy(pvDampRight, pVelocityModel + (nx-boundary) * nz, sizeof(double) * nz * boundary);
    pxDampRight = dampPml(puDampRight, pvDampRight, nz, boundary, boundary * dx);
    
    pxDamp = (double*)mxCalloc(nz * nx, sizeof(double));
    memcpy(pxDamp, pxDampLeft, sizeof(double) * nz * boundary);
    memcpy(pxDamp + (nx-boundary) * nz, pxDampRight, sizeof(double) * nz * boundary);
    
    pxb = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pxb[j * nz + i] = exp(-pxDamp[j * nz + i] * dt);
    
    
    puDampDown = (double*)mxCalloc(boundary * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for(i = 0; i < boundary; i++)
            puDampDown[j * boundary + i] = (i + 1) * dz;
    pvDampDown = (double*)mxCalloc(boundary * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for(i = 0; i < boundary; i++)
            pvDampDown[j * boundary + i] = pVelocityModel[j * nz + (nz - boundary + i)];
    pzDampDown = dampPml(puDampDown, pvDampDown, boundary, nx, boundary * dz);
    
    pzDamp = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = nz-boundary; i < nz; i++)
            pzDamp[j * nz + i] = pzDampDown[j * boundary + i-(nz-boundary)];
    
    pzb = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pzb[j * nz + i] = exp(-pzDamp[j * nz + i] * dt);
   
    pOldFdm = (double*)mxCalloc((nz+2*l) * (nx+2*l), sizeof(double));
    pCurFdm = (double*)mxCalloc((nz+2*l) * (nx+2*l), sizeof(double));
    pNewFdm = (double*)mxCalloc((nz+2*l) * (nx+2*l), sizeof(double));
    
    pzPhi = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
    pxPhi = (double*)mxCalloc(nz * (nx+2*l), sizeof(double));
    pzA = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
    pxA = (double*)mxCalloc(nz * (nx+2*l), sizeof(double));
    pzPsi = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pxPsi = (double*)mxCalloc(nz * (nx+l), sizeof(double));
    pzP = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pxP = (double*)mxCalloc(nz * (nx+l), sizeof(double));
    
    pVdtSq = (double*)mxCalloc(nz * nx, sizeof(double));
    for (j = 0; j < nx; j++)
        for (i = 0; i < nz; i++)
            pVdtSq[j * nz + i] = (pVelocityModel[j * nz + i] * dt) * (pVelocityModel[j * nz + i] * dt);
    
    pCurFdm_diffIn_zPhi = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pCurFdm_diffIn_xPhi = (double*)mxCalloc(nz * (nx+l), sizeof(double));
    pCurFdm_diffIn_zA = (double*)mxCalloc((nz+2*l) * nx, sizeof(double));
    pCurFdm_diffIn_xA = (double*)mxCalloc(nz * (nx+2*l), sizeof(double));
    pzA_diffIn = (double*)mxCalloc((nz+l) * nx, sizeof(double));
    pxA_diffIn = (double*)mxCalloc(nz * (nx+l), sizeof(double));

    for (t = 0; t < nt; t++)
    {
       
        for (j = l; j < nx+l; j++)
            for (i = diffOrder; i < nz+2*l-diffOrder+1; i++)
                pCurFdm_diffIn_zPhi[(j - l) * (nz+l) + (i-diffOrder)] = pCurFdm[j * (nz+2*l) + i];
        pCurFdm_diffOut_zPhi = diffOperator2d(pCurFdm_diffIn_zPhi, nz+l, nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzPhi[j * (nz+2*l) + i] = pzb[j * nz + (i - l)] * pzPhi[j * (nz+2*l) + i] +
                        (pzb[j * nz + (i - l)] - 1) * pCurFdm_diffOut_zPhi[j * nz + (i - l)];
        
        
        for (j = diffOrder; j < nx+2*l-diffOrder+1; j++)
            for (i = l; i < nz+l; i++)
                pCurFdm_diffIn_xPhi[(j-diffOrder) * nz + (i - l)] = pCurFdm[j * (nz+2*l) + i];
        pCurFdm_diffOut_xPhi = diffOperator2d(pCurFdm_diffIn_xPhi, nz, nx+l, pCoeff, diffOrder, dx, 2);
        
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPhi[j * nz + i] = pxb[(j - l) * nz + i] * pxPhi[j * nz + i] +
                        (pxb[(j - l) * nz + i] - 1) * pCurFdm_diffOut_xPhi[(j - l) * nz + i];
        
        
        memcpy(pCurFdm_diffIn_zA, pCurFdm + l * (nz+2*l), sizeof(double) * nx * (nz+2*l));
        pCurFdm_diffOut_zA = diffOperator2d(pCurFdm_diffIn_zA, nz+2*l, nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA[j * (nz+2*l) + i] = pCurFdm_diffOut_zA[j * (nz+l) + (i - (diffOrder - 1))] + pzPhi[j * (nz+2*l) + i];
        
        
        for (j = 0; j < nx+2*l; j++)
            for (i = l; i < nz+l; i++)
                pCurFdm_diffIn_xA[j * nz + (i - l)] = pCurFdm[j * (nz+2*l) + i];
        pCurFdm_diffOut_xA = diffOperator2d(pCurFdm_diffIn_xA, nz, nx+2*l, pCoeff, diffOrder, dx, 2);
        
        for (j = diffOrder - 1; j < nx+2*l-diffOrder; j++)
            for (i = 0; i < nz; i++)
                pxA[j * nz + i] = pCurFdm_diffOut_xA[(j - (diffOrder - 1)) * nz + i] + pxPhi[j * nz + i];
        
       
        for (j = 0; j < nx; j++)
            for (i = diffOrder - 1; i < nz+2*l-diffOrder; i++)
                pzA_diffIn[j * (nz+l) + (i - (diffOrder - 1))] = pzA[j * (nz+2*l) + i];
        pzA_diffOut = diffOperator2d(pzA_diffIn, nz+l, nx, pCoeff, diffOrder, dz, 1);
        
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzPsi[j * (nz+l) + i] = pzb[j * nz + (i - l)] * pzPsi[j * (nz+l) + i] +
                        (pzb[j * nz + (i - l)] - 1) * pzA_diffOut[j * nz + (i - l)];
        
       
        memcpy(pxA_diffIn, pxA + (diffOrder - 1) * nz, sizeof(double) * (nx+l) * nz);
        pxA_diffOut = diffOperator2d(pxA_diffIn, nz, nx+l, pCoeff, diffOrder, dx, 2);
        
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxPsi[j * nz + i] = pxb[(j - l) * nz + i] * pxPsi[j * nz + i] +
                        (pxb[(j - l) * nz + i] - 1) * pxA_diffOut[(j - l) * nz + i];
        
        
        for (j = 0; j < nx; j++)
            for (i = l; i < nz + l; i++)
                pzP[j * (nz+l) + i] = pzA_diffOut[j * nz + (i - l)] + pzPsi[j * (nz+l) + i];
        
       
        for (j = l; j < nx + l; j++)
            for (i = 0; i < nz; i++)
                pxP[j * nz + i] = pxA_diffOut[(j - l) * nz + i] + pxPsi[j * nz + i];
 
        for (j = l; j < nx + l; j++)
            for (i = l; i < nz + l; i++)
                pNewFdm[j * (nz+2*l) + i] = pVdtSq[(j - l) * nz + (i - l)] *
                        ( pzP[(j - l) * (nz+l) + i] + pxP[j * nz + (i - l)] + pSource[t * (nz * nx) + (j - l) * nz + (i - l)] ) +
                        2 * pCurFdm[j * (nz+2*l) + i] - pOldFdm[j * (nz+2*l) + i];
        
 
        memcpy(pOldFdm, pCurFdm, sizeof(double) * (nz+2*l) * (nx+2*l));
        memcpy(pCurFdm, pNewFdm, sizeof(double) * (nz+2*l) * (nx+2*l));

        for (i = 0; i < nx; i++)
            pData[t * nx + i] = pCurFdm[(i + l) * (nz+2*l) + l];

        for (j = 0; j < nx; j++)
            for (i = 0; i < nz; i++)
                pSnapshot[t * (nz * nx) + j * nz + i] = pCurFdm[(j + l) * (nz+2*l) + (i + l)];
        
   
        mxFree(pCurFdm_diffOut_zPhi);
        mxFree(pCurFdm_diffOut_xPhi);
        mxFree(pCurFdm_diffOut_zA);
        mxFree(pCurFdm_diffOut_xA);
        mxFree(pzA_diffOut);
        mxFree(pxA_diffOut);
    }
  
    mxFree(pCoeff);
    mxFree(pOldFdm);
    mxFree(pCurFdm);
    mxFree(pNewFdm);
    mxFree(puDampLeft);
    mxFree(pvDampLeft);
    mxFree(puDampRight);
    mxFree(pvDampRight);
    mxFree(puDampDown);
    mxFree(pvDampDown);
    mxFree(pxDampLeft);
    mxFree(pxDampRight);
    mxFree(pxDamp);
    mxFree(pxb);
    mxFree(pzDampDown);
    mxFree(pzDamp);
    mxFree(pzb);
    mxFree(pVdtSq);
    mxFree(pzPhi);
    mxFree(pxPhi);
    mxFree(pzA);
    mxFree(pxA);
    mxFree(pzPsi);
    mxFree(pxPsi);
    mxFree(pzP);
    mxFree(pxP);
    mxFree(pCurFdm_diffIn_zPhi);
    mxFree(pCurFdm_diffIn_xPhi);
    mxFree(pCurFdm_diffIn_zA);
    mxFree(pCurFdm_diffIn_xA);
    mxFree(pzA_diffIn);
    mxFree(pxA_diffIn);  
}
