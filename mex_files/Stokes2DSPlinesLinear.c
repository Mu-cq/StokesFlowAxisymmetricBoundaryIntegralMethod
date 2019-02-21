//
//  Stokes2DSPlinesLinear.c
//  
//
//  Created by Giacomo on 12/24/15.
//
//

#include "mex.h"
#include "math.h"
#include "stdio.h"

void Stokes2DSPlinesLinear(double *aax, double *bbx, double *ccx, double *ddx, double *aay, double *bby, double *ccy, double *ddy, double *x0, double *y0, double *GXX, double *GXY, double *GYY, double *TXXX, double *TXXY, double *TXYY, double *TYYY, int field, int elem)
{
  int i, j, dolog1, dolog2, dolog;
  double SXX, SXY, SYY;    //single layer
  double QXXX, QXXY, QXYY, QYYY;    //double layer
  double ax, bx, cx, dx, ay, by, cy, dy, h, h0, h1;
  double r, r0, d, d2, d3, d5, GPpar, phiA, phiB, eta, beta, deta, dbeta;
  double DX, DY, nx, ny;
  double distXA, distYA, distA, distXB, distYB, distB, Xleft, Yleft, Xright, Yright;
  double *GXXa, *GXYa, *GYXa, *GYYa, *TXXXa, *TXYXa, *TXXYa, *TXYYa, *TYXXa, *TYXYa, *TYYXa, *TYYYa;
  double *GXXb, *GXYb, *GYXb, *GYYb, *TXXXb, *TXYXb, *TXXYb, *TXYYb, *TYXXb, *TYXYb, *TYYXb, *TYYYb;
  
  //memory allocation
  GXXa = (double*)calloc(field*elem,sizeof(double));
  GXYa = (double*)calloc(field*elem,sizeof(double));
  GYYa = (double*)calloc(field*elem,sizeof(double));
  GXXb = (double*)calloc(field*elem,sizeof(double));
  GXYb = (double*)calloc(field*elem,sizeof(double));
  GYYb = (double*)calloc(field*elem,sizeof(double));
  TXXXa = (double*)calloc(field*elem,sizeof(double));
  TXXYa = (double*)calloc(field*elem,sizeof(double));
  TXYYa = (double*)calloc(field*elem,sizeof(double));
  TYYYa = (double*)calloc(field*elem,sizeof(double));
  TXXXb = (double*)calloc(field*elem,sizeof(double));
  TXXYb = (double*)calloc(field*elem,sizeof(double));
  TXYYb = (double*)calloc(field*elem,sizeof(double));
  TYYYb = (double*)calloc(field*elem,sizeof(double));
  
  //FILE *pFile1;
  //FILE *pFile2;
  //pFile1 = fopen ("myfile1.txt","w");
  //pFile2 = fopen ("myfile2.txt","w");
  
  //define pi
  double pi = 3.141592653589793;
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
  
  for (j=0; j<elem; j++) {  
      //elements variables      
      ax = aax[j];  bx = bbx[j];    cx = ccx[j];    dx = ddx[j];
      ay = aay[j];  by = bby[j];    cy = ccy[j];    dy = ddy[j];
      
      //metrics for first and last term
      h0 = sqrt(bx*bx+by*by);
      h1 =  sqrt((bx+2*cx+3*dx)*(bx+2*cx+3*dx)+(by+2*cy+3*dy)*(by+2*cy+3*dy));
      
      //first and last point of the element
      Xleft = ax;   Yleft = ay;     Xright = ax+bx+cx+dx;   Yright = ay+by+cy+dy;
      
      for (i=0; i<field; i++) {
          
          //decide if doing singular treatment
          //distXA = Xleft-x0[i];
          //distYA = Yleft-y0[i];
          //distA = sqrt(distXA*distXA+distYA*distYA);
          //dolog1 = (distA>10e-8);
          dolog1 = 1;
          
          //decide if doing singular treatment
          //distXB = Xright-x0[i];
          //distYB = Yright-y0[i];
          //distB = sqrt(distXB*distXB+distYB*distYB);
          //dolog2 = (distB>10e-8);
          dolog2 = 1;
          
          dolog = (dolog1+dolog2 > 1.9);
          
          //gauss integration
          for (int l=0; l<6; l++) {
                //gauss point in parametric space
                GPpar = (GP[l]+1)/2;
                                
                //gauss points in physical space
                eta = ax+bx*GPpar+cx*GPpar*GPpar+dx*GPpar*GPpar*GPpar;
                beta = ay+by*GPpar+cy*GPpar*GPpar+dy*GPpar*GPpar*GPpar;
                deta = bx+2*cx*GPpar+3*dx*GPpar*GPpar;
                dbeta = by+2*cy*GPpar+3*dy*GPpar*GPpar;
                
                //normal to the spline
                nx = dbeta/sqrt(deta*deta+dbeta*dbeta);
                ny = -deta/sqrt(deta*deta+dbeta*dbeta);
                
                //fprintf (pFile1, "%f \n",eta);
                //fprintf (pFile2, "%f \n",nx);
                
                //metric
                h = sqrt(deta*deta+dbeta*dbeta);
                
                //compute green function physical variables
                DX = eta-x0[i];
                DY = beta-y0[i];
                d2 = DX*DX+DY*DY;
                d = sqrt(d2);
                
                //hat function
                phiA = 1.0-(GP[l]+1)/2;
                phiB = (GP[l]+1)/2;
                
                //single layer
                SXX = -log(d)*dolog+(DX*DX)/pow(d,2);
                SXX += (-log(d/GPpar) - log(GPpar)*(1-h0/h))*(1-dolog1);     //remove singularity on the left
                SXX += (-log(d/(1-GPpar)) - log(1-GPpar)*(1-h1/h))*(1-dolog2);     //remove singularity on the right
                
                SXY = (DX*DY)/pow(d,2);
                
                SYY = -log(d)*dolog+(DY*DY)/pow(d,2);
                SYY += (-log(d/GPpar) - log(GPpar)*(1-h0/h))*(1-dolog1);     //remove singularity on the left
                SYY += (-log(d/(1-GPpar)) - log(1-GPpar)*(1-h1/h))*(1-dolog2);     //remove singularity on the right
                                
                //double layer
                QXXX = -4*(DX*DX*DX)/pow(d,4)*dolog;
                QXXY = -4*(DX*DX*DY)/pow(d,4)*dolog;
                QXYY = -4*(DX*DY*DY)/pow(d,4)*dolog;
                QYYY = -4*(DY*DY*DY)/pow(d,4)*dolog;
              
                //integration phiA SINGLE LAYER
                GXXa[j*field+i] += phiA*GW[l]*SXX/2 * h;
                GXYa[j*field+i] += phiA*GW[l]*SXY/2 * h;
                GYYa[j*field+i] += phiA*GW[l]*SYY/2 * h;
                
                //integration phiA DOUBLE LAYER
                TXXXa[j*field+i] += phiA*GW[l]*QXXX/2*nx * h;
                TXXYa[j*field+i] += phiA*GW[l]*QXXY/2*ny * h;
                TXYYa[j*field+i] += phiA*GW[l]*QXYY/2*ny * h;
                TYYYa[j*field+i] += phiA*GW[l]*QYYY/2*ny * h;
                
                //integration phiB SINGLE LAYER
                GXXb[j*field+i] += phiB*GW[l]*SXX/2 * h;
                GXYb[j*field+i] += phiB*GW[l]*SXY/2 * h;
                GYYb[j*field+i] += phiB*GW[l]*SYY/2 * h;
                
                //integration phiB DOUBLE LAYER
                TXXXb[j*field+i] += phiB*GW[l]*QXXX/2*nx * h;
                TXXYb[j*field+i] += phiB*GW[l]*QXXY/2*ny * h;
                TXYYb[j*field+i] += phiB*GW[l]*QXYY/2*ny * h;
                TYYYb[j*field+i] += phiB*GW[l]*QYYY/2*ny * h;
    
          }
                    
          //singular treatment: add part form analytical integration (as in research notes) IF I'M ON THE AXIS I DON'T NEED SINGULAR TREATMENT
          GXXa[j*field+i] += 0.75*h0*(1-dolog1);
          GYYa[j*field+i] += 0.75*h0*(1-dolog1);
          
          GXXb[j*field+i] += 0.25*h0*(1-dolog1);
          GYYb[j*field+i] += 0.25*h0*(1-dolog1);
          
          GXXa[j*field+i] += 0.25*h1*(1-dolog2);
          GYYa[j*field+i] += 0.25*h1*(1-dolog2);
          
          GXXb[j*field+i] += 0.75*h1*(1-dolog2);
          GYYb[j*field+i] += 0.75*h1*(1-dolog2);
           
          
          //build output matrix
          GXX[j*field+i] += GXXa[j*field+i];
          GXY[j*field+i] += GXYa[j*field+i];
          GYY[j*field+i] += GYYa[j*field+i];
          TXXX[j*field+i] += TXXXa[j*field+i];          
          TXXY[j*field+i] += TXXYa[j*field+i];         
          TXYY[j*field+i] += TXYYa[j*field+i];
          TYYY[j*field+i] += TYYYa[j*field+i];
          
          //exception if it is last element
          if (j==elem-1) {
                GXX[i] += GXXb[j*field+i];
                GXY[i] += GXYb[j*field+i];
                GYY[i] += GYYb[j*field+i];
                TXXX[i] += TXXXb[j*field+i];
                TXXY[i] += TXXYb[j*field+i];
                TXYY[i] += TXYYb[j*field+i];
                TYYY[i] += TYYYb[j*field+i];
          }else{
                GXX[(j+1)*field+i] += GXXb[j*field+i];
                GXY[(j+1)*field+i] += GXYb[j*field+i];
                GYY[(j+1)*field+i] += GYYb[j*field+i];
                TXXX[(j+1)*field+i] += TXXXb[j*field+i];
                TXXY[(j+1)*field+i] += TXXYb[j*field+i];
                TXYY[(j+1)*field+i] += TXYY[j*field+i];
                TYYY[(j+1)*field+i] += TYYYb[j*field+i];
          }
          
      }
  }
  
  
  
  //free memory
  free(GXXa);
  free(GXYa);
  free(GYYa);
  free(TXXXa);
  free(TXXYa);
  free(TXYYa);
  free(TYYYa);
  free(GXXb);
  free(GXYb);
  free(GYYb);
  free(TXXXb);
  free(TXXYb);
  free(TXYYb);
  free(TYYYb);
  
  //fclose (pFile1);
  //fclose (pFile2);
  
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
double *ax;       /* 1xN input matrix */
double *bx;       /* 1xN input matrix */
double *cx;       /* 1xN input matrix */
double *dx;       /* 1xN input matrix */
double *ay;       /* 1xN input matrix */
double *by;       /* 1xN input matrix */
double *cy;       /* 1xN input matrix */
double *dy;       /* 1xN input matrix */
mwSize nrows;           /* size of matrix */
double *inMatrixX0;       /* 1xM input matrix */
double *inMatrixY0;       /* 1xM input matrix */
mwSize ncols;           /* size of matrix */
double *GXX;      /* output matrix 1x(NxM) */
double *GXY;      /* output matrix 1x(NxM) */
double *GYY;      /* output matrix 1x(NxM) */
double *TXXX;      /* output matrix 1x(NxM) */
double *TXXY;      /* output matrix 1x(NxM) */
double *TXYY;      /* output matrix 1x(NxM) */
double *TYYY;      /* output matrix 1x(NxM) */


/* get the value of the scalar input  */
//multiplier = mxGetScalar(prhs[0]);

/* code here */
/* create a pointer to the real data in the input matrix  */
ax = mxGetPr(prhs[0]);
bx = mxGetPr(prhs[1]);
cx = mxGetPr(prhs[2]);
dx = mxGetPr(prhs[3]);
ay = mxGetPr(prhs[4]);
by = mxGetPr(prhs[5]);
cy = mxGetPr(prhs[6]);
dy = mxGetPr(prhs[7]);
inMatrixX0 = mxGetPr(prhs[8]);
inMatrixY0 = mxGetPr(prhs[9]);
    
/* get dimensions of the input matrix */
nrows = mxGetN(prhs[0]);
ncols = mxGetN(prhs[8]);
    
/* create the output matrix (their dimension depending on the inputs) */
plhs[0] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[4] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[5] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[6] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);


/* get a pointer to the real data in the output matrix */
GXX = mxGetPr(plhs[0]);
GXY = mxGetPr(plhs[1]);
GYY = mxGetPr(plhs[2]);
TXXX = mxGetPr(plhs[3]);
TXXY = mxGetPr(plhs[4]);
TXYY = mxGetPr(plhs[5]);
TYYY = mxGetPr(plhs[6]);

/* call the computational routine */
Stokes2DSPlinesLinear(ax, bx, cx, dx, ay, by, cy, dy, inMatrixX0, inMatrixY0, GXX, GXY, GYY, TXXX, TXXY, TXYY, TYYY, ncols, nrows);
}