//
//  Stokes2DAxisIntLinearProva.c
//  
//
//  Created by Giacomo on 7/1/15.
//
//

#include "mex.h"
#include "math.h"
#include "stdio.h"

void Stokes2DAxisIntLinearProva(double *x, double *y, double *x0, double *y0, double *GXX, double *GXY, double *GYX, double *GYY, double *TXXX, double *TXXY, double *TXYX, double *TXYY, double *TYXX, double *TYXY, double *TYYX, double *TYYY, double *T1, double *T2, double *D1, double *D2, int field, int elem)
{
  int i, j, dolog1, dolog2, dolog, axis;
  double SXX, SXY, SYX, SYY;    //single layer
  double QXXX, QXXY, QXYY, QYXX, QYXY, QYYY;    //double layer
  double PXX, PXY, PYX, PYY;    //desingularization double layer
  double F, E, G, B, D, C, k, k2;   //elliptic integrals
  double I10, I11, I30, I31, I32, I50, I51, I52, I53;
  double r, r0, d, d2, d3, d5, GPx, GPy, phiA, phiB;
  double k2p, k4, k5, k6, yy5;
  double FCTR, RL30, RL50, RL52, RL54, RL56;
  double DX, DY, DX2, DX3, r2, r3, r02, r03, nx, ny;
  double dx, dy, dl, xMid, yMid, distXA, distYA, distA, distXB, distYB, distB;
  double *GXXa, *GXYa, *GYXa, *GYYa, *TXXXa, *TXYXa, *TXXYa, *TXYYa, *TYXXa, *TYXYa, *TYYXa, *TYYYa;
  double *GXXb, *GXYb, *GYXb, *GYYb, *TXXXb, *TXYXb, *TXXYb, *TXYYb, *TYXXb, *TYXYb, *TYYXb, *TYYYb;
  
  //memory allocation
  GXXa = (double*)calloc(field*elem,sizeof(double));
  GXYa = (double*)calloc(field*elem,sizeof(double));
  GYXa = (double*)calloc(field*elem,sizeof(double));
  GYYa = (double*)calloc(field*elem,sizeof(double));
  GXXb = (double*)calloc(field*elem,sizeof(double));
  GXYb = (double*)calloc(field*elem,sizeof(double));
  GYXb = (double*)calloc(field*elem,sizeof(double));
  GYYb = (double*)calloc(field*elem,sizeof(double));
  TXXXa = (double*)calloc(field*elem,sizeof(double));
  TXXYa = (double*)calloc(field*elem,sizeof(double));
  TXYXa = (double*)calloc(field*elem,sizeof(double));
  TXYYa = (double*)calloc(field*elem,sizeof(double));
  TYXXa = (double*)calloc(field*elem,sizeof(double));
  TYXYa = (double*)calloc(field*elem,sizeof(double));
  TYYXa = (double*)calloc(field*elem,sizeof(double));
  TYYYa = (double*)calloc(field*elem,sizeof(double));
  TXXXb = (double*)calloc(field*elem,sizeof(double));
  TXXYb = (double*)calloc(field*elem,sizeof(double));
  TXYXb = (double*)calloc(field*elem,sizeof(double));
  TXYYb = (double*)calloc(field*elem,sizeof(double));
  TYXXb = (double*)calloc(field*elem,sizeof(double));
  TYXYb = (double*)calloc(field*elem,sizeof(double));
  TYYXb = (double*)calloc(field*elem,sizeof(double));
  TYYYb = (double*)calloc(field*elem,sizeof(double));
  
  //define pi and tol
  //double pi = 3.141592653589793;
  double pi = M_PI;
  double tol = 0.0000000000001;
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
  
  for (j=0; j<elem; j++) {  
      //elemnts variables
      dx = x[j+1]-x[j];
      dy = y[j+1]-y[j];
      dl = sqrt(dx*dx+dy*dy);
      xMid = (x[j+1]+x[j])/2;
      yMid = (y[j+1]+y[j])/2;
      
      //compute normal to the element
      nx = dy/dl;   ny = -dx/dl;
      
      for (i=0; i<field; i++) {
          
          //check if I'm on the axis
          axis = (y0[i] < 10e-8);
          
          //decide if doing singular treatment
          distXA = x[j]-x0[i];
          distYA = y[j]-y0[i];
          distA = sqrt(distXA*distXA+distYA*distYA);
          dolog1 = (distA>10e-8);
          
          //decide if doing singular treatment
          distXB = x[j+1]-x0[i];
          distYB = y[j+1]-y0[i];
          distB = sqrt(distXB*distXB+distYB*distYB);
          dolog2 = (distB>10e-8);
          
          dolog = (dolog1+dolog2 > 1.9);
          
          //gauss integration
          for (int l=0; l<6; l++) {
                //gauss point
                GPx = GP[l]/2*dx+xMid;
                GPy = GP[l]/2*dy+yMid;
                
                //compute green function physical variables
                DX = GPx-x0[i];
                DY = GPy-y0[i];
                d2 = DX*DX+DY*DY;
                d = sqrt(d2);
                r0 = y0[i];
                r = GPy;
                k2 = (4*r0*r)/(DX*DX+pow((r+r0),2));
                k = sqrt(k2);
                r2 = r*r;
                r3 = r*r*r;
                r02 = r0*r0;
                r03 = r0*r0*r0;
                
                //dummy variables
                k2p = 1-k2;
                k4 = k2*k2;
                k6 = k2*k2*k2;
                k5 = k4*k;
                yy5 = pow(r0*r,2.5);
                DX2 = DX*DX;
                DX3 = DX*DX*DX;
                
                //compute elliptic integral of first and second kind
                F = 0.5*pi;
                E = 1.0;
                G = 1.0;
                B = k;
                D = 10*tol;

                while (fabs(D) > tol) {
                  C = sqrt(1.0-B*B);
                  B = (1.0-C)/(1.0+C);
                  D = F*B;
                  F = F+D;
                  G = 0.50*G*B;
                  E = E+G;
                }

                E = F*(1.0-0.50*k2*E);
                
                //often used variables for single layer
                I10 = 2*k/sqrt(r0*r)*F;
                I11 = k/pow(r0*r,1.5)*((DX2+r0*r0+r2)*F-(DX2+r0*r0+r2+2*r*r0)*E);
                I30 = 2*k/pow(r0*r,0.5)*E/d2;
                I31 = k/pow(r0*r,1.5)*(-F+(r0*r0+r2+DX2)/d2*E);
                I32 = k/pow(r0*r,2.5)*(-F*(DX2+r0*r0+r2) + (DX2*DX2+2.0*DX2*r2+2.0*DX2*r0*r0+r2*r2+r0*r0*r0*r0)*E/d2);
                
                //dummy variables for double layer
                RL30 = E/k2p;
                RL50 =  (2.0 * (2.0 - k2)*RL30 - F) / (3.0*k2p);
                RL52 =  (RL50 - RL30)/k2;
                RL54 =  (F-RL50+2.0*k2*RL52)/k4;
                RL56 = -(E-RL50+3.0*k2*RL52-3.0*k4*RL54)/ k6;
                
                //often used variables for double layer
                FCTR = k5/(8.0*yy5);
                I50 = FCTR * RL50;
                I51 = FCTR * (2.0 * RL52-RL50);
                I52 = FCTR * (4.0 * (RL54-RL52)+RL50);
                I53 = FCTR * (8.0 * RL56 - 12.0 *RL54 +6.0 * RL52 - RL50);
                
                if (axis) {
                    
                     d3  = d2*d;
                     d5  = d3*d2;
                     I10 = 2*pi/d;
                     I11 = 0.0;
                     I30 = 2*pi/d3;
                     I31 = 0.0;
                     I32 = pi/d3;

                     I50 = 2*pi/d5;
                     I51 = 0.0;
                     I52 = pi/d5;
                     I53 = 0.0;
                }
                    
                
                //hat function
                phiA = 1.0-(GP[l]+1)/2;
                phiB = (GP[l]+1)/2;
                
                //single layer
                SXX = r*(I10+DX*DX*I30) + 2*log(d)*(1-dolog)*(1-axis);
                SXY = r*DX*(r*I30-r0*I31);
                SYX = r*DX*(r*I31-r0*I30);
                SYY = r*(I11+(r*r+r0*r0)*I31-r*r0*(I30+I32)) + 2*log(d)*(1-dolog)*(1-axis);
                
                //double layer
                QXXX = -6*r*DX3*I50;
                QXXY = -6*r*DX2*(r*I50-r0*I51);
                QXYY = -6*r*DX*(r0*r0*I52+r*r*I50-2*r*r0*I51);
                QYXX = -6*r*DX2*(r*I51-r0*I50);
                QYXY = -6*r*DX*((r*r+r0*r0)*I51-r*r0*(I50+I52));
                QYYY = -6*r*(r3*I51-r2*r0*(I50+2*I52)+r*r0*r0*(I53+2*I51)-r0*r0*r0*I52);
                
                //desingularization double layer
                PXX = QYXX;
                PXY = QYXY;
                PYX = - 6*r*DX * (r2*I52 - 2*r*r0*I51 + r02*I50);
                PYY = - 6*r*(r02*r * (2.0 * I52+I50) - r03 * I51 - r2*r0 * (2.0*I51+I53) + r3 * I52);
              
                //integration phiA SINGLE LAYER
                GXXa[j*field+i] += phiA*GW[l]*dl*SXX/2;
                GXYa[j*field+i] += phiA*GW[l]*dl*SXY/2;
                GYXa[j*field+i] += phiA*GW[l]*dl*SYX/2;
                GYYa[j*field+i] += phiA*GW[l]*dl*SYY/2;
                
                //integration phiA DOUBLE LAYER
                TXXXa[j*field+i] += phiA*GW[l]*dl*QXXX/2*nx;
                TXXYa[j*field+i] += phiA*GW[l]*dl*QXXY/2*ny;
                TXYXa[j*field+i] += phiA*GW[l]*dl*QXXY/2*nx;
                TXYYa[j*field+i] += phiA*GW[l]*dl*QXYY/2*ny;
                TYXXa[j*field+i] += phiA*GW[l]*dl*QYXX/2*nx;
                TYXYa[j*field+i] += phiA*GW[l]*dl*QYXY/2*ny;
                TYYXa[j*field+i] += phiA*GW[l]*dl*QYXY/2*nx;
                TYYYa[j*field+i] += phiA*GW[l]*dl*QYYY/2*ny;
                
                //integration phiB SINGLE LAYER
                GXXb[j*field+i] += phiB*GW[l]*dl*SXX/2;
                GXYb[j*field+i] += phiB*GW[l]*dl*SXY/2;
                GYXb[j*field+i] += phiB*GW[l]*dl*SYX/2;
                GYYb[j*field+i] += phiB*GW[l]*dl*SYY/2;
                
                //integration phiB DOUBLE LAYER
                TXXXb[j*field+i] += phiB*GW[l]*dl*QXXX/2*nx;
                TXXYb[j*field+i] += phiB*GW[l]*dl*QXXY/2*ny;
                TXYXb[j*field+i] += phiB*GW[l]*dl*QXXY/2*nx;
                TXYYb[j*field+i] += phiB*GW[l]*dl*QXYY/2*ny;
                TYXXb[j*field+i] += phiB*GW[l]*dl*QYXX/2*nx;
                TYXYb[j*field+i] += phiB*GW[l]*dl*QYXY/2*ny;
                TYYXb[j*field+i] += phiB*GW[l]*dl*QYXY/2*nx;
                TYYYb[j*field+i] += phiB*GW[l]*dl*QYYY/2*ny;
                
                //integration for desigularization
                T1[j*field+i] += (GW[l]*dl*QXXX/2*nx + GW[l]*dl*QXXY/2*ny);
                T2[j*field+i] += (GW[l]*dl*QYXX/2*nx + GW[l]*dl*QYXY/2*ny);
                
                D1[j*field+i] += (GW[l]*dl*PXX/2*nx + GW[l]*dl*PXY/2*ny);
                D2[j*field+i] += (GW[l]*dl*PYX/2*nx + GW[l]*dl*PYY/2*ny);
    
          }
                    
          //singular treatment: add part form analytical integration (as in research notes) IF I'M ON THE AXIS I DON'T NEED SINGULAR TREATMENT
          GXXa[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          GYYa[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          
          GXXb[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          GYYb[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          
          GXXa[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          GYYa[j*field+i] += (0.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          
          GXXb[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          GYYb[j*field+i] += (1.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          
          //build output matrix
          GXX[j*field+i] += GXXa[j*field+i];
          GXX[(j+1)*field+i] += GXXb[j*field+i];
          
          GXY[j*field+i] += GXYa[j*field+i];
          GXY[(j+1)*field+i] += GXYb[j*field+i];
          
          GYX[j*field+i] += GYXa[j*field+i];
          GYX[(j+1)*field+i] += GYXb[j*field+i];
          
          GYY[j*field+i] += GYYa[j*field+i];
          GYY[(j+1)*field+i] += GYYb[j*field+i];
          
          TXXX[j*field+i] += TXXXa[j*field+i];
          TXXX[(j+1)*field+i] += TXXXb[j*field+i];
          
          TXXY[j*field+i] += TXXYa[j*field+i];
          TXXY[(j+1)*field+i] += TXXYb[j*field+i];
          
          TXYX[j*field+i] += TXYXa[j*field+i];
          TXYX[(j+1)*field+i] += TXYXb[j*field+i];
          
          TXYY[j*field+i] += TXYYa[j*field+i];
          TXYY[(j+1)*field+i] += TXYYb[j*field+i];
          
          TYXX[j*field+i] += TYXXa[j*field+i];
          TYXX[(j+1)*field+i] += TYXXb[j*field+i];
          
          TYXY[j*field+i] += TYXYa[j*field+i];
          TYXY[(j+1)*field+i] += TYXYb[j*field+i];
          
          TYYX[j*field+i] += TYYXa[j*field+i];
          TYYX[(j+1)*field+i] += TYYXb[j*field+i];
          
          TYYY[j*field+i] += TYYYa[j*field+i];
          TYYY[(j+1)*field+i] += TYYYb[j*field+i];
          
      }
  }
  
  //free memory
  free(GXXa);
  free(GXYa);
  free(GYXa);
  free(GYYa);
  free(TXXXa);
  free(TXXYa);
  free(TXYXa);
  free(TXYYa);
  free(TYXXa);
  free(TYXYa);
  free(TYYXa);
  free(TYYYa);
  free(GXXb);
  free(GXYb);
  free(GYXb);
  free(GYYb);
  free(TXXXb);
  free(TXXYb);
  free(TXYXb);
  free(TXYYb);
  free(TYXXb);
  free(TYXYb);
  free(TYYXb);
  free(TYYYb);
  
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
double *inMatrixX;       /* 1xN input matrix */
double *inMatrixY;       /* 1xN input matrix */
mwSize nrows;           /* size of matrix */
double *inMatrixX0;       /* 1xM input matrix */
double *inMatrixY0;       /* 1xM input matrix */
mwSize ncols;           /* size of matrix */
double *GXX;      /* output matrix 1x(NxM) */
double *GXY;      /* output matrix 1x(NxM) */
double *GYX;      /* output matrix 1x(NxM) */
double *GYY;      /* output matrix 1x(NxM) */
double *TXXX;      /* output matrix 1x(NxM) */
double *TXXY;      /* output matrix 1x(NxM) */
double *TXYX;      /* output matrix 1x(NxM) */
double *TXYY;      /* output matrix 1x(NxM) */
double *TYXX;      /* output matrix 1x(NxM) */
double *TYXY;      /* output matrix 1x(NxM) */
double *TYYX;      /* output matrix 1x(NxM) */
double *TYYY;      /* output matrix 1x(NxM) */
double *T1;      /* output matrix 1x(NxM) */
double *T2;      /* output matrix 1x(NxM) */
double *D1;      /* output matrix 1x(NxM) */
double *D2;      /* output matrix 1x(NxM) */


/* get the value of the scalar input  */
//multiplier = mxGetScalar(prhs[0]);

/* code here */
/* create a pointer to the real data in the input matrix  */
inMatrixX = mxGetPr(prhs[0]);
inMatrixY = mxGetPr(prhs[1]);
inMatrixX0 = mxGetPr(prhs[2]);
inMatrixY0 = mxGetPr(prhs[3]);
    
/* get dimensions of the input matrix */
nrows = mxGetN(prhs[0]);
ncols = mxGetN(prhs[2]);
    
/* create the output matrix (their dimension depending on the inputs) */
plhs[0] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[4] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[5] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[6] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[7] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[8] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[9] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[10] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[11] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[12] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[13] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[14] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);
plhs[15] = mxCreateDoubleMatrix(1,nrows*ncols-ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
GXX = mxGetPr(plhs[0]);
GXY = mxGetPr(plhs[1]);
GYX = mxGetPr(plhs[2]);
GYY = mxGetPr(plhs[3]);
TXXX = mxGetPr(plhs[4]);
TXXY = mxGetPr(plhs[5]);
TXYX = mxGetPr(plhs[6]);
TXYY = mxGetPr(plhs[7]);
TYXX = mxGetPr(plhs[8]);
TYXY = mxGetPr(plhs[9]);
TYYX = mxGetPr(plhs[10]);
TYYY = mxGetPr(plhs[11]);
T1 = mxGetPr(plhs[12]);
T2 = mxGetPr(plhs[13]);
D1 = mxGetPr(plhs[14]);
D2 = mxGetPr(plhs[15]);

/* call the computational routine */
Stokes2DAxisIntLinearProva(inMatrixX, inMatrixY, inMatrixX0, inMatrixY0, GXX, GXY, GYX, GYY, TXXX, TXXY, TXYX, TXYY, TYXX, TYXY, TYYX, TYYY, T1, T2, D1, D2, ncols, nrows-1);
}
