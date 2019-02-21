//
//  Stokes2DAxisSPlinesLinear.c
//  
//
//  Created by Giacomo on 7/2/15.
//
//

#include "mex.h"
#include "math.h"
#include "stdio.h"

void Stokes2DAxisSPlinesLinear(double *aax, double *bbx, double *ccx, double *ddx, double *aay, double *bby, double *ccy, double *ddy, double *x0, double *y0, double *GXX, double *GXY, double *GYX, double *GYY, double *A11, double *A12, double *A21, double *A22, int field, int elem)
{
  int i, j, dolog1, dolog2, axis;
  double SXX, SXY, SYX, SYY;    //single layer
  double QXXX, QXXY, QXYX, QXYY, QYXX, QYXY, QYYX, QYYY;    //double layer
  double PXX, PXY, PYX, PYY;    //desingularization double layer
  double F, E, G, B, D, C, k, k2;   //elliptic integrals
  double I10, I11, I30, I31, I32, I50, I51, I52, I53;
  double ax, bx, cx, dx, ay, by, cy, dy, h, h0, h1;
  double r, r0, d, d2, d3, d5, GPpar, phiA, phiB, eta, beta, deta, dbeta;
  double k2p, k4, k5, k6, yy5;
  double FCTR, RL30, RL50, RL52, RL54, RL56;
  double DX, DY, DX2, DX3, r2, r3, r02, r03, nx, ny;
  double distXA, distYA, distA, distXB, distYB, distB, Xleft, Yleft, Xright, Yright;
  double *GXXa, *GXYa, *GYXa, *GYYa, *A11a, *A12a, *A21a, *A22a;
  double *GXXb, *GXYb, *GYXb, *GYYb, *A11b, *A12b, *A21b, *A22b;
  
  //memory allocation
  GXXa = (double*)calloc(field*elem,sizeof(double));
  GXYa = (double*)calloc(field*elem,sizeof(double));
  GYXa = (double*)calloc(field*elem,sizeof(double));
  GYYa = (double*)calloc(field*elem,sizeof(double));
  GXXb = (double*)calloc(field*elem,sizeof(double));
  GXYb = (double*)calloc(field*elem,sizeof(double));
  GYXb = (double*)calloc(field*elem,sizeof(double));
  GYYb = (double*)calloc(field*elem,sizeof(double));
  A11a = (double*)calloc(field*elem,sizeof(double));
  A12a = (double*)calloc(field*elem,sizeof(double));
  A21a = (double*)calloc(field*elem,sizeof(double));
  A22a = (double*)calloc(field*elem,sizeof(double));
  A11b = (double*)calloc(field*elem,sizeof(double));
  A12b = (double*)calloc(field*elem,sizeof(double));
  A21b = (double*)calloc(field*elem,sizeof(double));
  A22b = (double*)calloc(field*elem,sizeof(double));
  
  //FILE *pFile1;
  //FILE *pFile2;
  //pFile1 = fopen ("myfile1.txt","w");
  //pFile2 = fopen ("myfile2.txt","w");
  
  //define pi and tol
  double pi = M_PI;
  double tol = 0.0000000000001;
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
  
  for (j=0; j<elem; j++) {  
      //elemnts variables      
      ax = aax[j];  bx = bbx[j];    cx = ccx[j];    dx = ddx[j];
      ay = aay[j];  by = bby[j];    cy = ccy[j];    dy = ddy[j];
      
      //metrics for first and last term
      h0 = sqrt(bx*bx+by*by);
      h1 =  sqrt((bx+2*cx+3*dx)*(bx+2*cx+3*dx)+(by+2*cy+3*dy)*(by+2*cy+3*dy));
      
      //first and last point of the element
      Xleft = ax;   Yleft = ay;     Xright = ax+bx+cx+dx;   Yright = ay+by+cy+dy;
      
      for (i=0; i<field; i++) {
          
          //check if I'm on the axis
          axis = (y0[i] < 10e-8);
          
          //decide if doing singular treatment
          distXA = Xleft-x0[i];
          distYA = Yleft-y0[i];
          distA = sqrt(distXA*distXA+distYA*distYA);
          dolog1 = (distA>10e-8);
          //dolog1 = 1;
          
          //decide if doing singular treatment
          distXB = Xright-x0[i];
          distYB = Yright-y0[i];
          distB = sqrt(distXB*distXB+distYB*distYB);
          dolog2 = (distB>10e-8);
          //dolog2 = 1;
          
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
                r0 = y0[i];
                r = beta;
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
                SXX = r*(I10+DX*DX*I30);
                SXX += (2*log(d) - 2*(log(d/GPpar) + log(GPpar)*(1-h0/h)))*(1-dolog1)*(1-axis);     //remove singularity on the left
                SXX += (2*log(d) - 2*(log(d/(1-GPpar)) + log(1-GPpar)*(1-h1/h)))*(1-dolog2)*(1-axis);     //remove singularity on the right
                
                SXY = r*DX*(r*I30-r0*I31);
                SYX = r*DX*(r*I31-r0*I30);
                
                SYY = r*(I11+(r*r+r0*r0)*I31-r*r0*(I30+I32));
                SYY += (2*log(d) - 2*(log(d/GPpar) + log(GPpar)*(1-h0/h)))*(1-dolog1)*(1-axis);     //remove singularity on the left
                SYY += (2*log(d) - 2*(log(d/(1-GPpar)) + log(1-GPpar)*(1-h1/h)))*(1-dolog2)*(1-axis);     //remove singularity on the right
                                
                //double layer
                QXXX = -6*r*DX3*I50;
                QXXY = -6*r*DX2*(r*I50-r0*I51);
                QXYX = QXXY;
                QXYY = -6*r*DX*(r0*r0*I52+r*r*I50-2*r*r0*I51);
                QYXX = -6*r*DX2*(r*I51-r0*I50);
                QYXY = -6*r*DX*((r*r+r0*r0)*I51-r*r0*(I50+I52));
                QYYX = QYXY;
                QYYY = -6*r*(r3*I51-r2*r0*(I50+2*I52)+r*r0*r0*(I53+2*I51)-r0*r0*r0*I52);
              
                //integration phiA SINGLE LAYER
                GXXa[j*field+i] += phiA*GW[l]*SXX/2 * h;
                GXYa[j*field+i] += phiA*GW[l]*SXY/2 * h;
                GYXa[j*field+i] += phiA*GW[l]*SYX/2 * h;
                GYYa[j*field+i] += phiA*GW[l]*SYY/2 * h;
                
                //integration phiA DOUBLE LAYER
                A11a[j*field+i] += phiA*GW[l]*(QXXX*nx+QXXY*ny)*h/2;
                A12a[j*field+i] += phiA*GW[l]*(QXYX*nx+QXYY*ny)*h/2;
                A21a[j*field+i] += phiA*GW[l]*(QYXX*nx+QYXY*ny)*h/2;
                A22a[j*field+i] += phiA*GW[l]*(QYYX*nx+QYYY*ny)*h/2;
                
                //integration phiB SINGLE LAYER
                GXXb[j*field+i] += phiB*GW[l]*SXX/2 * h;
                GXYb[j*field+i] += phiB*GW[l]*SXY/2 * h;
                GYXb[j*field+i] += phiB*GW[l]*SYX/2 * h;
                GYYb[j*field+i] += phiB*GW[l]*SYY/2 * h;
                
                //integration phiB DOUBLE LAYER
                A11b[j*field+i] += phiB*GW[l]*(QXXX*nx+QXXY*ny)*h/2;
                A12b[j*field+i] += phiB*GW[l]*(QXYX*nx+QXYY*ny)*h/2;
                A21b[j*field+i] += phiB*GW[l]*(QYXX*nx+QYXY*ny)*h/2;
                A22b[j*field+i] += phiB*GW[l]*(QYYX*nx+QYYY*ny)*h/2;
    
          }
                    
          //singular treatment: add part from analytical integration (as in research notes) IF I'M ON THE AXIS I DON'T NEED SINGULAR TREATMENT
          GXXa[j*field+i] += 1.5*h0*(1-dolog1)*(1-axis);
          GYYa[j*field+i] += 1.5*h0*(1-dolog1)*(1-axis);
          
          GXXb[j*field+i] += 0.5*h0*(1-dolog1)*(1-axis);
          GYYb[j*field+i] += 0.5*h0*(1-dolog1)*(1-axis);
          
          GXXa[j*field+i] += 0.5*h1*(1-dolog2)*(1-axis);
          GYYa[j*field+i] += 0.5*h1*(1-dolog2)*(1-axis);
          
          GXXb[j*field+i] += 1.5*h1*(1-dolog2)*(1-axis);
          GYYb[j*field+i] += 1.5*h1*(1-dolog2)*(1-axis);
          
          //build output matrix
          GXX[j*field+i] += GXXa[j*field+i];
          GXX[(j+1)*field+i] += GXXb[j*field+i];
          
          GXY[j*field+i] += GXYa[j*field+i];
          GXY[(j+1)*field+i] += GXYb[j*field+i];
          
          GYX[j*field+i] += GYXa[j*field+i];
          GYX[(j+1)*field+i] += GYXb[j*field+i];
          
          GYY[j*field+i] += GYYa[j*field+i];
          GYY[(j+1)*field+i] += GYYb[j*field+i];
          
          A11[j*field+i] += A11a[j*field+i];
          A11[(j+1)*field+i] += A11b[j*field+i];
          
          A12[j*field+i] += A12a[j*field+i];
          A12[(j+1)*field+i] += A12b[j*field+i];
          
          A21[j*field+i] += A21a[j*field+i];
          A21[(j+1)*field+i] += A21b[j*field+i];
          
          A22[j*field+i] += A22a[j*field+i];
          A22[(j+1)*field+i] += A22b[j*field+i];
          
      }
  }
  
  //free memory
  free(GXXa);
  free(GXYa);
  free(GYXa);
  free(GYYa);
  free(A11a);
  free(A12a);
  free(A21a);
  free(A22a);
  free(GXXb);
  free(GXYb);
  free(GYXb);
  free(GYYb);
  free(A11b);
  free(A12b);
  free(A21b);
  free(A22b);
  
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
double *GYX;      /* output matrix 1x(NxM) */
double *GYY;      /* output matrix 1x(NxM) */
double *A11;      /* output matrix 1x(NxM) */
double *A12;      /* output matrix 1x(NxM) */
double *A21;      /* output matrix 1x(NxM) */
double *A22;      /* output matrix 1x(NxM) */

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
nrows = mxGetN(prhs[0])+1;
ncols = mxGetN(prhs[8]);
    
/* create the output matrix (their dimension depending on the inputs) */
plhs[0] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[4] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[5] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[6] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);
plhs[7] = mxCreateDoubleMatrix(1,nrows*ncols,mxREAL);

/* get a pointer to the real data in the output matrix */
GXX = mxGetPr(plhs[0]);
GXY = mxGetPr(plhs[1]);
GYX = mxGetPr(plhs[2]);
GYY = mxGetPr(plhs[3]);
A11 = mxGetPr(plhs[4]);
A12 = mxGetPr(plhs[5]);
A21 = mxGetPr(plhs[6]);
A22 = mxGetPr(plhs[7]);

/* call the computational routine */
Stokes2DAxisSPlinesLinear(ax, bx, cx, dx, ay, by, cy, dy, inMatrixX0, inMatrixY0, GXX, GXY, GYX, GYY, A11, A12, A21, A22, ncols, nrows-1);
}
