#include <stdio.h>
#include <math.h>
#include "mersenne_twister.h"
#include <complex.h>
#include "cmatrix.h"



int main(void){
    int i,j,m;/*必要なものの定義*/
    int seed = 1235;
    int k;
    init_genrand(seed);
    int N = 1000;
    double v = 1.0;
    int l = 5;
    double tau = 0.8*M_PI/(3*v*(double)l);/*0.24*/
    double T = tau*(double)N*(double)l;/*240定数系はあってそう(vを任意にとっていいのか)*/
    
    
    FILE *fp;
    fp = fopen("kadai3.dat","w");
    int d = 300;
    int E = d*d;
    int D = 2*E;
    
    double complex vecC[E],vecC0[E],vecB[E],vecB0[E];
    for(i=0;i<E;++i){
        vecC[i] = cexp(2*M_PI*genrand_real3()*1.0I) / sqrt((double)D);
        vecC0[i] = vecC[i];
        vecB[i] = cexp(2*M_PI*genrand_real3()*1.0I) / sqrt((double)D);
        vecB0[i] = vecB[i];
    }
    


    double DOS,DFT;
    double complex vecD[N+1],naiseki,C;/* vecD[j] = <phi| exp(-ijHt)|phi>*/
    vecD[0] = 1.0;
    for (m = 1;m<N+1;++m){
        for (i = 0;i<l;++i){
            for (j = 0; j < E; ++j){
                
                
                C = vecC[(j+d) % E];/*上*/
                vecC[(j+d) % E] = C * cos(tau*v) - 1.0I * vecB[j] * sin(tau*v);
                vecB[j] = vecB[j] * cos(tau*v) - 1.0I * C * sin(tau*v);
                
                
                if(j % d == d-1){/*右*/
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v) - 1.0I * vecB[j - (d-1)] * sin(tau*v);
                    vecB[j - (d-1)] = vecB[j - (d-1)] * cos(tau*v) - 1.0I * C * sin(tau*v);
                }else{
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v) - 1.0I * vecB[(j+1) % E] * sin(tau*v);
                    vecB[(j+1) % E] = vecB[(j+1) % E] * cos(tau*v) - 1.0I * C * sin(tau*v);
                }
                
                
                C = vecC[j];/*下*/
                vecC[j] = C * cos(tau*v) - 1.0I * vecB[(j+d) % E] * sin(tau*v);
                vecB[(j+d) % E] = vecB[(j+d) % E] * cos(tau*v) - 1.0I * C * sin(tau*v);
            }
            /*for (j = 0; j < E; ++j){
                C = vecC[(j+d) % E];
                vecC[(j+d) % E] = C * cos(tau*v/2) - 1.0I * vecB[j] * sin(tau*v/2);
                vecB[j] = vecB[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                if(j % d == d-1){
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[j - (d-1)] * sin(tau*v/2);
                    vecB[j - (d-1)] = vecB[j - (d-1)] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                }else{
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+1) % E] * sin(tau*v/2);
                    vecB[(j+1) % E] = vecB[(j+1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                }
                C = vecC[j];
                vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+d) % E] * sin(tau*v/2);
                vecB[(j+d) % E] = vecB[(j+d) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
            }*/
            /*for(j = D-1;0<=j;--j){
                C = vecC[j];
                vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+d) % E] * sin(tau*v/2);
                vecB[(j+d) % E] = vecB[(j+d) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                
                if(j % d == d-1){
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[j - (d-1)] * sin(tau*v/2);
                    vecB[j - (d-1)] = vecB[j - (d-1)] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                }else{
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+1) % E] * sin(tau*v/2);
                    vecB[(j+1) % E] = vecB[(j+1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                }
                C = vecC[(j+d) % E];
                vecC[(j+d) % E] = C * cos(tau*v/2) - 1.0I * vecB[j] * sin(tau*v/2);
                vecB[j] = vecB[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
            }*/
        }
        for (i = 0; i < E; ++i){
            naiseki = naiseki + conj(vecC0[i]) * vecC[i] + conj(vecB0[i]) * vecB[i] ;/*c0をcにしてunitaryは確認*/
        }
        vecD[m] = naiseki;
        /*fprintf(fp,"%d %lf %lf \n",m,creal(vecD[m]),cimag(vecD[m]));*/
        naiseki = 0;
    }
    
    

    
    for (k = -1000;k<3.5*T*v/M_PI;++k){/*DFT計算*/
        DFT = 1.0;
        for (j = 1;j<N;++j){
            DFT += creal(2*creal(vecD[j] * cexp(1.0I*M_PI*j*k/(double)N)));
        }
        if (k % 2 == 0){
            DFT += conj(vecD[N]);
        }else{
            DFT -= conj(vecD[N]);
        }
        
        DOS = ( T / (2 * M_PI * (double)N))*DFT;/*DOS計算*/
        fprintf(fp,"%lf  %lf\n",k*M_PI/(v*T),DOS);
    }
    fclose(fp);
    
}


