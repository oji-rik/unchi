#include <stdio.h>
#include <math.h>
#include "mersenne_twister.h"
#include <complex.h>
#include "cmatrix.h"



int main(void){
    int i,j,m;/*必要なものの定義*/
    int seed = 12345;
    int k;
    double K = 500.;
    init_genrand(seed);
    int N = 10000;
    double v = 1.0;
    int l = 5;
    double tau = 0.8*M_PI/(4*v*(double)l);
    double T = tau*(double)N*(double)l;
    
    
    FILE *fp;
    fp = fopen("kadai31.dat","w");
    int d = 500;
    int D = d*d;
    
    double complex vecC[D],vecC0[D];/*ランダム状態の生成(系数行列)*/
    for(i=0;i<D;++i){
        vecC[i] = cexp(2*M_PI*genrand_real3()*1.0I) / sqrt((double)D);
        vecC0[i] = vecC[i];
    }
    


    double DOS,DFT;
    double complex vecD[N+1],naiseki,C;/* vecD[j] = <phi| exp(-ijHt)|phi>*/
    vecD[0] = 1.0;
    for (m = 1;m<N+1;++m){
        for (i = 0;i<l;++i){
            for (j = 0; j < D; ++j){
                
                
                
                C = vecC[j];/*した*/
                vecC[j] = C * cos(tau*v/2) - 1.0I * vecC[(j+d) % D] * sin(tau*v/2);
                vecC[(j+d) % D] = vecC[(j+d) % D] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                
                
                
                if(j % d == d-1){/*右*/
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * vecC[j - (d-1)] * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                    vecC[j - (d-1)] = vecC[j - (d-1)] * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * C * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                }else{
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * vecC[(j+1) % D] * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                    vecC[(j+1) % D] = vecC[(j+1) % D] * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * C * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                }
            }
            for(j = D-1;0<=j;--j){
                if(j % d == d-1){/*右*/
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * vecC[j - (d-1)] * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                    vecC[j - (d-1)] = vecC[j - (d-1)] * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * C * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                }else{
                    C = vecC[j];
                    vecC[j] = C * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * vecC[(j+1) % D] * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                    vecC[(j+1) % D] = vecC[(j+1) % D] * cos(tau*v*(1-cos(M_PI*K*(j%d)/d))/4) - 1.0I * C * sin(tau*v*(1-cos(M_PI*K*(j%d)/d))/4);
                }
                C = vecC[j];/*した*/
                vecC[j] = C * cos(tau*v/2) - 1.0I * vecC[(j+d) % D] * sin(tau*v/2);
                vecC[(j+d) % D] = vecC[(j+d) % D] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                
            }
        }
        
        for (k = 0; k < D; ++k){
            naiseki += conj(vecC0[k]) * vecC[k];/*c0をcにしてunitaryは確認*/
        }
        vecD[m] = naiseki;
        /*fprintf(fp,"%d %lf %lf \n",m,creal(vecD[m]),cimag(vecD[m]));*/
        naiseki = 0;
    }
    
    

    
    for (k = -10000;k<5*T*v/M_PI;++k){
        DFT = 1.0;
        for (j = 1;j<N;++j){
            DFT += creal(2*creal(vecD[j] * cexp(1.0I*M_PI*j*k/(double)N)));
        }
        if (k % 2 == 0){
            DFT += conj(vecD[N]);
        }else{
            DFT -= conj(vecD[N]);
        }
        DOS = ( T / (2 * M_PI * (double)N))*DFT;
        fprintf(fp,"%lf  %lf\n",k*M_PI/(v*T),DOS);
    }
    fclose(fp);
    
}


