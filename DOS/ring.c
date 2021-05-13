#include <stdio.h>
#include <math.h>
#include "mersenne_twister.h"
#include <complex.h>
#include "cmatrix.h"



int main(void){
    int i,j,m;/*必要なものの定義*/
    int seed = 12345;
    int k;
    init_genrand(seed);
    int N = 1000;
    double v = 0.5;
    int l = 5;
    double tau = 0.8*M_PI/(2*v*(double)l);
    double T = tau*(double)N*(double)l;
    
    
    FILE *fp;
    fp = fopen("kadai1.dat","w");
    int D = 200000;
    
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
                C = vecC[j];
                vecC[j] = C * cos(tau*v/2) - 1.0I * vecC[(j+1) % D] * sin(tau*v/2);
                vecC[(j+1) % D] = vecC[(j+1) % D] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                }
            for (j = D-1; 0 <= j ; --j){
                C = vecC[j];
                vecC[j] = C * cos(tau*v/2) - 1.0I * vecC[(j+1) % D] * sin(tau*v/2);
                vecC[(j+1) % D] = vecC[(j+1) % D] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                }
        }
        for (i = 0; i < D; ++i){
            naiseki += conj(vecC0[i]) * vecC[i];/*c0をcにしてunitaryは確認*/
        }
        vecD[m] = naiseki;/*reは減衰振動(0.2)、imはランダム(0.02)*/
        naiseki = 0;
    }
    
    

    
    for (k = -1200;k<3*T*v/M_PI;++k){/*DFT計算*/
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


