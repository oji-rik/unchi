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
    double v = 1.0;
    int l = 5;
    double tau = 0.8*M_PI/(4*v*(double)l);
    double T = tau*(double)N*(double)l;
    
    
    FILE *fp;
    fp = fopen("kadai4.dat","w");
    int d = 250;
    int E = d*d;
    int D = 3*E;
    
    double complex vecC[E],vecC0[E],vecB[E],vecB0[E],vecA[E],vecA0[E];
    for(i=0;i<E;++i){
        vecC[i] = cexp(2*M_PI*genrand_real3()*1.0I) / sqrt((double)D);
        vecC0[i] = vecC[i];
        vecB[i] = cexp(2*M_PI*genrand_real3()*1.0I) / sqrt((double)D);
        vecB0[i] = vecB[i];
        vecA[i] =cexp(2*M_PI*genrand_real3()*1.0I) / sqrt((double)D);
        vecA0[i] = vecA[i];
    }
    


    double DOS,DFT;
    double complex vecD[N+1],naiseki,C;/* vecD[j] = <phi| exp(-ijHt)|phi>*/
    vecD[0] = 1.0;
    for (m = 1;m<N+1;++m){
        for (i = 0;i<l;++i){
            for (j = 0; j < E; ++j){
                C = vecB[j];/*同じBC*/
                vecB[j] = C * cos(tau*v/2) - 1.0I * vecC[j] * sin(tau*v/2);
                vecC[j] = vecC[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                C = vecA[j];/*同じAB*/
                vecA[j] = C * cos(tau*v/2) - 1.0I * vecB[j] * sin(tau*v/2);
                vecB[j] = vecB[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                C = vecA[j];/*同じAC*/
                vecA[j] = C * cos(tau*v/2) - 1.0I * vecC[j] * sin(tau*v/2);
                vecC[j] = vecC[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                if((j / d) % 2 == 0){/*奇数段か偶数弾で場合わけ*/
                    C = vecC[j];/*下CA*/
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecA[(j+d) % E] * sin(tau*v/2);
                    vecA[(j+d) % E] = vecA[(j+d) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                    /*printf("%d 下A%d ",j,(j+d) % E);*/
                
                    if(j % d == 0){
                        C = vecC[j];/*下CB*/
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+2*d - 1) % E] * sin(tau*v/2);
                        vecB[(j+2*d - 1) % E] = vecB[(j+2*d - 1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("下B%d ",(j+2*d - 1) % E);*/
                    }else{
                        C = vecC[j];/*下CB*/
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+d-1) % E] * sin(tau*v/2);
                        vecB[(j+d-1) % E] = vecB[(j+d-1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("下B%d ",(j+d-1) % E);*/
                    }
                }else{
                    
                    if(j % d == d-1){/*下CA*/
                        C = vecC[j];
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecA[(j+1) % E] * sin(tau*v/2);
                        vecA[(j+1) % E] = vecA[(j+1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("%d 下A%d ",j,(j+1) % E);*/
                    }else{
                        C = vecC[j];
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecA[(j+d+1) % E] * sin(tau*v/2);
                        vecA[(j+d+1) % E] = vecA[(j+d+1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("%d 下A%d ",j,(j+d+1) % E);*/
                    }
                    C = vecC[j];/*下CB*/
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+d) % E] * sin(tau*v/2);
                    vecB[(j+d) % E] = vecB[(j+d) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                    /*printf("下B%d ",(j+d) % E);*/
                }
                if(j % d == d-1){/*右BA*/
                    C = vecB[j];
                    vecB[j] = C * cos(tau*v) - 1.0I * vecA[j - (d-1)] * sin(tau*v);
                    vecA[j - (d-1)] = vecA[j - (d-1)] * cos(tau*v) - 1.0I * C * sin(tau*v);
                    /*printf("右A%d \n",j - (d-1));*/
                }else{
                    C = vecB[j];
                    vecB[j] = C * cos(tau*v) - 1.0I * vecA[(j+1) % E] * sin(tau*v);
                    vecA[(j+1) % E] = vecA[(j+1) % E] * cos(tau*v) - 1.0I * C * sin(tau*v);
                    /*printf("右A%d \n",(j+1) % E);*/
                }
                if((j / d) % 2 == 0){
                
                    if(j % d == 0){
                        C = vecC[j];/*下CB*/
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+2*d - 1) % E] * sin(tau*v/2);
                        vecB[(j+2*d - 1) % E] = vecB[(j+2*d - 1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("下B%d ",(j+2*d - 1) % E);*/
                    }else{
                        C = vecC[j];/*下CB*/
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+d-1) % E] * sin(tau*v/2);
                        vecB[(j+d-1) % E] = vecB[(j+d-1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("下B%d ",(j+d-1) % E);*/
                    }
                    C = vecC[j];/*下CA*/
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecA[(j+d) % E] * sin(tau*v/2);
                    vecA[(j+d) % E] = vecA[(j+d) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                }else{
                    C = vecC[j];/*下CB*/
                    vecC[j] = C * cos(tau*v/2) - 1.0I * vecB[(j+d) % E] * sin(tau*v/2);
                    vecB[(j+d) % E] = vecB[(j+d) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                    /*printf("下B%d ",(j+d) % E);*/
                    if(j % d == d-1){/*下CA*/
                        C = vecC[j];
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecA[(j+1) % E] * sin(tau*v/2);
                        vecA[(j+1) % E] = vecA[(j+1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("%d 下A%d ",j,(j+1) % E);*/
                    }else{
                        C = vecC[j];
                        vecC[j] = C * cos(tau*v/2) - 1.0I * vecA[(j+d+1) % E] * sin(tau*v/2);
                        vecA[(j+d+1) % E] = vecA[(j+d+1) % E] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                        /*printf("%d 下A%d ",j,(j+d+1) % E);*/
                    }
                }
                C = vecA[j];/*同じAC*/
                vecA[j] = C * cos(tau*v/2) - 1.0I * vecC[j] * sin(tau*v/2);
                vecC[j] = vecC[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                C = vecA[j];/*同じAB*/
                vecA[j] = C * cos(tau*v/2) - 1.0I * vecB[j] * sin(tau*v/2);
                vecB[j] = vecB[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
                C = vecB[j];/*同じBC*/
                vecB[j] = C * cos(tau*v/2) - 1.0I * vecC[j] * sin(tau*v/2);
                vecC[j] = vecC[j] * cos(tau*v/2) - 1.0I * C * sin(tau*v/2);
            }
        }
        for (i = 0; i < E; ++i){
            naiseki = naiseki + conj(vecC0[i]) * vecC[i] + conj(vecB0[i]) * vecB[i] + conj(vecA0[i]) * vecA[i] ;/*c0をcにしてunitaryは確認*/
        }
        vecD[m] = naiseki;
        /*fprintf(fp,"%d %lf %lf \n",m,creal(vecD[m]),cimag(vecD[m]));*/
        naiseki = 0;
    }
    
    

    
    for (k = -1000;k<5*T*v/M_PI;++k){
        DFT = 1.0;
        for (j = 1;j<N;++j){
            DFT += creal(2*creal(vecD[j] * cexp(-1.0I*M_PI*j*k/(double)N)));
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


