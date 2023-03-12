/************************************************************************
 ** file    : azimuth-tide.c
 ** summery : estimate 2D azimuth of tidal strain
 ** author  : T.Takano
 ** date    : 2017/04/25
 ** memo    : --
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "nrutil.h"


#define FREE_ARG char*
#define DEBUG

double pi;
int main(int argc, char *argv[])
{
  int num;
  int jj;
  int i,j; //loop counter
  double rad;
  double *exx_ptr=NULL;
  double *eyy_ptr=NULL;
  double *exy_ptr=NULL;
  double *ezz_ptr=NULL;
  double *eaz_ptr=NULL;
  double *eareal_ptr=NULL;
  double *ecubic_ptr=NULL;
  double dila;
  double shear;
  double *en_ptr=NULL;
  double *ep_ptr=NULL;
  double *emax_ptr=NULL;
  double *emin_ptr=NULL;
  double *thmax_ptr=NULL;
  double *thmin_ptr=NULL;
  double sumxx,sumyy,sumxy,avexx,aveyy,avexy;
  double dsumxx,dsumyy,dsumxy,davexx,daveyy,davexy;
  char yymmdd[100];
  char hhmm[100];
  FILE *fp;
  int nn=0;

  pi = M_PI;
  // file open
  if(argc != 3){
    printf("set the open file\n");
    printf("usage: ./a.out filelist total-num.  \n");
    exit(EXIT_FAILURE);
  }

  //input file
  if((fp = fopen(argv[1], "r")) == NULL){
    printf("ReadFile open error \n");
    exit(EXIT_FAILURE);
  }
  
  //input number
  sscanf(argv[2], "%d", &num);

  //memory allocation
  exx_ptr=(double *)dvector(0,num);
  eyy_ptr=(double *)dvector(0,num);
  exy_ptr=(double *)dvector(0,num);
  ezz_ptr=(double *)dvector(0,num);
  eaz_ptr=(double *)dvector(0,num);
  eareal_ptr=(double *)dvector(0,num);
  ecubic_ptr=(double *)dvector(0,num);
  emax_ptr=(double *)dvector(0,num);
  emin_ptr=(double *)dvector(0,num);
  en_ptr=(double *)dvector(0,360);
  ep_ptr=(double *)dvector(0,360);
  thmax_ptr=(double *)dvector(0,num);
  thmin_ptr=(double *)dvector(0,num);

  //read file
  for(i=0;i<num;++i){
    fscanf(fp, "%s %s %lf %lf %lf %lf %lf %lf", yymmdd, hhmm, eyy_ptr+i, exx_ptr+i, exy_ptr+i, eaz_ptr+i, eareal_ptr+i, ecubic_ptr+i);
    *(ezz_ptr+i) = *(ecubic_ptr+i) - *(eareal_ptr+i);
    dila = *(exx_ptr+i) + *(eyy_ptr+i);
    shear= sqrt(exy_ptr[i]*exy_ptr[i]+pow(exx_ptr[i]-eyy_ptr[i],2.0)*0.25);
    *(emax_ptr+i) = dila*0.5+shear;
    *(emin_ptr+i) = dila*0.5-shear;
    *(thmax_ptr+i) = atan((emax_ptr[i]-exx_ptr[i])/exy_ptr[i])*180/pi;
    *(thmin_ptr+i) = atan((emin_ptr[i]-exx_ptr[i])/exy_ptr[i])*180/pi;
  }  

  sumxx=0.0;
  sumyy=0.0;
  sumxy=0.0;
  nn=0;
  for(i=0;i<num;++i){
    if(emin_ptr[i]<-2E-8){
      sumxx+= exx_ptr[i];
      sumyy+= eyy_ptr[i];
      sumxy+= exy_ptr[i];
      nn++;
      //printf("%le\n",emin_ptr[i]);
    }
  }
  avexx = sumxx/nn;  
  aveyy = sumyy/nn;  
  avexy = sumxy/nn;  


  dsumxx=0.0;
  dsumyy=0.0;
  dsumxy=0.0;
  nn=0;
  for(i=0;i<num;++i){
    if(emax_ptr[i]>2E-8){
      dsumxx+= exx_ptr[i];
      dsumyy+= eyy_ptr[i];
      dsumxy+= exy_ptr[i];
      nn++;
      //printf("%le\n",emin_ptr[i]);
    }
  }
  davexx = dsumxx/nn;
  daveyy = dsumyy/nn;
  davexy = dsumxy/nn;


  for(j=0;j<360;++j){
    jj = 90-j;
    rad = (double)j*pi/180.0;
    en_ptr[j] = cos(rad)*cos(rad)*avexx + sin(rad)*sin(rad)*aveyy + 2.0*cos(rad)*sin(rad)*avexy;
  }

  for(j=0;j<360;++j){
    jj = 90-j;
    rad = (double)j*pi/180.0;
    ep_ptr[j] = cos(rad)*cos(rad)*davexx + sin(rad)*sin(rad)*daveyy + 2.0*cos(rad)*sin(rad)*davexy;
  }

  /*
  for(i=0;i<num;++i){
  for(j=0;j<360;++j){
  rad = (double)j*pi/180.0;
  en_ptr[i][j] = cos(rad)*cos(rad)*exx_ptr[i] + sin(rad)*sin(rad)*eyy_ptr[i] + 2.0*cos(rad)*sin(rad)*exy_ptr[i];
  }
  }*/
  
  for(j=0;j<360;++j){
    printf("%lf %le %le\n",(double)j, ep_ptr[j], en_ptr[j]);
  }

  free_dvector(en_ptr,0,360);
  free_dvector(ep_ptr,0,360);
  free_dvector(exx_ptr,0,num);
  free_dvector(eyy_ptr,0,num);
  free_dvector(exy_ptr,0,num);
  free_dvector(ezz_ptr,0,num);
  free_dvector(eaz_ptr,0,num);
  free_dvector(eareal_ptr,0,num);
  free_dvector(ecubic_ptr,0,num);
  free_dvector(emax_ptr,0,num);
  free_dvector(emin_ptr,0,num);
  free_dvector(thmax_ptr,0,num);
  free_dvector(thmin_ptr,0,num);

  return 0;
}
