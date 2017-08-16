//SixTrackLib
//
//Authors: R. De Maria, G. Iadarola, D. Pellegrini, H. Jasim
//
//Copyright 2017 CERN. This software is distributed under the terms of the GNU
//Lesser General Public License version 2.1, copied verbatim in the file
//`COPYING''.
//
//In applying this licence, CERN does not waive the privileges and immunities
//granted to it by virtue of its status as an Intergovernmental Organization or
//submit itself to any jurisdiction.

#ifndef _RFCAV_
#define _RFCAV_

#include "particle.h"
#include "constants.h"

/******************************************/
/*RF Cavity*/
/******************************************/

typedef struct {
  double volt;
  double freq;
  double lag;
} Cavity;

_CUDA_HOST_DEVICE_
int Cavity_track(CLGLOBAL Particle* p, CLGLOBAL Cavity *el){
  double volt = el->volt;
  double freq = el->freq;
  double lag = el->lag;
  double phase, pt, opd;
  phase=lag-2*M_PI/C_LIGHT*freq*p->sigma/p->beta0;
  //printf("ggg00 %e %e\n",p->psigma,p->psigma+p->chi*volt/(p->p0c));
  p->psigma+=p->chi*volt*sin(phase)/(p->p0c*p->beta0);
  pt=p->psigma * p->beta0;
  opd=sqrt( pt*pt+ 2*p->psigma + 1 );
  p->delta=opd - 1;
  p->beta=opd/(1/p->beta0+pt);
  //p->gamma=1/sqrt(1-p->beta*p->beta);
  p->gamma=(pt*p->beta0+1)*p->gamma0;
  p->rpp=1/opd;
  p->rvv=p->beta0/p->beta;
  //printf("ggg2 %e %e %e\n",pt,opd,p->delta);
  return 1;
}

#endif
