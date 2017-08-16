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

#ifndef _DRIFT_
#define _DRIFT_

#include "particle.h"
#include "constants.h"

/******************************************/
/*Drift*/
/******************************************/

typedef struct {
  double length;
} Drift ;

_CUDA_HOST_DEVICE_
int Drift_track(CLGLOBAL Particle* p, CLGLOBAL Drift *el){
  double xp, yp;
  double length=el->length;
  xp = p->px * p->rpp;
  yp = p->py * p->rpp;
  p->x += xp * length;
  p->y += yp * length;
  p->sigma += length * (1 - p->rvv*( 1 + (xp*xp+yp*yp)/2 ) );
  p->s+=length;
//  _DP("Drift_track: length=%g\n",length);
  return 1;
}

/******************************************/
/*Exact Drift*/
/******************************************/

typedef struct {
  double length;
} DriftExact ;

_CUDA_HOST_DEVICE_
int DriftExact_track(CLGLOBAL Particle* p, CLGLOBAL DriftExact *el){
  double lpzi, lbzi, px, py, opd;
  double length = el->length;
  opd=1+p->delta;
  px=p->px; py=p->py;
  lpzi= length/sqrt(opd*opd-px*px-py*py);
  lbzi=(p->beta0*p->beta0*p->psigma+1)*lpzi;
  p->x += px*lpzi ;
  p->y += py*lpzi ;
  p->sigma += length - lbzi;
  p->s += length ;
  return 1;
}

#endif
