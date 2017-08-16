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

#ifndef _ALIGN_
#define _ALIGN_

#include "particle.h"
#include "constants.h"

/******************************************/
/*Align*/
/******************************************/

typedef struct {
  double cz;
  double sz;
  double dx;
  double dy;
} Align;

_CUDA_HOST_DEVICE_
int Align_track(CLGLOBAL Particle* p, CLGLOBAL Align *el){
  double xn,yn;
  double cz = el->cz;
  double sz = el->sz;
  double dx = el->dx;
  double dy = el->dy;
  xn= cz*p->x-sz*p->y - dx;
  yn= sz*p->x+cz*p->y - dy;
  p->x=xn;
  p->y=yn;
  xn= cz*p->px+sz*p->py;
  yn=-sz*p->px+cz*p->py;
  p->px=xn;
  p->py=yn;
  return 1;
}

#endif
