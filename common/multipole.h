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

#ifndef _MULTIPOLE_
#define _MULTIPOLE_

#include "particle.h"
#include "constants.h"

/******************************************/
/*Multipole*/
/******************************************/

typedef struct {
  long int order;
  double l ;
  double hxl;
  double hyl;
  double bal[1];
} Multipole;

_CUDA_HOST_DEVICE_
int Multipole_track(CLGLOBAL Particle* p, CLGLOBAL Multipole *el){
  double x,y,chi,dpx,dpy,zre,zim,b1l,a1l,hxx,hyy;
  long int order=el->order;
  double hxl=el->hxl;
  double hyl=el->hyl;
  double l=el->l;
  CLGLOBAL double *bal = el->bal;
  dpx=bal[order*2];
  dpy=bal[order*2+1];
  x=p->x; y=p->y; chi=p->chi;
//  _DP("Multipole_track: dpx,y=%g %G\n",dpx,dpy);
  for (int ii=order-1;ii>=0;ii--){
    zre=(dpx*x-dpy*y);
    zim=(dpx*y+dpy*x);
//    _DP("Multipole_track: y,x=%g %G\n",x,y);
    dpx=bal[ii*2]+zre;
    dpy=bal[ii*2+1]+zim;
//    _DP("Multipole_track: dpx,y=%g %G\n",dpx,dpy);
  }
  dpx=-chi*dpx ;
  dpy=chi*dpy ;
//  _DP("Multipole_track: dpx,y=%g %G\n",dpx,dpy);
  if (l>0){
     b1l=chi*bal[0]; a1l=chi*bal[1];
     hxx=hxl/l*x; hyy=hyl/l*y;
     dpx+=hxl + hxl*p->delta - b1l*hxx;
     dpy-=hyl + hyl*p->delta - a1l*hyy;
     p->sigma-=chi*(hxx-hyy)*l*p->rvv;
  }
  p->px+=dpx ;  p->py+=dpy ;
  return 1 ;
}

#endif
