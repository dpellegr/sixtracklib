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

#ifndef _LINMAP_
#define _LINMAP_

#include "particle.h"
#include "constants.h"

/******************************************/
/*Linear Map*/
/******************************************/

typedef struct {
      double matrix[8];
      //double disp[4];
} LinMap_data;

_CUDA_HOST_DEVICE_
LinMap_data LinMap_init( double alpha_x_s0, double beta_x_s0, double alpha_x_s1, double beta_x_s1,
                         double alpha_y_s0, double beta_y_s0, double alpha_y_s1, double beta_y_s1,
                         double dQ_x, double dQ_y ) {
  LinMap_data res;
  double s,c;

  //sincos(dQ_x, &s, &c);
  s = sin(2.*M_PI*dQ_x); c = cos(2.*M_PI*dQ_x);
  res.matrix[0] = sqrt(beta_x_s1/beta_x_s0)*(c+alpha_x_s0*s);
  res.matrix[1] = sqrt(beta_x_s1*beta_x_s0)*s;
  res.matrix[2] = ((alpha_x_s0-alpha_x_s1)*c - (1.+alpha_x_s0*alpha_x_s1)*s)/sqrt(beta_x_s1*beta_x_s0);
  res.matrix[3] = sqrt(beta_x_s0/beta_x_s1)*(c-alpha_x_s1*s);

  //sincos(dQ_y, &s, &c);
  s = sin(dQ_y); c = cos(dQ_y);
  res.matrix[4] = sqrt(beta_y_s1/beta_y_s0)*(c+alpha_y_s0*s);
  res.matrix[5] = sqrt(beta_y_s1*beta_y_s0)*s;
  res.matrix[6] = ((alpha_y_s0-alpha_y_s1)*c - (1.+alpha_y_s0*alpha_y_s1)*s)/sqrt(beta_y_s1*beta_y_s0);
  res.matrix[7] = sqrt(beta_y_s0/beta_y_s1)*(c-alpha_y_s1*s);
  return res;
}

_CUDA_HOST_DEVICE_
int LinMap_track(CLGLOBAL Particle* p, CLGLOBAL LinMap_data *el){
  double M00 = el->matrix[0];
  double M01 = el->matrix[1];
  double M10 = el->matrix[2];
  double M11 = el->matrix[3];
  double M22 = el->matrix[4];
  double M23 = el->matrix[5];
  double M32 = el->matrix[6];
  double M33 = el->matrix[7];
  double x0  = p->x;
  double px0 = p->px;
  double y0  = p->y;
  double py0 = p->py;

  p->x  = M00*x0 + M01*px0;
  p->px = M10*x0 + M11*px0;
  p->y  = M22*y0 + M23*py0;
  p->py = M32*y0 + M33*py0;
  return 1;
}

#endif
