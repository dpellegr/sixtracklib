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

#ifndef _BB4D_
#define _BB4D_

#include "particle.h"
#include "constants.h"

/******************************************/
/*Beam-beam 4d*/
/******************************************/

#include "transverse_field_gauss_round.h"
#include "transverse_field_gauss_ellip.h"

typedef struct {
    double N_s; // Population strong beam
    double beta_s;
    double q_s;
    long int trasv_field_type; //1: round gaussian
    CLGLOBAL void* field_map_data;
} BB4D_data;



_CUDA_HOST_DEVICE_
int BB4D_track(CLGLOBAL Particle* p, CLGLOBAL BB4D_data *el){
  double Ex, Ey;

  #ifdef DATA_PTR_IS_OFFSET
    CLGLOBAL void * ptr = ((CLGLOBAL uint64_t*) (&(el->field_map_data))) + ((uint64_t) el->field_map_data) + 1;
  #else
    void * ptr = el->field_map_data;
  #endif

  switch(el->trasv_field_type){
    case 1:
      get_transv_field_gauss_round( (CLGLOBAL transv_field_gauss_round_data*) ptr, p->x, p->y, &Ex, &Ey);
      break;
    case 2:
      get_transv_field_gauss_ellip( (CLGLOBAL transv_field_gauss_ellip_data*) ptr, p->x, p->y, &Ex, &Ey);
      break;
    default:
      Ex = 1/0.;
      Ey = 1/0.;
  }

  double fact_kick = p->chi * el->N_s * el->q_s * p->q0 * (1. + p->beta * el->beta_s)/(p->p0c*QELEM*(p->beta + el->beta_s));

  p->px += fact_kick*Ex;
  p->py += fact_kick*Ey;
  return 1;
}


#endif
