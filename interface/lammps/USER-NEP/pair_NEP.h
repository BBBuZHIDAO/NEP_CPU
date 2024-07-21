/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Junjie Wang (Nanjing University)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(nep, PairNEP)

#else

#ifndef LMP_PAIR_NEP_H
#define LMP_PAIR_NEP_H

#include "nep.h"
#include "pair.h"
#include <string>

// #define FIX_TYPE_LIST 1
#define FIX_TYPE_LIST_BY_LAMMPS 1

namespace LAMMPS_NS
{
class PairNEP : public Pair
{
public:
  double cutoff;
  NEP3 nep_model;
#ifdef FIX_TYPE_LIST
  int *type_map;
#endif
#if defined(FIX_TYPE_LIST) || defined(FIX_TYPE_LIST_BY_LAMMPS)
  int *type_list;
  int type_list_allocated;
#endif
  PairNEP(class LAMMPS*);
  virtual ~PairNEP();
  virtual void coeff(int, char**);
  virtual void settings(int, char**);
  virtual double init_one(int, int);
  virtual void init_style();
  virtual void compute(int, int);

protected:
  bool inited;
  std::string model_filename;
  double cutoffsq;
  void allocate();
};
} // namespace LAMMPS_NS

#endif
#endif
