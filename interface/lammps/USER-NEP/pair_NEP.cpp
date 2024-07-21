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

#include "pair_NEP.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "nep.h"
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#define LAMMPS_VERSION_NUMBER 20231121 // use the new neighbor list starting from this version

using namespace LAMMPS_NS;

PairNEP::PairNEP(LAMMPS* lmp) : Pair(lmp)
{
#if LAMMPS_VERSION_NUMBER >= 20201130
  centroidstressflag = CENTROID_AVAIL;
#else
  centroidstressflag = 2;
#endif

  restartinfo = 0;
  manybody_flag = 1;

  single_enable = 0;
  // TODO: TEST
  one_coeff = 1;

  inited = false;
  allocated = 0;
#if defined(FIX_TYPE_LIST) || defined(FIX_TYPE_LIST_BY_LAMMPS)
  type_list_allocated = 0;
#endif
}

PairNEP::~PairNEP()
{
  if (copymode)
    return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
#ifdef FIX_TYPE_LIST
    memory->destroy(type_map);
#endif
#if defined(FIX_TYPE_LIST) || defined(FIX_TYPE_LIST_BY_LAMMPS)
    memory->destroy(type_list);
#endif
  }
}

void PairNEP::allocate()
{
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      setflag[i][j] = 1;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

#ifdef FIX_TYPE_LIST_BY_LAMMPS
  map = new int[n + 1];
#endif
  allocated = 1;
}

void PairNEP::coeff(int narg, char** arg)
{
  if (!allocated)
    allocate();
#ifdef FIX_TYPE_LIST_BY_LAMMPS
  // TODO
  map_element2type(narg-3, arg+3);
  model_filename = arg[2];
#endif
#ifdef FIX_TYPE_LIST
  // create a type map to translate the LAMMPS index to GPUMD index
  int ntype = atom->ntypes;
  int arg_n = narg - 2;
  if (ntype != arg_n)
    error->all(FLERR, "wrong type number: DATA file: {}, IN file: {}", ntype, arg_n);
  
  memory->create(type_map, ntype + 1, "pair_nep:type_map");

  for (int i = 0; i < ntype; ++i) {
    type_map[i + 1] = atoi(arg[i + 2]);
  }
#endif
}

void PairNEP::settings(int narg, char** arg)
{
#ifndef FIX_TYPE_LIST_BY_LAMMPS
  if (narg != 1) {
    error->all(FLERR, "Illegal pair_style command; nep requires 1 parameter");
  }
  model_filename = arg[0];
#endif
}

void PairNEP::init_style()
{
#if LAMMPS_VERSION_NUMBER >= 20220324
  neighbor->add_request(this, NeighConst::REQ_FULL);
#else
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
#endif

  bool is_rank_0 = (comm->me == 0);
  nep_model.init_from_file(model_filename, is_rank_0);
  inited = true;
  cutoff = nep_model.paramb.rc_radial;
  cutoffsq = cutoff * cutoff;
  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      cutsq[i][j] = cutoffsq;
}

double PairNEP::init_one(int i, int j) { return cutoff; }

void PairNEP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag, vflag);
  }
  double total_potential = 0.0;
  double total_virial[6] = {0.0};
  double* per_atom_potential = nullptr;
  double** per_atom_virial = nullptr;
  if (eflag_atom) {
    per_atom_potential = eatom;
  }
  if (cvflag_atom) {
    per_atom_virial = cvatom;
  }

#ifdef FIX_TYPE_LIST
  if (type_list_allocated == 0) {
    // create a type list store GPUMD index
    int atom_number = atom->nmax;
  
    type_list = nullptr;
    memory->create(type_list, atom_number, "pair_nep:type_list");

    int *raw_type = atom->type;
  
    for (int i = 0; i < atom_number; ++i) {
      type_list[i] = type_map[raw_type[i]];
    }

    // change the flag to avoid allocated again
    type_list_allocated = 0;
  }

  nep_model.compute_for_lammps(
    list->inum, list->ilist, list->numneigh, list->firstneigh, type_list,
#ifdef FIX_MOLECULAR
    atom->molecule,
#endif
    atom->x, total_potential,
    total_virial, per_atom_potential, atom->f, per_atom_virial);
#elif defined(FIX_TYPE_LIST_BY_LAMMPS)
  if (type_list_allocated == 0) {
    // create a type list store GPUMD index
    int atom_number = atom->nmax;
  
    type_list = nullptr;
    memory->create(type_list, atom_number, "pair_nep:type_list");

    int *raw_type = atom->type;
  
    for (int i = 0; i < atom_number; ++i) {
      // TODO: should change NEP.cpp to avoid to add 1 here
      type_list[i] = map[raw_type[i]] + 1;
    }

    // change the flag to avoid allocated again
    type_list_allocated = 0;
  }

  nep_model.compute_for_lammps(
    atom->nlocal, list->inum, list->ilist, list->numneigh, list->firstneigh, type_list, 
#ifdef FIX_MOLECULAR
    atom->molecule,
#endif
    atom->x, total_potential,
    total_virial, per_atom_potential, atom->f, per_atom_virial);
#else
  nep_model.compute_for_lammps(
    list->inum, list->ilist, list->numneigh, list->firstneigh, atom->type, atom->x, total_potential,
    total_virial, per_atom_potential, atom->f, per_atom_virial);
#endif

  if (eflag) {
    eng_vdwl += total_potential;
  }
  if (vflag) {
    for (int component = 0; component < 6; ++component) {
      virial[component] += total_virial[component];
    }
  }
#ifdef FIX_TYPE_LIST
  memory->destroy(type_list);
#endif
}
