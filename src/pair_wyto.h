/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Kentaro Hirai, Waseda University
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(wyto, PairWYTO);
// clang-format on
#else

#ifndef LMP_PAIR_WYTO_H
#define LMP_PAIR_WYTO_H

#include "pair.h"
#define PIVAL 3.1415926535898

namespace LAMMPS_NS {

class PairWYTO : public Pair {
 public:
  PairWYTO(class LAMMPS *);
  ~PairWYTO() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template <int SHIFT_FLAG, int EVFLAG, int EFLAG, int VFLAG_ATOM> void eval();

  static constexpr int NPARAMS_PER_LINE = 31;
  struct Param {
    double epsilon, sigma;
    double mu[2], nu[2], xi[2], z0[2];
    double aij[2], aik[2], cij, gammaij[2], gammaik[2], cos0[2], alpha[2];
    double biga, bigb;
    double powerp, powerq;
    double tol;
    double cutpair, cutij[2], cutik[2], cutpairsq, cutijsq[2], cutiksq[2];
    double sigma_gammaij[2], sigma_gammaik[2], lambda_epsilon[2];
    double c1, c2, c3, c4, c5, c6;
    int ielement, jelement, kelement;
  };

  struct Softparam {
    double p1, p2, p3, p4,  p5, p6, p7, p8, p9, bigr, bigd;
    int ielement, jelement;
  };

 protected:
  Param *params;      // parameter set for an I-J-K interaction
  double cutmax;      // max cutoff for all elements
  int maxshort;       // size of short neighbor list array

  Softparam *soft;    // parameter set for the bond softening
  int **softflag;     // flags showing exisistence
  int nsofts;           
  int maxsofts;
  double *numcoord;   // Array holding values ​​corresponding to the coordination number of different atoms
  int **elem2soft;    // mapping from element doubles to bond softening

  int type_index_Si;  // mapping type integers to Si atoms
  int type_index_O;   // mapping type integers to O atoms

  int nmax;

  void settings(int, char **) override;
  virtual void allocate();
  virtual void read_file(char *);
  virtual void setup_params();
  virtual void display_params();

  void twobody_regular(int, int, int, int,
                double, double, double, 
                double, double, double,
                double, double, double,
                double **, int *, int *, int *, int **,   
                Param *, double, double **, int, double &, 
                int, int);
  void twobody_external(int, int, int, int,
                double, double, double, 
                double, double, double,
                double **, int *, int *, int *, int **,   
                Param *, double, double **, int, double &, 
                int, int
  );
  double gSi(double, int, int);
  double gO(double, int, int);
  double d_gSi(double, int, int);
  double d_gO(double, int, int);

  void count_z(int, int, int, int *, double **, int *, int *, int *, int **);
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;  

  double lambda(double, double, double, double, double);
  double d_lambda(double, double, double, double, double);
  double fc(double, int, int);
  double d_fc(double, int, int);
  void threebody(int, int, int, int, 
                  int, int, int,
                  double, double, double,
                  double, double, double,
                  double, double, double,
                  double **, int *, int *, int *, int **,
                  Param *, double, double, double **,
                  int, double &, int, int);
  
};

}    // namespace LAMMPS_NS

#endif
#endif
