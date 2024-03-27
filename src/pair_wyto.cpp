// clang-format off
/* ----------------------------------------------------------------------
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

#include "pair_wyto.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;

#define DELTA 4
#define THREEBODY_NMAX 2
#define TWOBODY_REG 1
#define TWOBODY_EXT 1
#define THREEBODY 1

//#define DEBUG
/*
#ifdef DEBUG
#define LOG(message) std::cout << message << std::endl;
#define ROOT_LOG(proc_num, message) if (proc_num == 0) std::cout << message << std::endl;
#else
#define LOG(message)
#define ROOT_LOG(proc_num, messege)
#endif
#define PRINT(message) std::cout << message << std::endl;
#define ROOT_PRINT(proc_num, message) if (proc_num == 0) std::cout << message << std::endl;
*/
//using namespace std;


/* ---------------------------------------------------------------------- */

PairWYTO::PairWYTO(LAMMPS *lmp) : Pair(lmp)
{
  /* [Required] Constructor */

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  nmax = 0;
  comm_forward = 1; // flag for communicating coordination number

  params = nullptr;
  soft = nullptr;
  elem2soft = nullptr;

}

PairWYTO::~PairWYTO()
{
  /* [Required] Destructor*/

  //  check if allocated, since class can be destructed when incomplete

  if (copymode) return;

  memory->destroy(params);
  memory->destroy(soft);
  memory->destroy(elem3param);
  memory->destroy(elem2soft);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(softflag);
    memory->destroy(numcoord);
  }

  //LOG("Exit destructor");
}

/* ---------------------------------------------------------------------- */

void PairWYTO::compute(int eflag, int vflag)
{
    //LOG("Exit init_one");
  //ROOT_LOG(comm->me, "Start compute");

  int i, j, k, n, ii, jj, kk, inum, jnum, jnumm1, itag, jtag;
  int itype, jtype, ktype, ijparam, ikparam, ijkparam;
  double xtmpi, ytmpi, ztmpi, xtmpj, ytmpj, ztmpj, delx, dely, delz, evdwl, fpair;
  double rsq, rsq1, rsq2, bigr, bigd, gij, r;
  double delr1[3], delr2[3], delr3[3];
  int *ilist, *jlist, *numneigh, **firstneigh;
  
  evdwl = 0.0;
  ev_init(eflag,vflag);
  
  if (atom->nmax > nmax) {
    // allocate memory for numcoord array
    // size:The maximum number of atoms handled by each process

    nmax = atom->nmax;
    memory->create(numcoord, nmax, "pair:numcoord");
  }

  double **x = atom->x;                   // x[i][j]:coorddinate of local atom i (j component)
  double **f = atom->f;                   // f[i][j]:force of local atom i (j component)
  int *tag = atom->tag;                   // tag[i]:global index of atom i
  int *type = atom->type;                 // map[type[i]]:atom type integer of atom i
  int nlocal = atom->nlocal;              // number of atoms in responsible subdomain
  int nall = nlocal + atom->nghost;       // number of atoms in responsible subdomain and ghost atom
  
  // if newton_pair == 1: 
  //  Forces must be stored even with ghost atoms and after all forces are computed
  //  a “reverse communication” is performed to add those ghost atom forces to their 
  //  corresponding local atoms. (Default:newton_pair = 1)
  int newton_pair = force->newton_pair; 

  inum = list->inum;                // number of neiborlists in responsible subdomain
  ilist = list->ilist;              // array holding the local index of the atoms for which neighbor list was generated
  numneigh = list->numneigh;        // length of neighbor list of atom i
  firstneigh = list->firstneigh;    // firstneigh[i]:the first pointer of the neighbor list of atom i (local)

  // calculate coordination number of **local** atoms（update numcoord）
  count_z(nall, nmax, inum, ilist, x, tag, type, numneigh, firstneigh);

  for(ii = 0; ii < inum; ii++){
    i = ilist[ii];           // get index of a local atom i
    itag = tag[i];           // get global index of the atom i
    itype = map[type[i]];    // get atom type of the atom i
    
    // store coordinate of the atom i in buffer
    xtmpi = x[i][0];
    ytmpi = x[i][1];
    ztmpi = x[i][2];

    jlist = firstneigh[i];    // get neighbor list for the atom i
    jnum = numneigh[i];       // get length of the neighbor list

    // loop for twobody energy and regular force
    for(jj = 0; jj < jnum; jj++){
      j = jlist[jj];  // get atom index j from neighbot list of atom i
      j &= NEIGHMASK; // remove the flag indicating the state of the bond contained in the upper 2 bits of j

      if (j == i) continue;       // condtion for Σ

      // store coordinate of the atom j in buffer
      xtmpj = x[j][0];
      ytmpj = x[j][1];
      ztmpj = x[j][2];

      // calculate rij^2[Å^2]
      delx = xtmpi - xtmpj;
      dely = ytmpi - ytmpj;
      delz = ztmpi - ztmpj;
      rsq = delx*delx + dely*dely + delz*delz;

      jtype = map[type[j]]; // get atom type of the atom j

      // get the index of parameters for two-body interaction between ij
      ijparam = elem3param[itype][jtype][jtype];

      // Skip half of full neighbor list
      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (ztmpj < ztmpi) continue;
        if (ztmpj == ztmpi && ytmpj < ytmpi) continue;
        if (ztmpj == ztmpi && ytmpj == ytmpi && xtmpj < xtmpi) continue;
      }
      
      // squared cutoff 
      if (rsq >= params[ijparam].cutpairsq) continue;

      // for debug
      if (!TWOBODY_REG) continue;

      // calculate twobody energy and regular force
      twobody_regular(i, j, itype, jtype,
              delx, dely, delz,
              xtmpi, ytmpi, ztmpi,
              xtmpj, ytmpj, ztmpj,
              x, tag, type, numneigh, firstneigh,
              &params[ijparam], rsq, f, eflag, evdwl,
              nlocal, newton_pair
              );
      
    }

    // // loop for twobody external force
    for (jj = 0; jj < jnum; jj++){
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      if (softflag[itype][jtype] && TWOBODY_EXT) {

        // calculate rij^2[Å^2]
        delx = xtmpi - x[j][0];
        dely = ytmpi - x[j][1];
        delz = ztmpi - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        ijparam = elem3param[itype][jtype][jtype];

        if (rsq >= params[ijparam].cutpairsq) {
          continue;
        }

        // calculate twoboy external force
        twobody_external(i, j, itype, jtype,
              delx, dely, delz,
              xtmpi, ytmpi, ztmpi,
              x, tag, type, numneigh, firstneigh,
              &params[ijparam], rsq, f, eflag, evdwl,
              nlocal, newton_pair);
      }
    }

    // loop for threebody energy and force
    jnumm1 = jnum - 1;
    for (n = 0; n < THREEBODY_NMAX; n++){
      for (jj = 0; jj < jnumm1; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        if (j == i) continue; // condition for Σ

        jtype = map[type[j]];
        xtmpj = x[j][0];
        ytmpj = x[j][1];
        ztmpj = x[j][2];

        for (kk = jj + 1; kk < jnum; kk++) {
          k = jlist[kk];
          k &= NEIGHMASK;

          if (k == i) continue; // condition for Σ

          ktype = map[type[k]];
          ijkparam = elem3param[itype][jtype][ktype];

          // calculate rij^2[Å^2]
          delr1[0] = xtmpi - xtmpj;
          delr1[1] = ytmpi - ytmpj;
          delr1[2] = ztmpi - ztmpj;
          rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

          if (rsq1 >= params[ijkparam].cutijsq[n]) continue;

          // calculate rik^2[Å^2]
          delr2[0] = xtmpi - x[k][0];
          delr2[1] = ytmpi - x[k][1];
          delr2[2] = ztmpi - x[k][2];
          rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

          if (rsq2 >= params[ijkparam].cutiksq[n]) continue;
          
          // calculate each components of rjk
          delr3[0] = xtmpj - x[k][0];
          delr3[1] = ytmpj - x[k][1];
          delr3[2] = ztmpj - x[k][2];

          // for debug
          if (!THREEBODY) continue;

          // calculate threebody energy and force
          threebody(i, j, k, n, 
                  itype, jtype, ktype, 
                  delr1[0], delr1[1], delr1[2], 
                  delr2[0], delr2[1], delr2[2],
                  delr3[0], delr3[1], delr3[2],
                  x, tag, type, numneigh, firstneigh, 
                  &params[ijkparam], rsq1, rsq2, f, eflag, evdwl,
                  nlocal, newton_pair);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairWYTO::allocate()
{
  /* allocate memory (for "coeff" function) */

  allocated = 1;
  int n = atom->ntypes;           // # of atom type
  int num = ceil(atom->natoms);   // # of atom

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(softflag, n+1, n+1, "pair:softflag");
  map = new int[n+1];
}

void PairWYTO::settings(int narg, char **/*arg*/)
{
  /* [Required] processes the arguments to the pair_style command */
  
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

void PairWYTO::coeff(int narg, char **arg)
{
  /* [Required] set coeffs for one or more type pairs */

  if (!allocated) allocate();
  map_element2type(narg-3,arg+3);

  // read potential file and initialize potential parameters

  ROOT_LOG(comm->me, "#### DEBUG MODE ####");
  ROOT_LOG(comm->me, "  TWOBODY_REG:" << TWOBODY_REG);
  ROOT_LOG(comm->me, "  TWOBODY_EXT:" << TWOBODY_EXT);
  ROOT_LOG(comm->me, "  THREEBODY:" << THREEBODY);
  ROOT_LOG(comm->me, "####################");

  ROOT_PRINT(comm->me, "Reading paramter file ...");
  read_file(arg[2]);
  setup_params();
  display_params();
}

void PairWYTO::init_style()
{
  /* init specific to this pair style */

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style WYTO requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style WYTO requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this,NeighConst::REQ_FULL);

  //LOG("Exit init_style");
}

double PairWYTO::init_one(int i, int j)
{
  //LOG("Enter init_one");

  /* init for one type pair i,j and corresponding j,i */
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

void PairWYTO::read_file(char *file)
{
  /* read parameter file */

  memory->sfree(params);
  memory->sfree(soft);
  params = nullptr;
  nparams = maxparam = 0;
  soft = nullptr;
  nsofts = maxsofts = 0;

  // initialize softflag
  for (int y = 0; y < atom->ntypes; y++){
    for (int z = 0; z < atom->ntypes; z++){
      softflag[y][z] = 0;
    }
  }

  // open file on proc 0

  if (comm->me == 0) {
    // load parameter file with specified extension (wyto)

    // Parameters are read only by the root process and are later 
    // broadcast to other processes.

    PotentialFileReader reader(lmp, file, "wyto", unit_convert_flag);
    char *line;
    int param_type = 0;
    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,unit_convert);

    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        // determine whether the line contains the Soft parameter
        if(values.contains("SOFT")){
          param_type = 1;
          values.skip(1);  // Skip token "SOFT".
        }

        std::string iname = values.next_string();
        std::string jname = values.next_string();

        int ielement, jelement;

        if(param_type == 0){

          // load the parameters for twobody and threebody

          std::string kname = values.next_string();
          int kelement;

          // ielement,jelement,kelement = 1st args
          // if all 3 args are in element list, then parse this line
          // else skip to next entry in file
          for (ielement = 0; ielement < nelements; ielement++)
            if (iname == elements[ielement]) break;
          if (ielement == nelements) continue;
          for (jelement = 0; jelement < nelements; jelement++)
            if (jname == elements[jelement]) break;
          if (jelement == nelements) continue;
          for (kelement = 0; kelement < nelements; kelement++)
            if (kname == elements[kelement]) break;
          if (kelement == nelements) continue;

          // load up parameter settings and error check their values

          // Dynamically expand the memory area of ​​the soft array 
          //    according to the number of loaded parameter sets
          if (nparams == maxparam) {
            maxparam += DELTA;
            params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");

            // make certain all addional allocated storage is initialized
            // to avoid false positives when checking with valgrind

            memset(params + nparams, 0, DELTA*sizeof(Param));
          }

          // Load 31 tokens
          params[nparams].ielement    = ielement;
          params[nparams].jelement    = jelement;
          params[nparams].kelement    = kelement;
          params[nparams].epsilon     = values.next_double();
          params[nparams].sigma       = values.next_double();
          params[nparams].mu[0]         = values.next_double();
          params[nparams].mu[1]         = values.next_double();
          params[nparams].nu[0]        = values.next_double();
          params[nparams].nu[1]        = values.next_double();
          params[nparams].xi[0]         = values.next_double();
          params[nparams].xi[1]        = values.next_double();
          params[nparams].z0[0]       = values.next_double();
          params[nparams].z0[1]         = values.next_double();
          params[nparams].aij[0]        = values.next_double();
          params[nparams].aij[1]       = values.next_double();
          params[nparams].aik[0]       = values.next_double();
          params[nparams].aik[1]       = values.next_double();
          params[nparams].gammaij[0]   = values.next_double();
          params[nparams].gammaij[1]  = values.next_double();
          params[nparams].gammaik[0]    = values.next_double();
          params[nparams].gammaik[1]  = values.next_double();
          params[nparams].cos0[0]  = values.next_double();
          params[nparams].cos0[1]  = values.next_double();
          params[nparams].alpha[0]      = values.next_double();
          params[nparams].alpha[1]      = values.next_double();
          params[nparams].biga        = values.next_double();
          params[nparams].bigb        = values.next_double();
          params[nparams].powerp      = values.next_double();
          params[nparams].powerq      = values.next_double();
          params[nparams].cij         = values.next_double();
          params[nparams].tol         = values.next_double();

          if (unit_convert) {
            params[nparams].epsilon *= conversion_factor;
          }

          // check the sign of the parameters
          if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 ||
              params[nparams].mu[0] < 0.0 || params[nparams].mu[1] < 0.0 ||
              params[nparams].nu[0] < 0.0 || params[nparams].nu[1] < 0.0 ||
              params[nparams].xi[0] < 0.0 || params[nparams].xi[1] < 0.0 ||
              params[nparams].z0[0] < 0.0 || params[nparams].z0[1] < 0.0 ||
              params[nparams].aij[0] < 0.0 || params[nparams].aij[1] < 0.0 ||
              params[nparams].aik[0] < 0.0 || params[nparams].aik[1] < 0.0 ||
              params[nparams].gammaij[0] < 0.0 || params[nparams].gammaij[1] < 0.0 ||
              params[nparams].gammaik[0] < 0.0 || params[nparams].gammaik[1] < 0.0 ||
              params[nparams].cos0[0] > 0.0 || params[nparams].cos0[1] > 0.0 ||
              params[nparams].alpha[1] < 0.0 ||
              params[nparams].bigb < 0.0 ||
              params[nparams].powerp < 0.0 ||
              params[nparams].cij < 0.0
              ){
            error->one(FLERR,"Illegal ESW parameter");
          }

          nparams++;

        } else if(param_type == 1) {

          // Load soft paramters
          // load i-j pair
          for (ielement = 0; ielement < nelements; ielement++)
            if (iname == elements[ielement]) break;
          if (ielement == nelements) continue;
          for (jelement = 0; jelement < nelements; jelement++)
            if (jname == elements[jelement]) break;
          if (jelement == nelements) continue;
          
          // Dynamically expand the memory area of ​​the soft array 
          //    according to the number of loaded soft parameter sets
          if (nsofts == maxsofts) {
            maxsofts += DELTA;
            soft = (Softparam *) memory->srealloc(soft,
                                                  maxsofts * sizeof(Softparam),
                                                  "pair:soft");
            memset(soft + nsofts, 0, DELTA*sizeof(Softparam));
          }

          // Load 13 tokens ("SOFT" is skipped)
          soft[nsofts].ielement = ielement;
          soft[nsofts].jelement = jelement;
          soft[nsofts].p1 = values.next_double();
          soft[nsofts].p2 = values.next_double();
          soft[nsofts].p3 = values.next_double();
          soft[nsofts].p4 = values.next_double();
          soft[nsofts].p5 = values.next_double();
          soft[nsofts].p6 = values.next_double();
          soft[nsofts].p7 = values.next_double();
          soft[nsofts].p8 = values.next_double();
          soft[nsofts].p9 = values.next_double();
          soft[nsofts].bigr = values.next_double();
          soft[nsofts].bigd = values.next_double();

          nsofts++;

          // check the sign of the parameters
          if (soft[nsofts].p1 < 0.0 || soft[nsofts].p2 < 0.0 ||
              soft[nsofts].p3 < 0.0 || soft[nsofts].p4 < 0.0 ||
              soft[nsofts].p5 < 0.0 || soft[nsofts].p6 < 0.0 ||
              soft[nsofts].p7 < 0.0 || soft[nsofts].p8 < 0.0 ||
              soft[nsofts].p9 < 0.0 ||
              soft[nsofts].bigd < 0.0 || soft[nsofts].bigr < 0.0
              ){
            error->one(FLERR,"Illegal ESW soft parameter");
          }
        }   
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);
  MPI_Bcast(&nsofts, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxsofts, 1, MPI_INT, 0, world);

  if (comm->me != 0){
    // allocate memory in preparation for reception.
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
    soft = (Softparam *) memory->srealloc(soft, maxsofts*sizeof(Softparam), "pair:soft");
  }

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
  MPI_Bcast(soft, maxsofts*sizeof(Softparam), MPI_BYTE, 0, world);
}

void PairWYTO::display_params()
{
  /* display loaded parameters for Si-O systems*/
  int i;
  int Si;

  if (nelements == 2){
    int O;
    for (i = 0; i < nelements; i++){
      if (strcmp(elements[i], "Si") == 0){
        Si = i;
      } else if (strcmp(elements[i], "O") == 0){
        O = i;
      }
    }
    
    int SiSi, SiO, OSi, OO;
    int SiSiSi, SiSiO, OSiSi, SiOSi, OSiO; 

    // ij = elem3param[i][j][j]
    SiSi  = elem3param[Si][Si][Si];
    SiO   = elem3param[Si][O][O];
    OSi   = elem3param[O][Si][Si];
    OO    = elem3param[O][O][O];

    // jik = elem3param[i][j][k]
    SiSiSi = elem3param[Si][Si][Si];
    SiSiO = elem3param[Si][Si][O];
    OSiSi = elem3param[Si][O][Si];
    SiOSi = elem3param[O][Si][Si];
    OSiO = elem3param[Si][O][O];

    double A_SiSi           = params[SiSi].biga;     
    double B_SiSi           = params[SiSi].bigb;
    double p_SiSi           = params[SiSi].powerp;
    double q_SiSi           = params[SiSi].powerq;
    double a_SiSi           = params[SiSi].cij;
    double A_SiO            = params[SiO].biga;
    double B_SiO            = params[SiO].bigb;
    double p_SiO            = params[SiO].powerp;
    double q_SiO            = params[SiO].powerq;
    double a_SiO            = params[SiO].cij;
    double A_OO             = params[OO].biga;
    double B_OO             = params[OO].bigb;
    double q_OO             = params[OO].powerq;
    double a_OO             = params[OO].cij;
    double mu1_SiSiSi       = params[SiSiSi].mu[0];
    double gam1_SiSi_SiSiSi = params[SiSiSi].gammaij[0]; // = params[SiSiSi].gammaik[0] (1);

    double a1_SiSi_SiSiSi   = params[SiSiSi].aij[0];
    double cos01_SiSiSi     = params[SiSiSi].cos0[0];
    double alp1_SiSiSi      = params[SiSiSi].alpha[0];
    double mu2_SiSiSi       = params[SiSiSi].mu[1];
    double nu2_SiSiSi       = params[SiSiSi].nu[1];
    double xi2_SiSiSi       = params[SiSiSi].xi[1];
    double gam2_SiSi_SiSiSi = params[SiSiSi].gammaij[1]; // = params[SiSiSi].gammaik[1] (2);
    double z02_SiSiSi       = params[SiSiSi].z0[1];
    double a2_SiSi_SiSiSi   = params[SiSiSi].aij[1]; // = params[SiSiSi].aik[1] (3);
    double cos02_SiSiSi     = params[SiSiSi].cos0[1];
    double alp2_SiSiSi      = params[SiSiSi].alpha[1];
    double mu1_SiSiO        = params[SiSiO].mu[0];
    double nu1_SiSiO        = params[SiSiO].nu[0];
    double xi1_SiSiO        = params[SiSiO].xi[0];
    double z01_SiSiO        = params[SiSiO].z0[0];
    double gam1_SiSi_SiSiO  = params[SiSiO].gammaij[0]; // = params[OSiSi].gammaik[0] (4)

    double gam1_SiO_SiSiO  = params[SiSiO].gammaik[0]; // = params[OSiSi].gammaij[0] (4)
    double a1_SiSi_SiSiO   = params[SiSiO].aij[0]; // = params[OSiSi].aik[0] (5)
    double a1_SiO_SiSiO    = params[SiSiO].aik[0]; // = params[OSiSi].aij[0] (5)
    double cos01_SiSiO     = params[SiSiO].cos0[0];
    double alp1_SiSiO      = params[SiSiO].alpha[0];
    double mu2_SiSiO       = params[SiSiO].mu[1];
    double gam2_SiSi_SiSiO = params[SiSiO].gammaij[1]; // = params[OSiSi].gammaik[1] (6)
    double gam2_SiO_SiSiO  = params[SiSiO].gammaik[1]; // = params[OSiSi].gammaij[1] (6)
    double a2_SiSi_SiSiO   = params[SiSiO].aij[1]; // = params[OSiSi].aik[1] (7)
    double a2_SiO_SiSiO    = params[SiSiO].aik[1]; // = params[OSiSi].aij[1] (7)
    double cos02_SiSiO     = params[SiSiO].cos0[1];
    double alp2_SiSiO      = params[SiSiO].alpha[1];
    double mu1_SiOSi       = params[SiOSi].mu[0];
    double gam1_OSi_SiOSi  = params[SiOSi].gammaij[0]; // = params[SiOSi].gammaik[0] (8)
    double a1_OSi_SiOSi    = params[SiOSi].aij[0]; // params[SiOSi].aik[0]  (9)
    double cos01_SiOSi     = params[SiOSi].cos0[0];

    double alp1_SiOSi     = params[SiOSi].alpha[0];
    double mu1_OSiO       = params[OSiO].mu[0];
    double gam1_SiO_OSiO = params[OSiO].gammaij[0]; // = params[OSiO].gammaik[0] (10)
    double a1_SiO_OSiO = params[OSiO].aij[0]; // = params[OSiO].aik[0] (11)
    double cos01_OSiO = params[OSiO].cos0[0];

    SiO = elem2soft[Si][O];

    double p1 = soft[SiO].p1;
    double p2 = soft[SiO].p2;
    double p3 = soft[SiO].p3;
    double p4 = soft[SiO].p4;
    double p5 = soft[SiO].p5;
    double p6 = soft[SiO].p6;
    double p7 = soft[SiO].p7;
    double p8 = soft[SiO].p8;
    double p9 = soft[SiO].p9;
    double R  = soft[SiO].bigr;
    double D  = soft[SiO].bigd;

    // Check parameter consistency
    if (
      params[SiSiSi].gammaij[0] != params[SiSiSi].gammaik[0] ||
      params[SiSiSi].gammaij[1] != params[SiSiSi].gammaik[1] ||
      params[SiSiSi].aij[1] != params[SiSiSi].aik[1] ||
      params[SiSiO].gammaij[0] != params[OSiSi].gammaik[0] ||
      params[SiSiO].aij[0] != params[OSiSi].aik[0] ||
      params[SiSiO].gammaij[1] != params[OSiSi].gammaik[1] ||
      params[SiSiO].aij[1] != params[OSiSi].aik[1] ||
      params[SiOSi].gammaij[0] != params[SiOSi].gammaik[0] ||
      params[SiOSi].aij[0] != params[SiOSi].aik[0] ||
      params[OSiO].gammaij[0] != params[OSiO].gammaik[0] ||
      params[OSiO].aij[0] != params[OSiO].aik[0]
    ){
          ROOT_LOG(comm->me,
                    params[SiSiSi].gammaij[0] << "," << params[SiSiSi].gammaik[0] << "," <<
                    params[SiSiSi].gammaij[1] << "," << params[SiSiSi].gammaik[1] << "," <<
                    params[SiSiSi].aij[1] << "," << params[SiSiSi].aik[1] << "," <<
                    params[SiSiO].gammaij[0] << "," << params[OSiSi].gammaik[0] << "," <<
                    params[SiSiO].aij[0] << "," << params[OSiSi].aik[0] << "," <<
                    params[SiSiO].gammaij[1] << "," << params[OSiSi].gammaik[1] << "," <<
                    params[SiSiO].aij[1] << "," << params[OSiSi].aik[1] << "," <<
                    params[SiOSi].gammaij[0] << "," << params[SiOSi].gammaik[0] << "," <<
                    params[SiOSi].aij[0] << "," << params[SiOSi].aik[0] << "," <<
                    params[OSiO].gammaij[0] << "," << params[OSiO].gammaik[0] << "," <<
                    params[OSiO].aij[0] << "," << params[OSiO].aik[0]);
          error->all(FLERR, "Inconsistent parameters.");
    }

    if (comm->me == 0){
      printf("## Paramters ##\n");
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "A_SiSi",           A_SiSi,           "a_SiSi_SiSiSi",    a1_SiSi_SiSiSi,   "gam1_SiO_SiSiO", gam1_SiO_SiSiO, "alp1_SiOSi", alp1_SiOSi);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "B_SiSi",           B_SiSi,           "cos01_SiSiSi",     cos01_SiSiSi,     "a1_SiSi_SiSiO",  a1_SiSi_SiSiO,  "mu1_OSiO", mu1_OSiO);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "p_SiSi",           p_SiSi,           "alp1_SiSiSi",      alp1_SiSiSi,      "a1_SiO_SiSiO",   a1_SiO_SiSiO,   "gam1_SiO_OSiO", gam1_SiO_OSiO);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "q_SiSi",           q_SiSi,           "mu2_SiSiSi",       mu2_SiSiSi,       "cos01_SiSiO",    cos01_SiSiO,    "a1_SiO_OSiO", a1_SiO_OSiO);      
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "a_SiSi",           a_SiSi,           "nu2_SiSiSi",       nu2_SiSiSi,       "alp1_SiSiO",     alp1_SiSiO,     "cos01_OSiO", cos01_OSiO);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "A_SiO",            A_SiO,            "xi2_SiSiSi",       xi2_SiSiSi,       "mu2_SiSiO",      mu2_SiSiO,      "p1", p1);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "B_SiO",            B_SiO,            "gam2_SiSi_SiSiSi", gam2_SiSi_SiSiSi, "a1_SiO_SiSiO",   gam2_SiSi_SiSiO,"p2", p2);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "p_SiO",            p_SiO,            "z02_SiSiSi",       z02_SiSiSi,       "gam2_SiO_SiSiO", gam2_SiO_SiSiO, "p3", p3);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "q_SiO",            q_SiO,            "a2_SiSi_SiSiSi",   a2_SiSi_SiSiSi,   "a1_SiO_SiSiO",   a2_SiSi_SiSiO,  "p4", p4);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "a_SiO",            a_SiO,            "cos02_SiSiSi",     cos02_SiSiSi,     "a2_SiSi_SiSiO",  a2_SiO_SiSiO,   "p5", p5);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "A_OO",             A_OO,             "alp2_SiSiSi",      alp2_SiSiSi,      "cos02_SiSiO",    cos02_SiSiO,    "p6", p6);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "B_OO",             B_OO,             "mu1_SiSiO",        mu1_SiSiO,        "alp2_SiSiO",     alp2_SiSiO,     "p7", p7);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "q_OO",             q_OO,             "nu1_SiSiO",        nu1_SiSiO,        "mu1_SiOSi",      mu1_SiOSi,      "p8", p8);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "a_OO",             a_OO,             "xi1_SiSiO",        xi1_SiSiO,        "gam1_OSi_SiOSi", gam1_OSi_SiOSi, "p9", p9);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "mu1_SiSiSi",       mu1_SiSiSi,       "z01_SiSiO",        z01_SiSiO,        "a1_OSi_SiOSi",   a1_OSi_SiOSi,   "R", R);
      printf("%-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f   %-18s%-9.6f\n", 
              "gam1_SiSi_SiSiSi", gam1_SiSi_SiSiSi, "gam1_SiSi_SiSiO",  gam1_SiSi_SiSiO,  "cos01_SiOSi",   cos01_SiOSi,    "D", D);
      

      printf("## ##\n");
    }
  } else if (nelements == 1) {

  } 
}

void PairWYTO::setup_params()
{
  //LOG("Enter setup_params");

  /* set cutoff distance, elem2soft, elem3param */

  int i,j,k,m,n,o,p;

  // set elem3param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");
  memory->destroy(elem2soft);
  memory->create(elem2soft, nelements, nelements, "pair:elem2soft");

  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {

      // set elem2soft and softflag
      o = -1;
      for (p = 0; p < nsofts; p++) {
        if (i == soft[p].ielement && j == soft[p].jelement) {

          // When the atom type of the p-th soft parameter set are i,j
          softflag[i][j] = 1;   // flag i-j pair
          softflag[j][i] = 1;   // j-i pair as well
          if (o >= 0) 
            // A state in which the soft parameter set with atomic types i and j has already been loaded.
            error->all(FLERR,"  Potential file has duplicate entry");
          o = p;
          ROOT_PRINT(comm->me, "  Soft parameter for "<< elements[i] << "," << elements[j] << " found");
        }
      }

      elem2soft[i][j] = o;  // associate i-j pair with soft parameter index
      elem2soft[j][i] = o;  // j-i pair as well

      //LOG("Set elem2soft for a pair");
      
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
            ROOT_PRINT(comm->me, "  Parameter for "
                        << elements[i] << "," 
                        << elements[j] << "," 
                        << elements[k] << " found");
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem3param[i][j][k] = n;
      }
    } 
  }

  /*
  if (params[m].tol > 0.0){
      LOG("tol:" << params[m].tol);
      error->all(FLERR,"Potential not currently set to accept tol values");
  }
  */

  // set cutoff distance 
  double rtmp1, rtmp2, rtmp3, rtmp4, rtmp5;
  for (m = 0; m < nparams; m++) {

    // calculate each cutoff distance[Å]
    params[m].cutpair = params[m].sigma * params[m].cij;
    params[m].cutij[0] = params[m].sigma * params[m].aij[0];
    params[m].cutij[1] = params[m].sigma * params[m].aij[1];
    params[m].cutik[0] = params[m].sigma * params[m].aik[0];
    params[m].cutik[1] = params[m].sigma * params[m].aik[1];

    // store each cutoff distnace in buffer
    rtmp1 = params[m].cutpair;
    rtmp2 = params[m].cutij[0];
    rtmp3 = params[m].cutij[1];
    rtmp4 = params[m].cutik[0];
    rtmp5 = params[m].cutik[1];
      
    // calculate squared cutoff distnace[Å^2]
    params[m].cutpairsq = rtmp1 * rtmp1;
    params[m].cutijsq[0] = rtmp2 * rtmp2;
    params[m].cutijsq[1] = rtmp3 * rtmp3;
    params[m].cutiksq[0] = rtmp4 * rtmp4;
    params[m].cutiksq[1] = rtmp5 * rtmp5;

    ROOT_LOG(comm->me, "Square cutoff length ...");
    ROOT_LOG(comm->me, "m:"<< m);
    ROOT_LOG(comm->me, "  pairsq:" << params[m].cutpairsq);
    ROOT_LOG(comm->me, "  cutijsq1:" << params[m].cutijsq[0]);
    ROOT_LOG(comm->me, "  cutijsq2:" << params[m].cutijsq[1]);
    ROOT_LOG(comm->me, "  cutiksq1:" << params[m].cutiksq[0]);
    ROOT_LOG(comm->me, "  cutiksq2:" << params[m].cutiksq[1]);
  }

  // get maximum cutoff distance
  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp1 = sqrt(params[m].cutpairsq);
    rtmp2 = sqrt(params[m].cutijsq[0]);
    rtmp3 = sqrt(params[m].cutijsq[1]);
    rtmp4 = sqrt(params[m].cutiksq[0]);
    rtmp5 = sqrt(params[m].cutiksq[1]);

    if (rtmp1 > cutmax) cutmax = rtmp1;
    if (rtmp2 > cutmax) cutmax = rtmp2;
    if (rtmp3 > cutmax) cutmax = rtmp3;
    if (rtmp4 > cutmax) cutmax = rtmp4;
    if (rtmp5 > cutmax) cutmax = rtmp5;
  }
  //LOG("Exit setup_params");
}

void PairWYTO::count_z(int nall, int nmax, int inum, int *ilist, double **x, int *tag, int *type, int *numneigh, int **firstneigh){
  /* # calculate coordination number  */
  //  - update numcoord

  int i, j, ii, jj, itag, jtag, itype, jtype, jnum, ijparam, ijsoft;
  Softparam *softij;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz, rsq, r;
  int *jlist;

  // initialize numcoord
  for (i = 0; i < nall; i++) numcoord[i] = 0.0;

  // calculate coordination number for responsible subdmain
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];            // get index of target atom i
    itype = map[type[i]];     // get atom type of i 

    // store coordinate of atom i in buffer
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];  // get neighbor list of atom i
    jnum = numneigh[i];     // get length of neighbor list of atom i

    for (jj = 0; jj < jnum; jj++) {
      //ROOT_LOG(comm->me, "loop for jj");

      // get index of atom j
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = map[type[j]];  // get atom type of i 

      // count atom j whose atomic type is different from the i atom as coordinating atoms
      if (softflag[itype][jtype]) {
        //ROOT_LOG(comm->me, "softflag true");

        // calculate squared distance i-j[Å^2]
        delx = xtmp - x[j][0];                      
        dely = ytmp - x[j][1];                      
        delz = ztmp - x[j][2];                      
        rsq = delx*delx + dely*dely + delz*delz;

        // calculate SW distance i-j[σ=2.0951Å]
        ijparam = elem3param[itype][jtype][jtype];
        r = sqrt(rsq) / params[ijparam].sigma;

        // calculate coordination number of atom i using cutoff function
        numcoord[i] += fc(r, itype, jtype);
      }
    }
  }

  // receiving coordination number in ghost atom from neighbor processes
  comm->forward_comm(this);
}

int PairWYTO::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  /* function to send the coordination number to neighbor process */
  // Reference:pair_eam.cpp

  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = numcoord[j];
  }
  return m;
}

void PairWYTO::unpack_forward_comm(int n, int first, double *buf)
{
  /* function to reflect the received coordination number in the own process */
  // Reference:pair_eam.pp

  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) numcoord[i] = buf[m++];
}


double PairWYTO::gSi(double z, int itype, int jtype)
{
  double p1, p2, p3, p4, first, second;
  int ij = elem2soft[itype][jtype];
  p1 = soft[ij].p1;
  p2 = soft[ij].p2;
  p3 = soft[ij].p3;
  p4 = soft[ij].p4;
  if (z >= 4.0){
    return p4;
  } else {
    first = p1 * sqrt(z + p2) - p4;
    second = exp(p3 / (z - 4.0));
    return first * second + p4;
  }
}

double PairWYTO::gO(double z, int itype, int jtype)
{
  double p5, p6, p7, p8, p9, first, second;
  int ij = elem2soft[itype][jtype];

  p5 = soft[ij].p5;
  p6 = soft[ij].p6;
  p7 = soft[ij].p7;
  p8 = soft[ij].p8;
  p9 = soft[ij].p9;
  first = p5 / (exp((p6 - z) / p7) + 1.0);
  second = exp(p8 * (z - p9) * (z - p9));
  return first * second;
}

double PairWYTO::d_gSi(double z, int itype, int jtype){
  /* dgSi/dz */

  double p1, p2, p3, p4, tmp[7];
  int ij = elem2soft[itype][jtype];
  p1 = soft[ij].p1;
  p2 = soft[ij].p2;
  p3 = soft[ij].p3;
  p4 = soft[ij].p4; 
  if (z >= 4.0){
    return 0.0;
  } else {
    tmp[0] = p3 / (z - 4.0);
    tmp[1] = exp(tmp[0]);
    tmp[2] = sqrt(z + p2);
    tmp[3] = p1 / (2.0 * tmp[2]);
    tmp[4] = tmp[0] / (z - 4.0);
    tmp[5] = p1 * tmp[2] - p4;

    return tmp[1] * (tmp[3] - tmp[4] * tmp[5]);
  }
}

double PairWYTO::d_gO(double z, int itype, int jtype){
  /* dgO/dz */

  double p5, p6, p7, p8, p9, tmp[8];
  int ij = elem2soft[itype][jtype];
   
  p5 = soft[ij].p5;
  p6 = soft[ij].p6;
  p7 = soft[ij].p7;
  p8 = soft[ij].p8;
  p9 = soft[ij].p9;

  tmp[0] = (p6 - z) / p7;
  tmp[1] = exp(tmp[0]) + 1.0;
  tmp[2] = p8 * (z - p9);
  tmp[3] = tmp[2] * (z - p9);

  tmp[4] = p5 / (p7 * tmp[1] * tmp[1]);
  tmp[5] = exp(tmp[0] + tmp[3]);
  tmp[6] = 2.0 * p5 * tmp[2] / tmp[1];
  tmp[7] = exp(p8 * (z - p9) * (z - p9));

  return tmp[4] * tmp[5] + tmp[6] * tmp[7];

}


double PairWYTO::fc(double r, int itype, int jtype){
  /* cutof  function fc */

  double bigd, bigr;
  int ij = elem2soft[itype][jtype];

  bigd = soft[ij].bigd;
  bigr = soft[ij].bigr;

  if(r < bigr - bigd) return 1.0;
  else if((bigr - bigd <= r) && (r < bigr + bigd))
    return 1.0 - (r - bigr + bigd) / (2.0 * bigd) + sin(PIVAL * (r - bigr + bigd) / bigd)/(2.0 * PIVAL);
  else return 0.0;
}

double PairWYTO::d_fc(double r, int itype, int jtype){
  /* dfc/dr */

  double bigd, bigr;
  int ij = elem2soft[itype][jtype];

  bigd = soft[ij].bigd;
  bigr = soft[ij].bigr;

  if((bigr - bigd <= r) && (r < bigr + bigd)){
    return (0.5 / bigd) * cos(PIVAL / bigd * (r - bigr + bigd)) - 0.5 / bigd;
  } else {
    return 0.0;
  }
}

void PairWYTO::twobody_regular(
  int i, int j, int itype, int jtype,
  double dxij, double dyij, double dzij, 
  double xtmpi, double ytmpi, double ztmpi,
  double xtmpj, double ytmpj, double ztmpj,
  double **x, int *tag, int *type, int *numneigh, int **firstneigh,   
  Param *paramij, double rsq, double **force, int eflag, double &evdwl, 
  int nlocal, int newton_pair
){
  /* calculate twobody energy and regular force*/

  int k, kk, knum, ktype, ikparam, jkparam, ijsoft, iksoft, jksoft;
  int *klist;
  Param *paramik, *paramjk;
  Softparam *softij, *softik, *softjk;
  double tmp[7];
  double epsilon, A, B, p, q, a, sigma;
  double rij, rik, rjk, dxik, dyik, dzik, dxjk, dyjk, dzjk;
  std::string ielem, jelem, kelem;
  double gij, g1, g2, dg1, dg2, zi, zj;
  double fpair;
  const double ecoul = 0.0; // ecoul = 0.0 in this potential

  /* calculate energy */

  // Load SW paramters
  epsilon = paramij->epsilon;
  A = paramij->biga;
  B = paramij->bigb;
  p = paramij->powerp;
  q = paramij->powerq;
  a = paramij->cij;
  sigma = paramij->sigma;

  // calculate SW distance
  rij = sqrt(rsq) / sigma;
  
  // calculate SW twobody term
  tmp[0] = 1.0 / (rij - a);
  tmp[1] = exp(tmp[0]);                             // exp term
  tmp[2] = A * (B * pow(rij, - p) - pow(rij, -q));  // lj term
  tmp[3] = tmp[1] * tmp[2];                         // SW term

  // get coordination numbers
  zi = numcoord[i]; 
  zj = numcoord[j];

  // calculate gij
  if (strcmp(elements[itype], "Si") == 0  && 
      strcmp(elements[jtype], "O") == 0){
    //LOG("- i = Si, j = O");
    g1 = gSi(zi, itype, jtype);
    g2 = gO(zj, jtype, itype);
    gij = g1 * g2;

  } else if (strcmp(elements[itype], "O") == 0 && 
            strcmp(elements[jtype], "Si") == 0) {
    //LOG("- i = O, j = Si");
    g1 = gO(zi, itype, jtype);
    g2 = gSi(zj, jtype, itype);
    gij = g1 * g2;

  } else {
    //LOG("- i,j same");
    gij = 1.0;
  }
  
  // calculate energy
  evdwl = epsilon * gij * tmp[3];
  //ROOT_LOG(comm->me, "evdwl2=" << epsilon << "*" << gij << "*" << tmp[3] << "=" << evdwl);
  if (comm->me == 0 && !std::isfinite(evdwl)) exit(1);

  /* calculate regular force */

  // calculate the derivative of the SW twobody term
  // derivative of lj term 
  tmp[4] = A * (- p * B * pow(rij, - p - 1.0) + q * pow(rij, - q - 1.0));
  // deriative of exp term
  tmp[5] = - tmp[0] * tmp[0] * tmp[1];
  // derivative of SW twobody term
  tmp[6] = tmp[4] * tmp[1] + tmp[2] * tmp[5];
  
  // calculate force
  fpair = epsilon * gij * tmp[6] / (rij * sigma * sigma);

  // distribute force
  force[i][0] += - fpair * dxij;
  force[i][1] += - fpair * dyij;
  force[i][2] += - fpair * dzij;
  force[j][0] += fpair * dxij;
  force[j][1] += fpair * dyij;
  force[j][2] += fpair * dzij;

  // distribute energy and calculate virial(i-j)
  if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, ecoul, - fpair, dxij, dyij, dzij);
}

void PairWYTO::twobody_external(
  int i, int j, int itype, int jtype,
  double dxij, double dyij, double dzij, 
  double xtmpi, double ytmpi, double ztmpi,
  double **x, int *tag, int *type, int *numneigh, int **firstneigh,   
  Param *paramij, double rsq, double **force, int eflag, double &evdwl, 
  int nlocal, int newton_pair
){

  /* calculate twobody external force i-k*/
  // - for full neighbor list

  int k, kk, knum, ktype, ikparam;
  int *klist;
  Param *paramik;
  double tmp[4];
  double epsilon, A, B, p, q, a, sigma;
  double R, D;
  double rij, rik, dxik, dyik, dzik;
  double g2, dg1, zi, zj;
  double fpair;
  const double ecoul = 0.0; // ecoul = 0.0 in this pair_style
  
  // Load SW parameters
  epsilon = paramij->epsilon;
  A = paramij->biga;
  B = paramij->bigb;
  p = paramij->powerp;
  q = paramij->powerq;
  a = paramij->cij;
  sigma = paramij->sigma;

  // calculate SW distance rij
  rij = sqrt(rsq) / sigma;
  
  // calculate SW two-body term
  tmp[0] = 1.0 / (rij - a);
  tmp[1] = exp(tmp[0]);                             // exp term
  tmp[2] = A * (B * pow(rij, - p) - pow(rij, -q));  // lj term
  tmp[3] = tmp[1] * tmp[2];                         // SW term

  // get coordination numbers zi, zj
  zi = numcoord[i]; 
  zj = numcoord[j];
  
  // calculate g2(zj),dg1(zi)/dz
  if (strcmp(elements[itype], "Si") == 0  && 
      strcmp(elements[jtype], "O") == 0){
    g2 = gO(zj, jtype, itype);
    dg1 = d_gSi(zi, itype, jtype);
  } else if (strcmp(elements[itype], "O") == 0 && 
            strcmp(elements[jtype], "Si") == 0) {
    g2 = gSi(zj, jtype, itype);
    dg1 = d_gO(zi, itype, jtype);
  } else {
    g2 = 0.0;
    dg1 = 0.0;
  }

  // calculate external force (i-k)

  klist = firstneigh[i];
  knum = numneigh[i];
  for (kk = 0; kk < knum; kk++) {
    k = klist[kk];
    k &= NEIGHMASK;
    ktype = map[type[k]];
    if (softflag[itype][ktype]) {
      // i and k are different atoms

      // calculate rik
      dxik = xtmpi - x[k][0];
      dyik = ytmpi - x[k][1];
      dzik = ztmpi - x[k][2];
      rik = sqrt(dxik * dxik + dyik * dyik + dzik * dzik) / sigma;

      // calculate force
      fpair = epsilon * g2 * tmp[3] * dg1 * d_fc(rik, itype, ktype) / rik; // [eV*(σ/A)^2 /atom]
      fpair /= sigma * sigma;                                              // [eV*/A^2 /atom]

      // distribute force
      force[i][0] += - fpair * dxik; // [eV*/A /atom]
      force[i][1] += - fpair * dyik;
      force[i][2] += - fpair * dzik;
      force[k][0] += fpair * dxik;
      force[k][1] += fpair * dyik;
      force[k][2] += fpair * dzik;
      
      // Calcurate virial (i-k)
      double vecrik[3] = {dxik, dyik, dzik};
      if(vflag_either) v_tally2(i, k, - fpair, vecrik);
    }
  }
}

double PairWYTO::lambda(double z, double mu, double nu, double xi, double z0){
  return mu * (1.0 + nu * exp(-xi * (z - z0) * (z - z0)));
}

double PairWYTO::d_lambda(double z, double mu, double nu, double xi, double z0){
  /* dλ/dz */

  double tmp[4];
  tmp[0] = -2.0 * mu * nu;
  tmp[1] = z - z0;
  tmp[2] = exp(- xi * (tmp[1] * tmp[1]));
  tmp[3] = xi * tmp[1];
  return tmp[0] * tmp[2] * tmp[3];
}

void PairWYTO::threebody(
  int i, int j, int k, int n, 
  int itype, int jtype, int ktype,
  double dxij, double dyij, double dzij, 
  double dxik, double dyik, double dzik,
  double dxjk, double dyjk, double dzjk,
  double **x, int *tag, int *type, int *numneigh, int **firstneigh,
  Param *paramijk, double rsq1, double rsq2, double **force, 
  int eflag, double &evdwl, int nlocal, int newton_pair
){
  /* 
  - Calculate the energy and force of three-body terms
    - Calculate Λn(i,j,k)Θn(θjik)
  */

  double rij, rik, rjk;
  double epsilon;
  double zi;
  double mu, nu, xi, z0;
  double gammaij, gammaik;
  double aij, aik;
  double alpha, cos0, cosjik;
  double sigma;
  double tmp[5];
  double lam;
  double biglam, biglam_p, bigthe;
  double d_biglam, d_biglam_p, d_bigthe, dcosjik, d_bl_bt;
  double fij, fik, fjk;
  
  // for ev_tally
  double fj[3] = {0.0, 0.0, 0.0}, fk[3] = {0.0, 0.0, 0.0};
  double vecij[3] = {- dxij, - dyij, - dzij};
  double vecik[3] = {- dxik, - dyik, - dzik};
  const double ecoul = 0.0; 
  
  // for external force calculation
  int l, ll, lnum, ilparam, ltype;
  int *llist;
  double xtmpi, ytmpi, ztmpi, ril;
  double dxil, dyil, dzil;
  double dlam, dfc, fil;
  double vecil[3];

  // Parameter loading
  epsilon = paramijk->epsilon;
  mu = paramijk->mu[n];
  nu = paramijk->nu[n];
  xi = paramijk->xi[n];
  z0 = paramijk->z0[n];
  gammaij = paramijk->gammaij[n];
  gammaik = paramijk->gammaik[n];
  aij = paramijk->aij[n];
  aik = paramijk->aik[n];
  alpha = paramijk->alpha[n];
  cos0 = paramijk->cos0[n];
  sigma = paramijk->sigma;

  // Calculate SW distance
  rij = sqrt(rsq1) / sigma;
  rik = sqrt(rsq2) / sigma;
  rjk = sqrt(dxjk*dxjk + dyjk*dyjk + dzjk*dzjk) / sigma;

  // Calculate Λ term
  zi = numcoord[i];
  tmp[0] = 1.0 / (rij - aij);
  tmp[1] = gammaij * tmp[0];
  tmp[2] = 1.0 / (rik - aik);
  tmp[3] = gammaik * tmp[2];
  biglam_p = exp(tmp[1] + tmp[3]);  // Λ' tern(exp term)
  lam = lambda(zi, mu, nu, xi, z0); // λ term
  biglam = lam * biglam_p;          // Λ term

  // Calculate cosΘjik
  cosjik = (dxij * dxik + dyij * dyik + dzij * dzik) / (rij * rik) / sigma / sigma;

  // Calculate Θ term
  tmp[4] = cosjik - cos0;
  bigthe = tmp[4] * tmp[4] + alpha * tmp[4] * tmp[4] * tmp[4];

  // Energy calculation
  evdwl = epsilon * biglam * bigthe;
  if (comm->me == 0 && !std::isfinite(evdwl)) exit(1);

  /* Regular force calculation */
  
  // Force i-j

  // dΛ'/drij
  d_biglam_p = -gammaij * biglam_p * tmp[0] * tmp[0];
  // dcosθjik/drij
  dcosjik = (rij - rik * cosjik) / rij / rik;
  // dΘ/drij
  d_bigthe = (2.0 * tmp[4] + 3.0 * alpha * tmp[4] * tmp[4]) * dcosjik;
  // d(Λ'Θ)/drij
  d_bl_bt = d_biglam_p * bigthe + biglam_p * d_bigthe;

  // Force distribution
  fij = epsilon * lam * d_bl_bt / rij;  // [eV*(σ/A)^2 /atom]
  fij /= sigma * sigma;                 // [eV/A^2 /atom]
  if (comm->me == 0 && !std::isfinite(fij)) exit(1);
  force[i][0] += - fij * dxij; // [eV/A /atom]
  force[i][1] += - fij * dyij;
  force[i][2] += - fij * dzij;
  force[j][0] += fij * dxij;
  force[j][1] += fij * dyij;
  force[j][2] += fij * dzij;

  // update fj(for ev_tally3)
  fj[0] += fij * dxij;
  fj[1] += fij * dyij;
  fj[2] += fij * dzij;

  // Force i-k 

  // dΛ'/drik
  d_biglam_p = -gammaik * biglam_p *  tmp[2] * tmp[2];
  // dcosθjik/drik
  dcosjik = (rik - rij * cosjik) / rij / rik;
  // dΘ/drik
  d_bigthe = (2.0 * tmp[4] + 3.0 * alpha * tmp[4] * tmp[4]) * dcosjik;
  // d(Λ'Θ)/drik
  d_bl_bt = d_biglam_p * bigthe + biglam_p * d_bigthe;

  // force distribution
  fik = epsilon * lam * d_bl_bt / rik;  // [eV*(σ/Å)^2 /atom]
  fik /= sigma * sigma;                 // [eV/Å^2 /atom]
  if (comm->me == 0 && !std::isfinite(fik)) exit(1);
  force[i][0] += - fik * dxik;  // [eV/Å /atom]
  force[i][1] += - fik * dyik;
  force[i][2] += - fik * dzik;
  force[k][0] += fik * dxik;
  force[k][1] += fik * dyik;
  force[k][2] += fik * dzik;
  
  // update fk (for ev_tally3)
  fk[0] += fik * dxik;
  fk[1] += fik * dyik;
  fk[2] += fik * dzik;

  // distribute energy and calculate virial(i-j, i-k)
  if(evflag) ev_tally3(i, j, k, evdwl, ecoul, fj, fk, vecij, vecik);

  // Calculate force (j-k)
  
  // dΛ'/drjk
  d_biglam_p = 0;
  // dcosθjik/drjk
  dcosjik = -rjk / rij / rik;
  // dΘ/drjk
  d_bigthe = (2.0 * tmp[4] + 3.0 * alpha * tmp[4] * tmp[4]) * dcosjik;
  // d(Λ'Θ)/drjk
  d_bl_bt = d_biglam_p * bigthe + biglam_p * d_bigthe;

  // distribute force
  fjk = epsilon * lam * d_bl_bt / rjk; // [eV*(σ/Å)^2 /atom]
  fjk /= sigma * sigma;                // [eV/Å^2 /atom]
  if (comm->me == 0 && !std::isfinite(fjk)) exit(1);
  force[j][0] += - fjk * dxjk; // [eV/Å /atom]
  force[j][1] += - fjk * dyjk;
  force[j][2] += - fjk * dzjk;
  force[k][0] += fjk * dxjk;
  force[k][1] += fjk * dyjk;
  force[k][2] += fjk * dzjk;

  // Virial calculation
  double vecjk[3] = {dxjk, dyjk, dzjk}; // dxjk = x[j][0] - x[k][0]
  if (vflag_either) v_tally2(j, k, - fjk, vecjk);

  /* External force calculation */

  xtmpi = x[i][0];
  ytmpi = x[i][1];
  ztmpi = x[i][2];

  // loop for neighbor of i
  lnum = numneigh[i];   
  llist = firstneigh[i];
  for (ll = 0; ll < lnum; ll++) {
    l = llist[ll];
    l &= NEIGHMASK;
    ltype = map[type[l]];
    if (softflag[itype][ltype]){
      
      // calcurate ril
      ilparam = elem3param[itype][ltype][ltype];
      sigma = params[ilparam].sigma;
      dxil = xtmpi - x[l][0]; // [Å]
      dyil = ytmpi - x[l][1];
      dzil = ztmpi - x[l][2];
      ril = sqrt(dxil * dxil + dyil * dyil + dzil * dzil) / sigma; // [Å/σ]
      
      // dλ/dzi
      dlam = d_lambda(zi, mu, nu, xi, z0);
      
      // dfc/ril
      dfc = d_fc(ril, itype, ltype);
      
      // distribute　force
      fil = epsilon * biglam_p * bigthe * dlam * dfc / ril; // [eV*(σ/Å)^2 /atom]
      fil /= sigma * sigma;                                 // [eV/Å^2 /atom]
      if (comm->me == 0 && !std::isfinite(fil)) exit(1);
      force[i][0] += - fil *dxil;                           // [eV/Å /atom]
      force[i][1] += - fil *dyil;
      force[i][2] += - fil *dzil;
      force[l][0] += fil *dxil;
      force[l][1] += fil *dyil;
      force[l][2] += fil *dzil;

      // calculate virial(i-l)
      vecil[0] = dxil;
      vecil[1] = dyil;
      vecil[2] = dzil;
      if(vflag_either) v_tally2(i, l, - fil, vecil);
    }
  }
}
