/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_MDTYPES_INPUTREC_H
#define GMX_MDTYPES_INPUTREC_H

#include <cstdio>

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#define EGP_EXCL  (1<<0)
#define EGP_TABLE (1<<1)

struct pull_params_t;

namespace gmx
{
class Awh;
struct AwhParams;
class KeyValueTreeObject;
}

typedef struct t_grpopts {
    int       ngtc;           /* # T-Coupl groups                        */
    int       nhchainlength;  /* # of nose-hoover chains per group       */
    int       ngacc;          /* # Accelerate groups                     */
    int       ngfrz;          /* # Freeze groups                         */
    int       ngener;         /* # Ener groups			    */
    Real     *nrdf;           /* Nr of degrees of freedom in a group	    */
    Real     *ref_t;          /* Coupling temperature	per group   */
    int      *annealing;      /* No/simple/periodic SA for each group    */
    int      *anneal_npoints; /* Number of annealing time points per grp */
    Real    **anneal_time;    /* For ea. group: Time points              */
    Real    **anneal_temp;    /* For ea. grp: Temperature at these times */
                              /* Final temp after all intervals is ref_t */
    Real     *tau_t;          /* Tau coupling time              */
    rvec     *acc;            /* Acceleration per group		    */
    ivec     *nFreeze;        /* Freeze the group in each direction ?    */
    int      *egp_flags;      /* Exclusions/tables of energy group pairs */

    /* QMMM stuff */
    int          ngQM;         /* nr of QM groups                              */
    int         *QMmethod;     /* Level of theory in the QM calculation        */
    int         *QMbasis;      /* Basisset in the QM calculation               */
    int         *QMcharge;     /* Total charge in the QM region                */
    int         *QMmult;       /* Spin multiplicicty in the QM region          */
    gmx_bool    *bSH;          /* surface hopping (diabatic hop only)          */
    int         *CASorbitals;  /* number of orbiatls in the active space       */
    int         *CASelectrons; /* number of electrons in the active space      */
    Real        *SAon;         /* at which gap (A.U.) the SA is switched on    */
    Real        *SAoff;
    int         *SAsteps;      /* in how many steps SA goes from 1-1 to 0.5-0.5*/
} t_grpopts;

typedef struct t_simtemp {
    int   eSimTempScale; /* simulated temperature scaling; linear or exponential */
    Real  simtemp_low;   /* the low temperature for simulated tempering  */
    Real  simtemp_high;  /* the high temperature for simulated tempering */
    Real *temperatures;  /* the range of temperatures used for simulated tempering */
} t_simtemp;

typedef struct t_lambda {
    int    nstdhdl;                 /* The frequency for calculating dhdl           */
    double init_lambda;             /* fractional value of lambda (usually will use
                                       init_fep_state, this will only be for slow growth,
                                       and for legacy free energy code. Only has a
                                       valid value if positive)   */
    int      init_fep_state;        /* the initial number of the state                   */
    double   delta_lambda;          /* change of lambda per time step (fraction of (0.1) */
    int      edHdLPrintEnergy;      /* print no, total or potential energies in dhdl    */
    int      n_lambda;              /* The number of foreign lambda points               */
    double **all_lambda;            /* The array of all lambda values                    */
    int      lambda_neighbors;      /* The number of neighboring lambda states to
                                       calculate the energy for in up and down directions
                                       (-1 for all) */
    int      lambda_start_n;        /* The first lambda to calculate energies for */
    int      lambda_stop_n;         /* The last lambda +1 to calculate energies for */
    Real     sc_alpha;              /* free energy soft-core parameter                   */
    int      sc_power;              /* lambda power for soft-core interactions           */
    Real     sc_r_power;            /* r power for soft-core interactions                */
    Real     sc_sigma;              /* free energy soft-core sigma when c6 or c12=0      */
    Real     sc_sigma_min;          /* free energy soft-core sigma for ?????             */
    gmx_bool bScCoul;               /* use softcore for the coulomb portion as well (default FALSE) */
    gmx_bool separate_dvdl[efptNR]; /* whether to print the dvdl term associated with
                                       this term; if it is not specified as separate,
                                       it is lumped with the FEP term */
    int    separate_dhdl_file;      /* whether to write a separate dhdl.xvg file
                                       note: NOT a gmx_bool, but an enum */
    int    dhdl_derivatives;        /* whether to calculate+write dhdl derivatives
                                       note: NOT a gmx_bool, but an enum */
    int    dh_hist_size;            /* The maximum table size for the dH histogram */
    double dh_hist_spacing;         /* The spacing for the dH histogram */
} t_lambda;

typedef struct t_expanded {
    int      nstexpanded;         /* The frequency of expanded ensemble state changes */
    int      elamstats;           /* which type of move updating do we use for lambda monte carlo (or no for none) */
    int      elmcmove;            /* what move set will be we using for state space moves */
    int      elmceq;              /* the method we use to decide of we have equilibrated the weights */
    int      equil_n_at_lam;      /* the minumum number of samples at each lambda for deciding whether we have reached a minimum */
    Real     equil_wl_delta;      /* WL delta at which we stop equilibrating weights */
    Real     equil_ratio;         /* use the ratio of weights (ratio of minimum to maximum) to decide when to stop equilibrating */
    int      equil_steps;         /* after equil_steps steps we stop equilibrating the weights */
    int      equil_samples;       /* after equil_samples total samples (steps/nstfep), we stop equilibrating the weights */
    int      lmc_seed;            /* random number seed for lambda mc switches */
    gmx_bool minvar;              /* whether to use minumum variance weighting */
    int      minvarmin;           /* the number of samples needed before kicking into minvar routine */
    Real     minvar_const;        /* the offset for the variance in MinVar */
    int      c_range;             /* range of cvalues used for BAR */
    gmx_bool bSymmetrizedTMatrix; /* whether to print symmetrized matrices */
    int      nstTij;              /* How frequently to print the transition matrices */
    int      lmc_repeats;         /* number of repetitions in the MC lambda jumps */  /*MRS -- VERIFY THIS */
    int      lmc_forced_nstart;   /* minimum number of samples for each state before free sampling */ /* MRS -- VERIFY THIS! */
    int      gibbsdeltalam;       /* distance in lambda space for the gibbs interval */
    Real     wl_scale;            /* scaling factor for wang-landau */
    Real     wl_ratio;            /* ratio between largest and smallest number for freezing the weights */
    Real     init_wl_delta;       /* starting delta for wang-landau */
    gmx_bool bWLoneovert;         /* use one over t convergece for wang-landau when the delta get sufficiently small */
    gmx_bool bInit_weights;       /* did we initialize the weights? TODO: REMOVE FOR 5.0, no longer needed with new logic */
    Real     mc_temp;             /* To override the main temperature, or define it if it's not defined */
    Real    *init_lambda_weights; /* user-specified initial weights to start with  */
} t_expanded;


/* Abstract types for enforced rotation only defined in pull_rotation.c       */
typedef struct gmx_enfrot *gmx_enfrot_t;
typedef struct gmx_enfrotgrp *gmx_enfrotgrp_t;

typedef struct {
    int         eType;             /* Rotation type for this group                  */
    int         bMassW;            /* Use mass-weighed positions?                   */
    int         nat;               /* Number of atoms in the group                  */
    int        *ind;               /* The global atoms numbers                      */
    rvec       *x_ref;             /* The reference positions                       */
    rvec        vec;               /* The normalized rotation vector                */
    Real        rate;              /* Rate of rotation (degree/ps)                  */
    Real        k;                 /* Force constant (kJ/(mol nm^2)                 */
    rvec        pivot;             /* Pivot point of rotation axis (nm)             */
    int         eFittype;          /* Type of fit to determine actual group angle   */
    int         PotAngle_nstep;    /* Number of angles around the reference angle
                                      for which the rotation potential is also
                                      evaluated (for fit type 'potential' only)     */
    Real            PotAngle_step; /* Distance between two angles in degrees (for
                                      fit type 'potential' only)                    */
    Real            slab_dist;     /* Slab distance (nm)                            */
    Real            min_gaussian;  /* Minimum value the gaussian must have so that
                                      the force is actually evaluated               */
    Real            eps;           /* Additive constant for radial motion2 and
                                      flexible2 potentials (nm^2)                   */
    gmx_enfrotgrp_t enfrotgrp;     /* Stores non-inputrec rotation data per group   */
} t_rotgrp;

typedef struct t_rot {
    int          ngrp;       /* Number of rotation groups                     */
    int          nstrout;    /* Output frequency for main rotation outfile    */
    int          nstsout;    /* Output frequency for per-slab data            */
    t_rotgrp    *grp;        /* Groups to rotate                              */
    gmx_enfrot_t enfrot;     /* Stores non-inputrec enforced rotation data    */
} t_rot;

/* Abstract type for IMD only defined in IMD.c */
struct t_gmx_IMD;

typedef struct t_IMD {
    int               nat;   /* Number of interactive atoms                   */
    int              *ind;   /* The global indices of the interactive atoms   */
    struct t_gmx_IMD *setup; /* Stores non-inputrec IMD data                  */
} t_IMD;

/* Abstract types for position swapping only defined in swapcoords.cpp */
typedef struct t_swap *gmx_swapcoords_t;

typedef struct t_swapGroup {
    char            *molname;             /* Name of the swap group, e.g. NA, CL, SOL       */
    int              nat;                 /* Number of atoms in this group                  */
    int             *ind;                 /* The global ion group atoms numbers             */
    int              nmolReq[eCompNR];    /* Requested number of molecules of this type
                                             per compartment                                */
} t_swapGroup;

typedef struct t_swapcoords {
    int               nstswap;             /* Every how many steps a swap is attempted?    */
    gmx_bool          massw_split[2];      /* Use mass-weighted positions in split group?  */
    Real              cyl0r, cyl1r;        /* Split cylinders defined by radius, upper and */
    Real              cyl0u, cyl1u;        /* ... lower extension. The split cylinders de- */
    Real              cyl0l, cyl1l;        /* ... fine the channels and are each anchored  */
                                           /* ... in the center of the split group         */
    int               nAverage;            /* Coupling constant (nr of swap attempt steps) */
    Real              threshold;           /* Ion counts may deviate from the requested
                                              values by +-threshold before a swap is done  */
    Real              bulkOffset[eCompNR]; /* Offset of the swap layer (='bulk') w.r.t.
                                              the compartment-defining layers              */
    int               ngrp;                /* Number of groups to be controlled            */
    t_swapGroup      *grp;                 /* All swap groups, including split and solvent */
    gmx_swapcoords_t  si_priv;             /* swap private data accessible in
                                            * swapcoords.cpp                               */
} t_swapcoords;

struct t_inputrec
{
    t_inputrec();
    explicit t_inputrec(const t_inputrec &) = delete;
    t_inputrec &operator=(const t_inputrec &) = delete;
    ~t_inputrec();

    int             eI;                      /* Integration method                 */
    gmx_int64_t     nsteps;                  /* number of steps to be taken			*/
    int             simulation_part;         /* Used in checkpointing to separate chunks */
    gmx_int64_t     init_step;               /* start at a stepcount >0 (used w. convert-tpr)    */
    int             nstcalcenergy;           /* frequency of energy calc. and T/P coupl. upd.	*/
    int             cutoff_scheme;           /* group or verlet cutoffs     */
    int             ns_type;                 /* which ns method should we use?               */
    int             nstlist;                 /* number of steps before pairlist is generated	*/
    int             ndelta;                  /* number of cells per rlong			*/
    int             nstcomm;                 /* number of steps after which center of mass	*/
                                             /* motion is removed				*/
    int             comm_mode;               /* Center of mass motion removal algorithm      */
    int             nstlog;                  /* number of steps after which print to logfile	*/
    int             nstxout;                 /* number of steps after which X is output	*/
    int             nstvout;                 /* id. for V					*/
    int             nstfout;                 /* id. for F					*/
    int             nstenergy;               /* number of steps after which energies printed */
    int             nstxout_compressed;      /* id. for compressed trj (.xtc,.tng)           */
    double          init_t;                  /* initial time (ps)              */
    double          delta_t;                 /* time step (ps)				*/
    Real            x_compression_precision; /* precision of x in compressed trajectory file */
    Real            fourier_spacing;         /* requested fourier_spacing, when nk? not set  */
    int             nkx, nky, nkz;           /* number of k vectors in each spatial dimension*/
                                             /* for fourier methods for long range electrost.*/
    int             pme_order;               /* interpolation order for PME                  */
    Real            ewald_rtol;              /* Real space tolerance for Ewald, determines   */
                                             /* the Real/reciprocal space relative weight    */
    Real            ewald_rtol_lj;           /* Real space tolerance for LJ-Ewald            */
    int             ewald_geometry;          /* normal/3d ewald, or pseudo-2d LR corrections */
    Real            epsilon_surface;         /* Epsilon for PME dipole correction            */
    int             ljpme_combination_rule;  /* Type of combination rule in LJ-PME          */
    int             ePBC;                    /* Type of periodic boundary conditions		*/
    int             bPeriodicMols;           /* Periodic molecules                           */
    gmx_bool        bContinuation;           /* Continuation run: starting state is correct	*/
    int             etc;                     /* temperature coupling               */
    int             nsttcouple;              /* interval in steps for temperature coupling   */
    gmx_bool        bPrintNHChains;          /* whether to print nose-hoover chains        */
    int             epc;                     /* pressure coupling                            */
    int             epct;                    /* pressure coupling type			*/
    int             nstpcouple;              /* interval in steps for pressure coupling      */
    Real            tau_p;                   /* pressure coupling time (ps)			*/
    tensor          ref_p;                   /* reference pressure (kJ/(mol nm^3))		*/
    tensor          compress;                /* compressability ((mol nm^3)/kJ)        */
    int             refcoord_scaling;        /* How to scale absolute reference coordinates  */
    rvec            posres_com;              /* The COM of the posres atoms                  */
    rvec            posres_comB;             /* The B-state COM of the posres atoms          */
    int             andersen_seed;           /* Random seed for Andersen thermostat (obsolete) */
    Real            verletbuf_tol;           /* Per atom pair energy drift tolerance (kJ/mol/ps/atom) for list buffer  */
    Real            rlist;                   /* short range pairlist cut-off (nm)		*/
    Real            rtpi;                    /* Radius for test particle insertion           */
    int             coulombtype;             /* Type of electrostatics treatment             */
    int             coulomb_modifier;        /* Modify the Coulomb interaction              */
    Real            rcoulomb_switch;         /* Coulomb switch range start (nm)		*/
    Real            rcoulomb;                /* Coulomb cutoff (nm)		                */
    Real            epsilon_r;               /* relative dielectric constant                 */
    Real            epsilon_rf;              /* relative dielectric constant of the RF       */
    int             implicit_solvent;        /* No (=explicit water), or GBSA solvent models */
    int             gb_algorithm;            /* Algorithm to use for calculation Born radii  */
    int             nstgbradii;              /* Frequency of updating Generalized Born radii */
    Real            rgbradii;                /* Cutoff for GB radii calculation              */
    Real            gb_saltconc;             /* Salt concentration (M) for GBSA models       */
    Real            gb_epsilon_solvent;      /* dielectric coeff. of implicit solvent     */
    Real            gb_obc_alpha;            /* 1st scaling factor for Bashford-Case GB      */
    Real            gb_obc_beta;             /* 2nd scaling factor for Bashford-Case GB      */
    Real            gb_obc_gamma;            /* 3rd scaling factor for Bashford-Case GB      */
    Real            gb_dielectric_offset;    /* Dielectric offset for Still/HCT/OBC     */
    int             sa_algorithm;            /* Algorithm for SA part of GBSA                */
    Real            sa_surface_tension;      /* Energy factor for SA part of GBSA */
    int             vdwtype;                 /* Type of Van der Waals treatment              */
    int             vdw_modifier;            /* Modify the VdW interaction                   */
    Real            rvdw_switch;             /* Van der Waals switch range start (nm)        */
    Real            rvdw;                    /* Van der Waals cutoff (nm)	        */
    int             eDispCorr;               /* Perform Long range dispersion corrections    */
    Real            tabext;                  /* Extension of the table beyond the cut-off,   *
                                              * as well as the table length for 1-4 interac. */
    Real            shake_tol;               /* tolerance for shake				*/
    int             efep;                    /* free energy calculations                     */
    t_lambda       *fepvals;                 /* Data for the FEP state                       */
    gmx_bool        bSimTemp;                /* Whether to do simulated tempering            */
    t_simtemp      *simtempvals;             /* Variables for simulated tempering            */
    gmx_bool        bExpanded;               /* Whether expanded ensembles are used          */
    t_expanded     *expandedvals;            /* Expanded ensemble parameters              */
    int             eDisre;                  /* Type of distance restraining                 */
    Real            dr_fc;                   /* force constant for ta_disre			*/
    int             eDisreWeighting;         /* type of weighting of pairs in one restraints	*/
    gmx_bool        bDisreMixed;             /* Use comb of time averaged and instan. viol's	*/
    int             nstdisreout;             /* frequency of writing pair distances to enx   */
    Real            dr_tau;                  /* time constant for memory function in disres    */
    Real            orires_fc;               /* force constant for orientational restraints  */
    Real            orires_tau;              /* time constant for memory function in orires    */
    int             nstorireout;             /* frequency of writing tr(SD) to enx           */
    Real            em_stepsize;             /* The stepsize for updating			*/
    Real            em_tol;                  /* The tolerance				*/
    int             niter;                   /* Number of iterations for convergence of      */
                                             /* steepest descent in relax_shells             */
    Real            fc_stepsize;             /* Stepsize for directional minimization        */
                                             /* in relax_shells                              */
    int             nstcgsteep;              /* number of steps after which a steepest       */
                                             /* descents step is done while doing cg         */
    int             nbfgscorr;               /* Number of corrections to the hessian to keep */
    int             eConstrAlg;              /* Type of constraint algorithm                 */
    int             nProjOrder;              /* Order of the LINCS Projection Algorithm      */
    Real            LincsWarnAngle;          /* If bond rotates more than %g degrees, warn   */
    int             nLincsIter;              /* Number of iterations in the final Lincs step */
    gmx_bool        bShakeSOR;               /* Use successive overrelaxation for shake      */
    Real            bd_fric;                 /* Friction coefficient for BD (amu/ps)         */
    gmx_int64_t     ld_seed;                 /* Random seed for SD and BD                    */
    int             nwall;                   /* The number of walls                          */
    int             wall_type;               /* The type of walls                            */
    Real            wall_r_linpot;           /* The potentail is linear for r<=wall_r_linpot */
    int             wall_atomtype[2];        /* The atom type for walls                      */
    Real            wall_density[2];         /* Number density for walls                     */
    Real            wall_ewald_zfac;         /* Scaling factor for the box for Ewald         */

    /* COM pulling data */
    gmx_bool              bPull;             /* Do we do COM pulling?                        */
    struct pull_params_t *pull;              /* The data for center of mass pulling          */
    // TODO: Remove this by converting pull into a ForceProvider
    struct pull_t        *pull_work;         /* The COM pull force calculation data structure */

    /* AWH bias data */
    gmx_bool                 bDoAwh;    /* Use awh biasing for PMF calculations?        */
    gmx::AwhParams          *awhParams; /* AWH biasing parameters                       */
    // TODO: Remove this by converting AWH into a ForceProvider
    gmx::Awh                *awh;       /* AWH work object */

    /* Enforced rotation data */
    gmx_bool                 bRot;           /* Calculate enforced rotation potential(s)?    */
    t_rot                   *rot;            /* The data for enforced rotation potentials    */

    int                      eSwapCoords;    /* Do ion/water position exchanges (CompEL)?    */
    t_swapcoords            *swap;

    gmx_bool                 bIMD;           /* Allow interactive MD sessions for this .tpr? */
    t_IMD                   *imd;            /* Interactive molecular dynamics               */

    Real                     cos_accel;      /* Acceleration for viscosity calculation       */
    tensor                   deform;         /* Triclinic deformation velocities (nm/ps)     */
    int                      userint1;       /* User determined parameters                   */
    int                      userint2;
    int                      userint3;
    int                      userint4;
    Real                     userreal1;
    Real                     userreal2;
    Real                     userreal3;
    Real                     userreal4;
    t_grpopts                opts;          /* Group options				*/
    gmx_bool                 bQMMM;         /* QM/MM calculation                            */
    int                      QMconstraints; /* constraints on QM bonds                      */
    int                      QMMMscheme;    /* Scheme: ONIOM or normal                      */
    Real                     scalefactor;   /* factor for scaling the MM charges in QM calc.*/

    /* Fields for removed features go here (better caching) */
    gmx_bool                 bAdress;      // Whether AdResS is enabled - always false if a valid .tpr was read
    gmx_bool                 useTwinRange; // Whether twin-range scheme is active - always false if a valid .tpr was read

    gmx::KeyValueTreeObject *params;
};

int ir_optimal_nstcalcenergy(const t_inputrec *ir);

int tcouple_min_integration_steps(int etc);

int ir_optimal_nsttcouple(const t_inputrec *ir);

int pcouple_min_integration_steps(int epc);

int ir_optimal_nstpcouple(const t_inputrec *ir);

/* Returns if the Coulomb force or potential is switched to zero */
gmx_bool ir_coulomb_switched(const t_inputrec *ir);

/* Returns if the Coulomb interactions are zero beyond the rcoulomb.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_coulomb_is_zero_at_cutoff(const t_inputrec *ir);

/* As ir_coulomb_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rcoulomb.
 */
gmx_bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec *ir);

/* Returns if the Van der Waals force or potential is switched to zero */
gmx_bool ir_vdw_switched(const t_inputrec *ir);

/* Returns if the Van der Waals interactions are zero beyond the rvdw.
 * Note: always returns TRUE for the Verlet cut-off scheme.
 */
gmx_bool ir_vdw_is_zero_at_cutoff(const t_inputrec *ir);

/* As ir_vdw_is_zero_at_cutoff, but also returns TRUE for user tabulated
 * interactions, since these might be zero beyond rvdw.
 */
gmx_bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec *ir);

/*! \brief Free memory from input record.
 *
 * All arrays and pointers will be freed.
 *
 * \param[in] ir The data structure
 */
void done_inputrec(t_inputrec *ir);

void pr_inputrec(FILE *fp, int indent, const char *title, const t_inputrec *ir,
                 gmx_bool bMDPformat);

void cmp_inputrec(FILE *fp, const t_inputrec *ir1, const t_inputrec *ir2, Real ftol, Real abstol);

void comp_pull_AB(FILE *fp, pull_params_t *pull, Real ftol, Real abstol);


gmx_bool inputrecDeform(const t_inputrec *ir);

gmx_bool inputrecDynamicBox(const t_inputrec *ir);

gmx_bool inputrecPreserveShape(const t_inputrec *ir);

gmx_bool inputrecNeedMutot(const t_inputrec *ir);

gmx_bool inputrecTwinRange(const t_inputrec *ir);

gmx_bool inputrecExclForces(const t_inputrec *ir);

gmx_bool inputrecNptTrotter(const t_inputrec *ir);

gmx_bool inputrecNvtTrotter(const t_inputrec *ir);

gmx_bool inputrecNphTrotter(const t_inputrec *ir);

/*! \brief Return true if the simulation is 2D periodic with two walls. */
bool     inputrecPbcXY2Walls(const t_inputrec *ir);

/* Returns true for MD integator with T and/or P-coupling that supports
 * calculating the conserved energy quantity.
 */
bool integratorHasConservedEnergyQuantity(const t_inputrec *ir);

/*! \brief Return the number of bounded dimensions
 *
 * \param[in] ir The input record with MD parameters
 * \return the number of dimensions in which
 * the coordinates of the particles are bounded, starting at X.
 */
int inputrec2nboundeddim(const t_inputrec *ir);

/*! \brief Returns the number of degrees of freedom in center of mass motion
 *
 * \param[in] ir the inputrec structure
 * \return the number of degrees of freedom of the center of mass
 */
int ndof_com(const t_inputrec *ir);

#endif /* GMX_MDTYPES_INPUTREC_H */
