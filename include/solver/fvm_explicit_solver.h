/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/



#ifndef __fvm_explicit_solver_h__
#define __fvm_explicit_solver_h__

#include "enum_petsc_type.h"
#include "fvm_pde_solver.h"
#include "petscvec.h"


/**
 * The explicit solver contex.
 */
class FVM_ExplicitSolver : public FVM_PDESolver
{
  public:

    /**
     * construct the explicit solver contex:
     * using petsc vec for storage
     * here we build vec as well as parallel scatter
     */
    FVM_ExplicitSolver ( SimulationSystem & system );

    /**
     * free all the contex
     */
    virtual ~FVM_ExplicitSolver();

    /**
     * setup explicit context
     */
    void setup_explicit_data();

    /**
     * clear explicit data
     */
    void clear_explicit_data();

  protected:

    /**
     * the global solution vector
     */
    Vec            x;

    /**
     * the previous global solution vector
     */
    Vec            x_prev;

    /**
     * the local solution vector
     */
    Vec            lx;

    /**
     * the global volume vector
     */
    Vec            vol;

    /**
     * the global flux vector
     */
    Vec            f;


    /**
     * the global max time step vector
     */
    Vec            t;


    /**
     * the local time step vector
     */
    Vec            lt;


    /**
     * the global working (temp) vector
     */
    Vec            w;

    /**
     * the global index set of solution vector this processor needs
     */
    IS             gis;

    /**
     * the local index set of solution vector this processor needs
     */
    IS             lis;

    /**
     * data structure used to map global solution vector to local vector
     */
    VecScatter     scatter;
};

#endif

