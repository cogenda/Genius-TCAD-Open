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

//  $Id: petsc_utils.h,v 1.3 2008/07/09 05:58:16 gdiso Exp $

#ifndef __petsc_utils_h__
#define __petsc_utils_h__



#include "petscmat.h"

#include "genius_common.h"
#include "dense_vector.h"
#include "dense_matrix.h"

namespace PetscUtils
{

  /**
   * wrap to petsc MatZeroRows
   */
  extern PetscErrorCode  MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag);


  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  vec        Petsc Vector
   * @param  rows       number of operation rows
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier array.
   *
   * @note   If vec is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  VecAddRowToRow(Vec vec, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha[]);

  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  vec        Petsc Vector
   * @param  rows       number of operation rows
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier, default is 1.0.
   *
   * @note   If vec is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  VecAddRowToRow(Vec vec, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha=1.0);


  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  vec        Petsc Vector
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier, default is 1.0.
   *
   * @note   If vec is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  VecAddRowToRow(Vec vec, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, PetscScalar alpha=1.0);


  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  vec        Petsc Vector
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier array.
   *
   * @note   If vec is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  VecAddRowToRow(Vec vec, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, std::vector<PetscScalar> & alpha);


 /**-------------------------------------------------------------------
  * @brief add source rows to destination rows, and clear some rows
  *
  * @param  vec        Petsc Vector
  * @param  src_rows   source rows
  * @param  dst_rows   the destination rows will be added to
  * @param  clear_rows the rows to be cleared .
  *
  * @note   If vec is not assembled, this function will do it for you. the src row should on local processor.
  *
  */
  extern PetscErrorCode  VecAddClearRow(Vec vec, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, std::vector<PetscInt> & clear_rows);

  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  mat        Petsc Matrix
   * @param  rows       number of operation rows
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier array.
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  MatAddRowToRow(Mat mat, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha[]);

  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  mat        Petsc Matrix
   * @param  rows       number of operation rows
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier, default is 1.0.
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  MatAddRowToRow(Mat mat, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha=1.0);

  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  mat        Petsc Matrix
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier, default is 1.0.
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  MatAddRowToRow(Mat mat, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, PetscScalar alpha=1.0);

  /**
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  mat        Petsc Matrix
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier array
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  extern PetscErrorCode  MatAddRowToRow(Mat mat, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, std::vector<PetscScalar> & alpha);


  /**
   * @brief add source rows to destination rows, and clear some rows
   *
   * @param  mat        Petsc Matrix
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  clear_rows the rows to be cleared
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  //extern PetscErrorCode  MatAddClearRow(Mat mat, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, std::vector<PetscInt> & clear_rows);

  /**
   * @brief add real DenseVector to PetscVec by dof_indices
   *
   * @param  vec         Petsc Vector
   * @param  real_vec    small vector
   * @param  dof_indices location for add
   *
   */
  extern PetscErrorCode  VecAdd(Vec vec, const DenseVector<PetscScalar> &real_vec, const std::vector<PetscInt> & dof_indices);

  /**
   * @brief add complex DenseVector to PetscVec by dof_indices
   *
   * @param  vec         Petsc Vector
   * @param  complex_vec small vector
   * @param  dof_indices location for add
   *
   */
  //extern PetscErrorCode  VecAdd(Vec vec, const DenseVector<Complex> &complex_vec, const std::vector<PetscInt> & dof_indices);

  /**
   * @brief add real DenseMatrix to PetscMat by dof_indices
   *
   * @param  mat         Petsc Mat
   * @param  real_mat    small matrix
   * @param  dof_indices location for add
   *
   */
  extern PetscErrorCode  MatAdd(Mat mat, const DenseMatrix<PetscScalar> &real_mat, const std::vector<PetscInt> & dof_indices);

  /**
   * @brief add complex DenseMatrix to PetscMat by dof_indices
   *
   * @param  mat         Petsc Mat
   * @param  complex_mat small matrix
   * @param  dof_indices location for add
   *
   */
  //extern PetscErrorCode  MatAdd(Mat mat, const DenseMatrix<Complex> &complex_mat, const std::vector<PetscInt> & dof_indices);

}

#endif //#define __petsc_utils_h__
