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

//  $Id: petsc_utils.cc,v 1.5 2008/07/09 05:58:16 gdiso Exp $

#include <map>
#include <vector>

#include "genius_petsc.h"
#include "petsc_utils.h"


namespace PetscUtils
{
  /**
   * wrap to petsc MatZeroRows
   */
  PetscErrorCode MatZeroRows(Mat mat,PetscInt numRows,const PetscInt rows[],PetscScalar diag)
  {
#if PETSC_VERSION_GE(3,2,0)
      return ::MatZeroRows(mat, numRows, rows, diag, PETSC_NULL, PETSC_NULL);
#else
      return ::MatZeroRows(mat, numRows, rows, diag);
#endif
  }

  /*-------------------------------------------------------------------
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
  PetscErrorCode  VecAddRowToRow(Vec vec, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha[])
  {
    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    if( rows>0 )
    {
      std::vector<PetscScalar> y(rows);
      VecGetValues(vec, rows, src_rows, &y[0]);

      for(PetscInt nrow=0; nrow<rows; ++nrow)
        y[nrow] *= alpha[nrow];

      // add value to f
      VecSetValues(vec, rows, dst_rows, &y[0], ADD_VALUES);
    }

    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    return 0;
  }


  /*-------------------------------------------------------------------
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  vec        Petsc Vector
   * @param  rows       number of operation rows
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier.
   *
   * @note   If vec is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  PetscErrorCode  VecAddRowToRow(Vec vec, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha)
  {
    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    if( rows>0 )
    {
      std::vector<PetscScalar> y(rows);
      VecGetValues(vec, rows, src_rows, &y[0]);

      for(PetscInt nrow=0; nrow<rows; ++nrow)
        y[nrow] *= alpha;

      // add value to f
      VecSetValues(vec, rows, dst_rows, &y[0], ADD_VALUES);
    }

    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    return 0;
  }


  /*-------------------------------------------------------------------
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  vec        Petsc Vector
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier.
   *
   * @note   If vec is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  PetscErrorCode  VecAddRowToRow(Vec vec, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, PetscScalar alpha)
  {
    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    if( src_rows.size() )
    {
      assert(src_rows.size() == dst_rows.size());

      std::vector<PetscScalar> y(src_rows.size());
      VecGetValues(vec, src_rows.size(), &src_rows[0], &y[0]);

      for(unsigned int nrow=0; nrow<src_rows.size(); ++nrow)
        y[nrow] *= alpha;

      // add value to f
      VecSetValues(vec, dst_rows.size(), &dst_rows[0], &y[0], ADD_VALUES);
    }

    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    return 0;
  }


  /*-------------------------------------------------------------------
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
  PetscErrorCode  VecAddRowToRow(Vec vec, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, std::vector<PetscScalar> & alpha)
  {
    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    if( src_rows.size() )
    {
      assert(src_rows.size() == dst_rows.size() && src_rows.size() == alpha.size() );

      std::vector<PetscScalar> y(src_rows.size());
      VecGetValues(vec, src_rows.size(), &src_rows[0], &y[0]);

      for(unsigned int nrow=0; nrow<src_rows.size(); ++nrow)
        y[nrow] *= alpha[nrow];

      // add value to f
      VecSetValues(vec, dst_rows.size(), &dst_rows[0], &y[0], ADD_VALUES);
    }

    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    return 0;
  }


 /*-------------------------------------------------------------------
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
  PetscErrorCode  VecAddClearRow(Vec vec, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, std::vector<PetscInt> & clear_rows)
  {
    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    assert(src_rows.size() == dst_rows.size());
    if( src_rows.size() )
    {
      // get source value from vec
      std::vector<PetscScalar> y(src_rows.size());
      VecGetValues(vec, src_rows.size(), &src_rows[0], &y[0]);
      // add value to vec destination
      VecSetValues(vec, dst_rows.size(), &dst_rows[0], &y[0], ADD_VALUES);
    }

    // we assembly vec for INSERT_VALUES!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    if( clear_rows.size() )
    {
      std::vector<PetscScalar> y(clear_rows.size(), 0.0);
      VecSetValues(vec, clear_rows.size(), &clear_rows[0], &y[0], INSERT_VALUES);
    }

    // we assembly vec for any case!
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    return 0;
  }

  /*-------------------------------------------------------------------
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
  PetscErrorCode  MatAddRowToRow(Mat mat, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha[])
  {

    // test if the matrix is assembled
    // note: the test is not work properly! if it is a bug...

    //PetscBool assembled;
    //MatAssembled(mat, &assembled);
    //if( !assembled )
    {
      MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    }

    std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > > row_add_to_buffer;

    // read the row
    for(PetscInt nrow=0; nrow<rows; nrow++)
    {
      PetscInt ncols;
      const PetscInt * row_cols_pointer;
      const PetscScalar * row_vals_pointer;

      MatGetRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);

      // save the values
      std::vector<PetscInt>    row_cols;
      std::vector<PetscScalar> row_vals;
      for(PetscInt i=0;i <ncols; i++)
      {
        row_cols.push_back(row_cols_pointer[i]);
        row_vals.push_back(row_vals_pointer[i]*alpha[nrow]);
      }

      row_add_to_buffer.insert(std::pair< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >
                               (dst_rows[nrow], std::pair<std::vector<PetscInt>, std::vector<PetscScalar> >(row_cols, row_vals)));

      // restore pointers
      MatRestoreRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);
    }

    //ok, we add rows to destination
    std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >::iterator row_add_it =  row_add_to_buffer.begin();
    for(; row_add_it !=  row_add_to_buffer.end(); ++row_add_it)
    {
      MatSetValues(mat, 1, &((*row_add_it).first) ,
                   (*row_add_it).second.first.size(), &((*row_add_it).second.first[0]) ,
                   &((*row_add_it).second.second[0]), ADD_VALUES);

    }

    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

    return 0;
  }


  /*-------------------------------------------------------------------
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  mat        Petsc Matrix
   * @param  rows       number of operation rows
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier, default is 1.0
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  PetscErrorCode  MatAddRowToRow(Mat mat, PetscInt rows, PetscInt src_rows[], PetscInt dst_rows[], PetscScalar alpha)
  {

    // test if the matrix is assembled
    // note: the test is not work properly! if it is a bug...

    //PetscBool assembled;
    //MatAssembled(mat, &assembled);

    //if( !assembled )
    {
      MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    }


    std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > > row_add_to_buffer;

    // read the row
    for(PetscInt nrow=0; nrow<rows; nrow++)
    {
      PetscInt ncols;
      const PetscInt * row_cols_pointer;
      const PetscScalar * row_vals_pointer;

      MatGetRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);

      // save the values
      std::vector<PetscInt>    row_cols;
      std::vector<PetscScalar> row_vals;
      for(PetscInt i=0;i <ncols; i++)
      {
        row_cols.push_back(row_cols_pointer[i]);
        row_vals.push_back(row_vals_pointer[i]*alpha);
      }

      row_add_to_buffer.insert(std::pair< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >
                               (dst_rows[nrow], std::pair<std::vector<PetscInt>, std::vector<PetscScalar> >(row_cols, row_vals)));

      // restore pointers
      MatRestoreRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);
    }

    //ok, we add rows to destination
    std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >::iterator row_add_it =  row_add_to_buffer.begin();
    for(; row_add_it !=  row_add_to_buffer.end(); ++row_add_it)
    {
      MatSetValues(mat, 1, &((*row_add_it).first) ,
                   (*row_add_it).second.first.size(), &((*row_add_it).second.first[0]) ,
                   &((*row_add_it).second.second[0]), ADD_VALUES);

    }

    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

    return 0;
  }



  /*-------------------------------------------------------------------
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  mat        Petsc Matrix
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier, default is 1.0
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  PetscErrorCode  MatAddRowToRow(Mat mat, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, PetscScalar alpha)
  {

    // test if the matrix is assembled
    // note: the test is not work properly! if it is a bug...

    //PetscBool assembled;
    //MatAssembled(mat, &assembled);

    //if( !assembled )
    {
      MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    }

    if( src_rows.size() )
    {
      assert(src_rows.size() == dst_rows.size());

      std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > > row_add_to_buffer;

      // read the row
      for(unsigned int nrow=0; nrow<src_rows.size(); nrow++)
      {
        PetscInt ncols;
        const PetscInt * row_cols_pointer;
        const PetscScalar * row_vals_pointer;

        MatGetRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);

        // save the values
        std::vector<PetscInt>    row_cols;
        std::vector<PetscScalar> row_vals;
        for(PetscInt i=0;i <ncols; i++)
        {
          row_cols.push_back(row_cols_pointer[i]);
          row_vals.push_back(row_vals_pointer[i]*alpha);
        }

        row_add_to_buffer.insert(std::pair< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >
                                 (dst_rows[nrow], std::pair<std::vector<PetscInt>, std::vector<PetscScalar> >(row_cols, row_vals)));

        // restore pointers
        MatRestoreRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);
      }

      //ok, we add rows to destination
      std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >::iterator row_add_it =  row_add_to_buffer.begin();
      for(; row_add_it !=  row_add_to_buffer.end(); ++row_add_it)
      {
        MatSetValues(mat, 1, &((*row_add_it).first) ,
                     (*row_add_it).second.first.size(), &((*row_add_it).second.first[0]) ,
                     &((*row_add_it).second.second[0]), ADD_VALUES);

      }

    }

    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

    return 0;
  }


  /*-------------------------------------------------------------------
   * @brief add source rows to destination rows, multi by alpha
   *
   * @param  mat        Petsc Matrix
   * @param  src_rows   source rows
   * @param  dst_rows   the destination rows will be added to
   * @param  alpha      the scalar multiplier, default is 1.0
   *
   * @note   If mat is not assembled, this function will do it for you. the src row should on local processor.
   *
   */
  PetscErrorCode  MatAddRowToRow(Mat mat, std::vector<PetscInt> & src_rows, std::vector<PetscInt> & dst_rows, std::vector<PetscScalar> & alpha)
  {

    // test if the matrix is assembled
    // note: the test is not work properly! if it is a bug...

    //PetscBool assembled;
    //MatAssembled(mat, &assembled);

    //if( !assembled )
    {
      MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    }

    if( src_rows.size() )
    {
      assert(src_rows.size() == dst_rows.size() && src_rows.size() == alpha.size());

      std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > > row_add_to_buffer;

      // read the row
      for(unsigned int nrow=0; nrow<src_rows.size(); nrow++)
      {
        PetscInt ncols;
        const PetscInt * row_cols_pointer;
        const PetscScalar * row_vals_pointer;

        MatGetRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);

        // save the values
        std::vector<PetscInt>    row_cols;
        std::vector<PetscScalar> row_vals;
        for(PetscInt i=0;i <ncols; i++)
        {
          row_cols.push_back(row_cols_pointer[i]);
          row_vals.push_back(row_vals_pointer[i]*alpha[nrow]);
        }

        row_add_to_buffer.insert(std::pair< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >
                                 (dst_rows[nrow], std::pair<std::vector<PetscInt>, std::vector<PetscScalar> >(row_cols, row_vals)));

        // restore pointers
        MatRestoreRow(mat, src_rows[nrow], &ncols, &row_cols_pointer, &row_vals_pointer);
      }

      //ok, we add rows to destination
      std::multimap< PetscInt, std::pair<std::vector<PetscInt>, std::vector<PetscScalar> > >::iterator row_add_it =  row_add_to_buffer.begin();
      for(; row_add_it !=  row_add_to_buffer.end(); ++row_add_it)
      {
        MatSetValues(mat, 1, &((*row_add_it).first) ,
                     (*row_add_it).second.first.size(), &((*row_add_it).second.first[0]) ,
                     &((*row_add_it).second.second[0]), ADD_VALUES);

      }

    }

    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

    return 0;
  }



  /*-------------------------------------------------------------------
   * @brief add real DenseVector to PetscVec by dof_indices
   *
   * @param  vec         Petsc Vector
   * @param  real_vec    small vector
   * @param  dof_indices location for add
   *
   */
  PetscErrorCode  VecAdd(Vec vec, const DenseVector<PetscScalar> &real_vec, const std::vector<PetscInt> & dof_indices)
  {
    assert(real_vec.size()==dof_indices.size());
    VecSetValues(vec, dof_indices.size(), &dof_indices[0], &((real_vec.get_values())[0]), ADD_VALUES);
    return 0;
  }



  /*-------------------------------------------------------------------
   * @brief add complex DenseVector to PetscVec by dof_indices
   *
   * @param  vec         Petsc Vector
   * @param  complex_vec small vector
   * @param  dof_indices location for add
   *
   */
//   PetscErrorCode  VecAdd(Vec vec, const DenseVector<Complex> &complex_vec, const std::vector<PetscInt> & dof_indices)
//   {
//     assert(complex_vec.size()*2==dof_indices.size());
//     //convert complex indices to real indices?
//     std::vector<Real> real_vec(complex_vec.size()*2);
//     for(unsigned int n=0; n<complex_vec.size(); ++n)
//     {
//       real_vec[2*n]   = complex_vec(n).real();
//       real_vec[2*n+1] = complex_vec(n).imag();
//     }
//     VecSetValues(vec, dof_indices.size(), &dof_indices[0], &(real_vec[0]), ADD_VALUES);
//     return 0;
//   }



  /*-------------------------------------------------------------------
   * @brief add real DenseMatrix to PetscMat by dof_indices
   *
   * @param  mat         Petsc Mat
   * @param  real_mat    small matrix
   * @param  dof_indices location for add
   *
   */
  PetscErrorCode  MatAdd(Mat mat, const DenseMatrix<PetscScalar> &real_mat, const std::vector<PetscInt> & dof_indices)
  {
    assert(real_mat.m()==dof_indices.size());
    assert(real_mat.n()==dof_indices.size());
    MatSetValues(mat, dof_indices.size(), &dof_indices[0], dof_indices.size(), &dof_indices[0], &((real_mat.get_values())[0]), ADD_VALUES);
    return 0;
  }



  /*-------------------------------------------------------------------
   * @brief add complex DenseMatrix to PetscMat by dof_indices
   *
   * @param  mat         Petsc Mat
   * @param  complex_mat small matrix
   * @param  dof_indices location for add
   *
   */
//   PetscErrorCode  MatAdd(Mat mat, const DenseMatrix<Complex> &complex_mat, const std::vector<PetscInt> & dof_indices)
//   {
//     assert(complex_mat.m()*2==dof_indices.size());
//     assert(complex_mat.n()*2==dof_indices.size());
//
//     //convert complex indices to real indices?
//     DenseMatrix<Real> real_mat(complex_mat.m()*2, complex_mat.n()*2);
//     for(unsigned int m=0; m<complex_mat.m(); ++m)
//     {
//       for(unsigned int n=0; n<complex_mat.n(); ++n)
//       {
//         real_mat.el(2*m,  2*n)   =   complex_mat.el(m,n).real();
//         real_mat.el(2*m,  2*n+1) = - complex_mat.el(m,n).imag();
//         real_mat.el(2*m+1,2*n)   =   complex_mat.el(m,n).imag();
//         real_mat.el(2*m+1,2*n+1) =   complex_mat.el(m,n).real();
//       }
//     }
//
//     MatSetValues(mat, dof_indices.size(), &dof_indices[0], dof_indices.size(), &dof_indices[0], &((real_mat.get_values())[0]), ADD_VALUES);
//     return 0;
//   }

}

