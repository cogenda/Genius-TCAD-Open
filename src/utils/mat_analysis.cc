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

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

#include "petscmat.h"
#include "parallel.h"




void  mat_analysis(const Mat mat)
{
  PetscBool assembled;
  MatAssembled(mat, &assembled);
  if( !assembled )
  {
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  }

  PetscInt m;     //  - the number of local rows
  PetscInt n;     //  - the number of local columns
  MatGetLocalSize(mat, &m, &n);
  PetscInt M, N;      // - the number of global rows and columns
  MatGetSize(mat, &M, &N);

  PetscInt row_begin ;  // - the global index of the first local row
  PetscInt row_end;   // - one more than the global index of the last local row
  MatGetOwnershipRange(mat, &row_begin,  &row_end);

  // matrix entry (row, col) pairs
  std::vector<PetscInt> rows;
  std::vector<PetscInt> cols;
  std::vector<PetscScalar> values;

  // read the row
  for(PetscInt row=row_begin; row<row_end; row++)
  {
    PetscInt ncols;
    const PetscInt * row_cols_pointer;
    const PetscScalar * row_vals_pointer;

    MatGetRow(mat, row, &ncols, &row_cols_pointer, &row_vals_pointer);

    for(PetscInt c=0; c<ncols; c++)
    {
      PetscInt col = row_cols_pointer[c];
      rows.push_back(row);
      cols.push_back(col);
      values.push_back(fabs(row_vals_pointer[c]));
    }

    // restore pointers
    MatRestoreRow(mat, row, &ncols, &row_cols_pointer, &row_vals_pointer);
  }

  // gather from other processor
  Parallel::gather(0, rows);
  Parallel::gather(0, cols);
  Parallel::gather(0, values);

  if( Genius::is_first_processor() )
  {
    std::cout<< "Matrix has dimension of " << M << "x" << N <<" with " << values.size() << " fill in values" << std::endl;
    // find the max abs value of matrix
    std::cout<< "Max matrix entry is: " << *( std::max_element( values.begin(), values.end() ) ) << std::endl;

    // find the unsymmetric rate of the matrix
    std::set< std::pair<PetscInt, PetscInt> > mat_entry_set;
    for(unsigned int n=0; n<values.size(); ++n)
      mat_entry_set.insert(std::make_pair(rows[n], cols[n]));

    unsigned int n_unsymmetric_entry = 0;
    for(unsigned int n=0; n<values.size(); ++n)
    {
      if( mat_entry_set.find(std::make_pair(cols[n], rows[n])) == mat_entry_set.end() )
        n_unsymmetric_entry++;
    }
    std::cout<< "Unsymmetric matrix entries: " << n_unsymmetric_entry
    << " , rate " << n_unsymmetric_entry/double(values.size()) << std::endl;
  }
}


#ifdef HAVE_TIFF

#include <tiffio.h>

void mat_to_image(const Mat mat, const std::string &image_file)
{
  PetscBool assembled;
  MatAssembled(mat, &assembled);
  if( !assembled )
  {
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
  }

  PetscInt m;     //  - the number of local rows
  PetscInt n;     //  - the number of local columns
  MatGetLocalSize(mat, &m, &n);
  PetscInt M, N;      // - the number of global rows and columns
  MatGetSize(mat, &M, &N);

  PetscInt row_begin ;  // - the global index of the first local row
  PetscInt row_end;   // - one more than the global index of the last local row
  MatGetOwnershipRange(mat, &row_begin,  &row_end);

  // matrix entry (row, col) pairs
  std::vector<PetscInt> rows;
  std::vector<PetscInt> cols;
  std::vector<PetscScalar> values;

  // read the row
  for(PetscInt row=row_begin; row<row_end; row++)
  {
    PetscInt ncols;
    const PetscInt * row_cols_pointer;
    const PetscScalar * row_vals_pointer;

    MatGetRow(mat, row, &ncols, &row_cols_pointer, &row_vals_pointer);

    for(PetscInt c=0; c<ncols; c++)
    {
      PetscInt col = row_cols_pointer[c];
      rows.push_back(row);
      cols.push_back(col);
      values.push_back(fabs(row_vals_pointer[c]));
    }

    // restore pointers
    MatRestoreRow(mat, row, &ncols, &row_cols_pointer, &row_vals_pointer);
  }

  // gather from other processor
  Parallel::gather(0, rows);
  Parallel::gather(0, cols);
  Parallel::gather(0, values);

  const PetscScalar max_value = *( std::max_element( values.begin(), values.end() ) );
  std::vector< std::vector< std::pair<PetscInt, PetscScalar> > > entry(N);
  for(unsigned int n=0; n<rows.size(); ++n)
  {
    PetscInt row = rows[n];
    PetscInt col = cols[n];
    entry[row].push_back( std::make_pair(col, values[n]) );
  }
  std::vector< std::vector< std::pair<PetscInt, PetscScalar> > > pad_entry(N);
  {
    std::set< std::pair<PetscInt, PetscInt> > mat_entry_set;
    for(unsigned int n=0; n<values.size(); ++n)
      mat_entry_set.insert(std::make_pair(rows[n], cols[n]));

    for(unsigned int n=0; n<values.size(); ++n)
    {
      if( mat_entry_set.find(std::make_pair(cols[n], rows[n])) == mat_entry_set.end() )
        pad_entry[cols[n]].push_back( std::make_pair(rows[n], 1.0) );
    }
  }

  rows.clear();
  cols.clear();
  values.clear();


  TIFF *fout= TIFFOpen(image_file.c_str(), "w");
  int sampleperpixel = 3;    //RGB channel

  TIFFSetField(fout, TIFFTAG_IMAGEWIDTH, M);  // set the width of the image
  TIFFSetField(fout, TIFFTAG_IMAGELENGTH, N);    // set the height of the image
  TIFFSetField(fout, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
  TIFFSetField(fout, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
  TIFFSetField(fout, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
  TIFFSetField(fout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(fout, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(fout, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);

  tsize_t linebytes = sampleperpixel*M;     // length in memory of one row of pixel in the image.
  // We set the strip size of the file to be size of one row of pixels
  TIFFSetField(fout, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(fout, linebytes));

  unsigned int rainbow[20][3] =
    {
      0, 0, 65535,
      0, 29298, 65535,
      0, 41377, 65535,
      0, 50629, 65535,
      0, 58596, 65535,
      0, 65535, 65535,
      0, 65535, 56797,
      0, 65535, 46003,
      0, 65535, 32639,
      0, 65535, 0,
      32125, 65535, 0,
      46260, 65535, 0,
      56540, 65535, 0,
      65535, 65535, 0,
      65535, 59881, 0,
      65535, 53199, 0,
      65535, 45746, 0,
      65535, 38036, 0,
      65535, 26471, 0,
      65535, 0, 0
    };

  //Now writing image to the file one strip at a time
  for (PetscInt row = 0; row < N; row++)
  {
    std::vector<unsigned char> image(linebytes, 255);
    const std::vector< std::pair<PetscInt, PetscScalar> > & one_row = entry[row];
    for(unsigned int n=0; n<one_row.size(); ++n)
    {
      int col = one_row[n].first;
      PetscScalar value = one_row[n].second;
      int level = 19*(value/max_value);//map to rainbow
      image[3*col + 0] = rainbow[level][0]/255;
      image[3*col + 1] = rainbow[level][1]/255;
      image[3*col + 2] = rainbow[level][2]/255;
    }
    TIFFWriteScanline(fout, &image[0], row, 0);
  }
  TIFFClose(fout);
}

#else
void mat_to_image(const Mat , const std::string &)
{}
#endif
