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



#ifndef __mesh_generation_quad4_h__
#define __mesh_generation_quad4_h__

#include "mesh_generation_struct.h"

/**
 * the simple quad4 mesh generator
 */
class MeshGeneratorQuad4 : public MeshGeneratorStruct
{
public:
  /**
   * constructor
   */
  MeshGeneratorQuad4(MeshBase& mesh, Parser::InputParser & decks):
  MeshGeneratorStruct(mesh),_decks(decks)  {}


  /**
   * distructor
   */
  ~MeshGeneratorQuad4()
  {
    if(point_array3d)
    {
      for(unsigned int j=0;j<IY;j++)
        delete [] point_array3d[0][j];
      delete [] point_array3d[0];
      delete [] point_array3d;
    }
  }

  /**
   * create mesh
   */
  int do_mesh();


  /**
   * refine existing mesh
   */
  int do_refine(MeshRefinement &);

  /**
   * return megh generator magic munber
   */
  unsigned int magic_num()
  { return static_cast<unsigned int>(975); }

private:
  /**
   * Reference to the input decks
   */
  Parser::InputParser & _decks;

  /**
   * the data structure for reserved boundary edge in xy plane
   */
  std::map<SkeletonEdge, int, lt_edge>        edge_table;

 /**
   * the user specified face information
   */
  std::vector<SkeletonFace>                   face_array1d;

  /**
   * the user specified region information
   */
  std::vector<SkeletonRegion2D>               region_array1d;


  /**
   * build the rectangle Skeleton mesh
   */
  void build_rectangle_mesh();

  /**
   * read the region information from card "REGION", and fills into region_array1d
   */
  int  set_region(const Parser::Card &c);

  /**
   * read the boundary face information from card "FACE", and fills into face_array1d
   */
  int  set_face(const Parser::Card &c);

  /**
   * set region boundary with label "Neumann", this may be overwrite
   * by make_face() function
   */
  int  make_region_boundary();

  /**
   * set user defined boundary face, user should give a label to these faces
   */
  int  make_face(unsigned int ixmin,unsigned int ixmax,
                 unsigned int iymin,unsigned int iymax,
                 const std::string &label);

  /**
   * get the region index of a quad4 cell by its bound
   */
  int  get_cell_region_id(unsigned int ixmin, unsigned int ixmax,
                          unsigned int iymin, unsigned int iymax);

};


#endif // #define __mesh_generation_quad4_h__
