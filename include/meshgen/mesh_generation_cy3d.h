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



#ifndef __mesh_generation_cy3d_h__
#define __mesh_generation_cy3d_h__

#include "mesh_generation_cy.h"

/**
 * the simple cy3d mesh generator
 */
class MeshGeneratorCylinder3D : public MeshGeneratorCylinder
{
public:
  /**
   * constructor
   */
  MeshGeneratorCylinder3D(MeshBase& mesh, Parser::InputParser & decks):
  MeshGeneratorCylinder(mesh),_decks(decks)  {}


  /**
   * distructor
   */
  ~MeshGeneratorCylinder3D()
  {
  }

  /**
   * create mesh
   */
  int do_mesh();


  /**
   * refine existing mesh
   */
  int do_refine(MeshRefinement &) {}

  /**
   * return megh generator magic munber
   */
  unsigned int magic_num()
  { return static_cast<unsigned int>(2200); }

private:

  /**
   * Reference to the input decks
   */
  Parser::InputParser & _decks;


  /**
   * the user specified face information
   */
  std::vector<SkeletonSector>          face_r_array1d;

  /**
   * the user specified face information
   */
  std::vector<SkeletonSector>          face_theta_array1d;

  /**
   * the user specified face information
   */
  std::vector<SkeletonSector>          face_z_array1d;

  /**
   * map segment label to id
   */
  std::map<std::string, int>                  face_id_map;

  /**
   * the user specified region information
   */
  std::vector<SkeletonRegion3D>               region_array1d;


  /**
   * read the region information from card "REGION", and fills into region_array1d
   */
  int  set_region(const Parser::Card &c);

  /**
   * read the boundary face information from card "FACE", and fills into face_array1d
   */
  int  set_face(const Parser::Card &c);


  /**
   * set user defined boundary face by mid point of edge
   */
  int  get_face_id_theta(double rmid, const Point &mid);

  /**
   * set user defined boundary face by mid point of edge
   */
  int  get_face_id_r(const Point &p1, const Point &p2, const Point &p3, const Point &p4);

  /**
   * set user defined boundary face by mid point of edge
   */
  int  get_face_id_z(double rmid, const Point &mid);

  /**
   * get the region index of a cell by its inner point
   */
  int  get_cell_region_id(const Point &in);

};


#endif // #define __mesh_generation_cy3d_h__
