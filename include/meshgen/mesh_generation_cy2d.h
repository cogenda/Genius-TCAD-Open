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



#ifndef __mesh_generation_cy2d_h__
#define __mesh_generation_cy2d_h__

#include "mesh_generation_cy.h"

/**
 * the simple cy2d mesh generator
 */
class MeshGeneratorCylinder2D : public MeshGeneratorCylinder
{
public:
  /**
   * constructor
   */
  MeshGeneratorCylinder2D(MeshBase& mesh, Parser::InputParser & decks):
  MeshGeneratorCylinder(mesh),_decks(decks)  {}


  /**
   * distructor
   */
  ~MeshGeneratorCylinder2D()
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
  { return static_cast<unsigned int>(220); }

private:
  
  /**
   * Reference to the input decks
   */
  Parser::InputParser & _decks;


  /**
   * the user specified segment information
   */ 
  std::vector<SkeletonSector>          segment_r_array1d;
  
   /**
   * the user specified segment information
   */ 
  std::vector<SkeletonSector>          segment_theta_array1d;
  
  /**
   * map segment label to id
   */  
  std::map<std::string, int>                  segment_id_map;

  /**
   * the user specified region information
   */
  std::vector<SkeletonRegion2D>               region_array1d;
  

  /**
   * read the region information from card "REGION", and fills into region_array1d
   */
  int  set_region(const Parser::Card &c);

  /**
   * read the boundary face information from card "FACE", and fills into face_array1d
   */
  int  set_face(const Parser::Card &c);


  /**
   * set user defined boundary segment by mid point of edge
   */
  int  get_segment_id_theta(double rmid, const Point &mid);
  
  /**
   * set user defined boundary segment by mid point of edge
   */
  int  get_segment_id_r(const Point &p1, const Point &p2);

  /**
   * get the region index of a cell by its inner point
   */
  int  get_cell_region_id(const Point &in);

};


#endif // #define __mesh_generation_cy2d_h__
