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

//  $Id: mesh_generation_base.h,v 1.7 2008/07/09 05:58:16 gdiso Exp $



#ifndef __mesh_generation_base_h__
#define __mesh_generation_base_h__

#include "mesh_base.h"
#include "mesh_modification.h"
#include "mesh_refinement.h"


class MeshGeneratorBase
{
public:

  /**
   * Constructor.
   */
  MeshGeneratorBase (MeshBase& mesh):_mesh(mesh) {}

  /**
   * destructor, do nothing
   */
  virtual ~MeshGeneratorBase() {}
  /** 
   * the virtual do_mesh method
   */
  virtual int do_mesh()=0;

  /**
   * refine existing mesh 
   */
  virtual int do_refine(MeshRefinement & )=0;
    
  /**
   * @return megh generator magic munber
   */
  virtual unsigned int magic_num()=0;
  
protected:
 /**
  * Reference to the mesh.
  */
  MeshBase & _mesh;
  
};




#endif // #define __mesh_generation_base_h__
