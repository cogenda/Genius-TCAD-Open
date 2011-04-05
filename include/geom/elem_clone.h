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


#ifndef __elem_clone_h__
#define __elem_clone_h__


// ------------------------------------------------------------
//ElemClone class definition
template <class ElemType>
class ElemClone : public ElemType
{
public:

  /**
   * Constructor.  Creates a ElemClone from an element.
   */
  ElemClone (const Elem* parent=NULL): ElemType(const_cast<Elem*>(parent))
  {}

  /**
   * Destructor, delete my nodes
   */
  ~ElemClone()
  {
    for(unsigned int n=0; n<ElemType::n_nodes(); ++n)
      delete ElemType::_nodes[n];
  }

};


#endif

