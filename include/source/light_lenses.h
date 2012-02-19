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
/*  Author: xianghua zhang   zhangxih@163.com                                   */
/*                                                                              */
/********************************************************************************/

#ifndef __light_lenses_h__
#define __light_lenses_h__

#include <vector>
#include <map>

#include "sphere.h"

namespace Parser{
  class InputParser;
}
class LightThread;

class LightLenses
{
public:
  /**
   * construct pre-defined lenses
   */
  LightLenses(Parser::InputParser & decks);

  /**
   * free pre-defined lenses
   */
  ~LightLenses();

  /**
   * set active lenses to be used
   */
  void set_lenses(const std::vector<std::string> & lens);

  /**
   * @return true when no lens is active
   */
  bool empty() const {return _effect_lenses.empty();}

  /**
   * @return the bounding sphere of active lenses
   */
  Sphere bounding_sphere() const;

  /**
   * apply lens to light thread
   */
  LightThread * operator << (LightThread *) const;

private:

  struct Lens
  {
    std::string id;
    Point center;
    Point norm;
    double radius;
    double A,B,C,D;
  };

  std::vector<std::string> _effect_lenses;
  std::map<std::string, Lens *> _lenses;
};

#endif
