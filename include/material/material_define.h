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

//  $Id: material_define.h,v 1.6 2008/07/09 05:58:16 gdiso Exp $

#ifndef __material_define_h_
#define __material_define_h_

#include <string>
#include <vector>

namespace Material {

/**
 * the material type Genius supproted now
 */
enum  MaterialType
{
 Semiconductor                  ,
 SingleCompoundSemiconductor    ,
 ComplexCompoundSemiconductor   ,
 Insulator                      ,
 Conductor                      ,
 Resistance                     ,
 PML                            ,
 Vacuum                         ,
 INVALID_MATERIAL_TYPE            // should always be last
};

//----------------------------------------------------------
// call this function before any other function!
//----------------------------------------------------------
/**
 * load the material definitions from file
 */
extern void init_material_define(const std::string &fname);


//----------------------------------------------------------
// the following functions are used to judge material type
//----------------------------------------------------------


/**
 * @return true if mat_name match any semiconductor material name
 */
extern bool IsSemiconductor(const std::string & mat_name);

/**
 * @return true if mat_name match any single compound semiconductor material name
 */
extern bool IsSingleCompSemiconductor(const std::string & mat_name);


/**
 * @return true if mat_name match any complex compound semiconductor material name
 */
extern bool IsComplexCompSemiconductor(const std::string & mat_name);

/**
 * @return true if mat_name match a material name of insulator
 */
extern bool IsInsulator(const std::string & mat_name);

/**
 * @return true if mat_name match a material name of conductor
 */
extern bool IsConductor(const std::string & mat_name);

/**
 * @return true if mat_name match a material name of resistance material
 */
extern bool IsResistance(const std::string & mat_name);


/**
 * @return true if mat_name is Vacuum
 */
extern bool IsVacuum(const std::string & mat_name);

/**
 * @return true if mat_name is MPL
 */
extern bool IsPML(const std::string & mat_name);

/**
 * @return enum MaterialType by material name
 */
extern MaterialType material_type(const std::string & mat_name);

/**
 * @return a weight factor of material, the semiconductor weight is 3
 * others 1
 */
extern int material_weight(const std::string & mat_name);

//----------------------------------------------------------
// since some matrial has alias ( such as "Ox" and "SiO2" )
// format them to unique name
//----------------------------------------------------------


/**
 * convert material name to an unique name
 */
extern std::string FormatMaterialString(const std::string & mat_name);



//----------------------------------------------------------
// convert material string <--> id
//----------------------------------------------------------

extern std::string get_material_by_id(unsigned int id);
extern unsigned int get_id_by_material(const std::string & mat_name);

extern std::vector<unsigned int> get_material_ids();


//----------------------------------------------------------
// get material color by id
//----------------------------------------------------------

extern void get_material_color(unsigned int id, double &r, double &g, double &b, double &alpha);

}

#endif
