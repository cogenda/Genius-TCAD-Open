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

//  $Id: material_define.cc,v 1.2 2008/07/09 05:58:16 gdiso Exp $

#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <numeric>
#include <cstdlib>

#include "genius_common.h"
#include "material_define.h"

namespace Material
{

  struct Material
  {
    std::string   name;
    std::string   property;
    std::vector<std::string> alias;
    unsigned int  id;
    unsigned int  color;
  };

  std::vector<Material> materials;

  // read predefined material information from $GENIUS_DIR/lib/material.def
  static void _read_predefined_material_information(const std::string &fname)
  {
    std::ifstream in(fname.c_str());
    genius_assert(in.good());

    std::string content;

    //first, we read whole file into a std::string
    while(!in.eof())
    {
      std::string line;
      std::getline(in, line);
      // change cr/newline chars to blank
      {
        std::string filt_elems("\r\n");
        std::string::size_type pos = 0;
        while (( pos = line.find_first_of( filt_elems, pos )) != std::string::npos )
        {
          line.replace(pos, 1, " ");
        }
      }
      // skip the line begin with '#'
      if(line.find('#')==0) continue;
      content += line;
    }
    in.close();

    // add blank before and after "{}=", then we can use stringstream to get each token
    {
      std::string filt_elems("{}=");
      std::string::size_type pos = 0;
      while (( pos = content.find_first_of( filt_elems, pos )) != std::string::npos )
      {
        content.insert(pos,   " ");
        content.insert(pos+2, " ");
        pos+=3;
      }
    }


    std::stringstream ss(content);
    while(1)
    {
      std::string token1, token2, token3;

      Material material;
      // the beginning of each material define in file material.def
      ss >> material.name;
      if(ss.fail()) break;
      // should have a '{'
      ss >> token1;
      genius_assert(token1=="{");

      while(1)
      {
        // do until we meet a '}'
        ss >> token1;
        if(token1=="}") break;

        // token2 should be '='
        ss >> token2; genius_assert(token2=="=");
        ss >> token3;

        if(token1=="property")
          material.property = token3;
        else if(token1=="alias")
          material.alias.push_back(token3);
        else if(token1=="color")
        {
          std::stringstream to_uint(token3);
          to_uint << std::hex;
          to_uint >> material.color;
        }
      }
      materials.push_back(material);
    }
  }


  std::map<const std::string, MaterialType> material_type_string_to_enum;
  static void _build_material_type_string_to_enum_map()
  {
    material_type_string_to_enum["Semiconductor"               ] = Semiconductor;
    material_type_string_to_enum["SingleCompoundSemiconductor" ] = SingleCompoundSemiconductor;
    material_type_string_to_enum["ComplexCompoundSemiconductor"] = ComplexCompoundSemiconductor;
    material_type_string_to_enum["Conductor"                   ] = Conductor;
    material_type_string_to_enum["Insulator"                   ] = Insulator;
    material_type_string_to_enum["Vacuum"                      ] = Vacuum;
    material_type_string_to_enum["PML"                         ] = PML;
  }


  std::map<const std::string, MaterialType> material_name_to_material_type;
  static void _init_material_name_to_material_type_map()
  {
    for(unsigned int n=0; n<materials.size(); ++n)
    {
      std::string property = materials[n].property;
      material_name_to_material_type[materials[n].name] = material_type_string_to_enum[property];
      for(unsigned int k=0; k<materials[n].alias.size(); ++k)
        material_name_to_material_type[materials[n].alias[k]] = material_type_string_to_enum[property];
    }
  }


  std::map<const std::string, std::string> material_name_conversion;
  static void _init_material_name_conversion_map()
  {
    for(unsigned int n=0; n<materials.size(); ++n)
    {
      material_name_conversion[materials[n].name] = materials[n].name;
      for(unsigned int k=0; k<materials[n].alias.size(); ++k)
        material_name_conversion[materials[n].alias[k]] = materials[n].name;
    }
  }


  std::map<std::string, unsigned int> material_name_to_material_id;
  std::map<unsigned int, std::string> material_id_to_material_name;
  static void _init_material_name_to_material_id_map()
  {
     for(unsigned int n=0; n<materials.size(); ++n)
     {
        std::string mat_name = materials[n].name;

        // build a unique id for the material
        unsigned int sum = std::accumulate(mat_name.begin(), mat_name.end(), 0);
        MaterialType material_type = material_name_to_material_type[mat_name];
        unsigned int id = 10000*static_cast<unsigned int>(material_type) + sum;
        genius_assert(material_id_to_material_name.find(id)==material_id_to_material_name.end());
        materials[n].id = id;
        material_name_to_material_id[mat_name] = id;
        material_id_to_material_name[id] = mat_name;
     }
  }

  //----------------------------------------------------------
  // setup material information
  //----------------------------------------------------------

  void init_material_define(const std::string &fname)
  {
      _read_predefined_material_information(fname);
      _build_material_type_string_to_enum_map();
      _init_material_name_to_material_type_map();
      _init_material_name_conversion_map();
      _init_material_name_to_material_id_map();
  }

  //----------------------------------------------------------
  // the following functions are used to judge material type
  //----------------------------------------------------------



  bool IsSemiconductor(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());

    if ( material_name_to_material_type[mat_name] == Semiconductor                  ||
         material_name_to_material_type[mat_name] == SingleCompoundSemiconductor    ||
         material_name_to_material_type[mat_name] == ComplexCompoundSemiconductor
       )
      return true;

    return false;
  }



  bool IsSingleCompSemiconductor(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());

    if ( material_name_to_material_type[mat_name] == SingleCompoundSemiconductor )
      return true;

    return false;
  }


  bool IsComplexCompSemiconductor(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());

    if ( material_name_to_material_type[mat_name] == ComplexCompoundSemiconductor )
      return true;

    return false;
  }


  bool IsInsulator(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());

    if ( material_name_to_material_type[mat_name] == Insulator )
      return true;

    return false;
  }



  bool IsConductor(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());

    if ( material_name_to_material_type[mat_name] == Conductor )
      return true;

    return false;
  }



  bool IsVacuum(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());

    if ( material_name_to_material_type[mat_name] == Vacuum )
      return true;

    return false;
  }



  bool IsPML(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());

    if ( material_name_to_material_type[mat_name] == PML )
      return true;

    return false;
  }



  MaterialType material_type(const std::string & mat_name)
  {
    genius_assert(material_name_to_material_type.size());
    return material_name_to_material_type[mat_name];
  }




  int material_weight(const std::string & mat_name)
  {
    if( IsSemiconductor(mat_name) ) return 6;
    return 1;
  }


  //----------------------------------------------------------
  // since some matrial has alias ( such as "Ox" and "SiO2" )
  // format them to unique name
  //----------------------------------------------------------

  std::string FormatMaterialString(const std::string & mat_name)
  {
    genius_assert(material_name_conversion.size());

    if( material_name_conversion.find(mat_name) != material_name_conversion.end() )
      return material_name_conversion[mat_name];
    return mat_name;
  }

  //----------------------------------------------------------
  // convert material string <--> id
  //----------------------------------------------------------

  std::string get_material_by_id(unsigned int id)
  {
    genius_assert(material_id_to_material_name.find(id)!=material_id_to_material_name.end());
    return material_id_to_material_name[id];
  }

  unsigned int get_id_by_material(const std::string & mat_name)
  {
    std::string unique_name = FormatMaterialString(mat_name);
    genius_assert(material_name_to_material_id.find(unique_name)!=material_name_to_material_id.end());
    return material_name_to_material_id[unique_name];
  }


  std::vector<unsigned int> get_material_ids()
  {
    std::vector<unsigned int> material_ids;
    for(unsigned int n=0; n<materials.size(); ++n)
      material_ids.push_back(materials[n].id);
    return material_ids;
  }

//----------------------------------------------------------
// get material color by id
//----------------------------------------------------------

  void get_material_color(unsigned int id, double &r, double &g, double &b, double &alpha)
  {
     //default value
     r = 1.0;
     g = 1.0;
     b = 1.0;
     alpha = 1.0;

     for(unsigned int n=0; n<materials.size(); ++n)
       if(materials[n].id == id)
       {
         unsigned int color = materials[n].color;
	 r = ( (color & 0xff000000) >> 24)/255.0;
	 g = ( (color & 0x00ff0000) >> 16)/255.0;
	 b = ( (color & 0x0000ff00) >> 8 )/255.0;
	 alpha = (color & 0x000000ff)/255.0;
       }
  }

}


