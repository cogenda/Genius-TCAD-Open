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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "dfise_grid.h"

//---------------------------------------------------
// functions for class GRID
namespace DFISE
{

  GRID::GRID()
  {
    for(int i=0; i<3; ++i) translate[i] = 0;
    for(int i=0; i<9; ++i) transform[i] = 0;
    transform[0] = transform[4] = transform[8] = 1.0;
  }

  void GRID::clear()
  {
    for(int i=0; i<3; ++i) translate[i] = 0;
    for(int i=0; i<9; ++i) transform[i] = 0;
    transform[0] = transform[4] = transform[8] = 1.0;
    Vertices.clear();
    Edges.clear();
    Faces.clear();
    Locations.clear();
    Elements.clear();
    region_elements.clear();
    _edge_set.clear();
  }

  std::vector<int> GRID::get_face_nodes(int face_index)
  {
    bool face_inverse = face_index < 0;
    face_index = face_index < 0 ? -face_index-1 : face_index;
    assert(static_cast<unsigned int>(face_index) < Faces.size());

    std::vector<int> nodes;

    const std::vector<int> & face = Faces[face_index];
    for(unsigned int n=0; n<face.size(); ++n)
    {
      int edge_index = face[n];
      bool edge_inverse = edge_index < 0;
      edge_index = edge_index < 0 ? -edge_index-1 : edge_index;
      if(!edge_inverse)
      {
        if(!nodes.empty() && nodes.back() == Edges[edge_index].second)
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
        else
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
      }
      else
      {
        if(!nodes.empty() && nodes.back() == Edges[edge_index].first)
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
        else
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
      }
    }
    if(nodes.front() != nodes.back() ) std::swap(nodes[0], nodes[1]);

    std::unique(nodes.begin(), nodes.end());
    nodes.resize(face.size());

    if(face_inverse)
    {
      nodes.push_back(nodes.front());
      std::reverse(nodes.begin(), nodes.end());
      nodes.resize(face.size());
    }

    return nodes;
  }


  void GRID::build_node(Element & elem)
  {
    switch(elem.elem_code)
    {
      case 1:  build_edge2_node(elem); break;
      case 2:  build_tri3_node(elem);  break;
      case 3:  build_quad4_node(elem); break;
      case 5:  build_tet4_node(elem);  break;
      case 6:  build_pyramid5_node(elem);  break;
      case 7:  build_prism6_node(elem);    break;
      case 8:  build_hex8_node(elem);      break;
    }
  }


  void GRID::build_edge2_node(Element & elem)
  {
    std::vector<int> nodes;
    nodes.push_back(elem.faces[0]);
    nodes.push_back(elem.faces[1]);
    elem.vertices = nodes;
    assert(elem.vertices.size()==2);
  }

  void GRID::build_tri3_node(Element & elem)
  {
    std::vector<int> nodes;
    for(unsigned int n=0; n<elem.faces.size(); ++n)
    {
      int edge_index = elem.faces[n];
      bool inverse = edge_index < 0;
      edge_index = edge_index < 0 ? -edge_index-1 : edge_index;
      if(!inverse)
      {
        if(!nodes.empty() && nodes.back() == Edges[edge_index].second)
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
        else
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
      }
      else
      {
        if(!nodes.empty() && nodes.back() == Edges[edge_index].first)
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
        else
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
      }
    }

    if(nodes.front() != nodes.back() ) std::swap(nodes[0], nodes[1]);
    std::unique(nodes.begin(), nodes.end());
    nodes.resize(elem.faces.size());

    // check for ccw, 2D only
    if( dimension==2 && !is_ccw_2d(Vertices[nodes[0]], Vertices[nodes[1]], Vertices[nodes[2]]) )
      std::reverse(nodes.begin(), nodes.end());

    elem.vertices = nodes;
  }


  void GRID::build_quad4_node(Element & elem)
  {
    std::vector<int> nodes;
    for(unsigned int n=0; n<elem.faces.size(); ++n)
    {
      int edge_index = elem.faces[n];
      bool inverse = edge_index < 0;
      edge_index = edge_index < 0 ? -edge_index-1 : edge_index;
      if(!inverse)
      {
        if(!nodes.empty() && nodes.back() == Edges[edge_index].second)
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
        else
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
      }
      else
      {
        if(!nodes.empty() && nodes.back() == Edges[edge_index].first)
        {
          nodes.push_back(Edges[edge_index].first);
          nodes.push_back(Edges[edge_index].second);
        }
        else
        {
          nodes.push_back(Edges[edge_index].second);
          nodes.push_back(Edges[edge_index].first);
        }
      }
    }

    if(nodes.front() != nodes.back() ) std::swap(nodes[0], nodes[1]);
    std::unique(nodes.begin(), nodes.end());
    nodes.resize(elem.faces.size());

    // check for ccw, 2D only
    if( dimension==2 && !is_ccw_2d(Vertices[nodes[0]], Vertices[nodes[1]], Vertices[nodes[2]]) )
      std::reverse(nodes.begin(), nodes.end());

    elem.vertices = nodes;
  }

  void GRID::build_tet4_node(Element & elem)
  {
    std::vector<int> face3_nodes = get_face_nodes(elem.faces[3]);
    std::vector<int> face1_nodes = get_face_nodes(elem.faces[1]);

    int p4=-1;
    for(unsigned int m=0; m<face1_nodes.size(); ++m)
    {
      if( std::find(face3_nodes.begin(), face3_nodes.end(), face1_nodes[m])==face3_nodes.end() )
      { p4 = face1_nodes[m]; break;}
    }
    assert(p4>=0);

    if(!is_above(Vertices[face3_nodes[0]], Vertices[face3_nodes[1]], Vertices[face3_nodes[2]], Vertices[p4]))
    {
      std::reverse(face3_nodes.begin(), face3_nodes.end());
    }

    elem.vertices = face3_nodes;
    elem.vertices.push_back(p4);
  }

  void GRID::build_pyramid5_node(Element & elem)
  {
    std::vector<int> face4_nodes = get_face_nodes(elem.faces[4]);
    std::vector<int> face1_nodes = get_face_nodes(elem.faces[1]);

    int p5=-1;
    for(unsigned int m=0; m<face1_nodes.size(); ++m)
    {
      if( std::find(face4_nodes.begin(), face4_nodes.end(), face1_nodes[m]) == face4_nodes.end() )
      {  p5 = face1_nodes[m];break; }
    }
    assert(p5>=0);

    if(!is_above(Vertices[face4_nodes[0]], Vertices[face4_nodes[1]], Vertices[face4_nodes[2]], Vertices[p5]))
    {
      std::reverse(face4_nodes.begin(), face4_nodes.end());
    }

    elem.vertices = face4_nodes;
    elem.vertices.push_back(p5);
  }

  void GRID::build_prism6_node(Element & elem)
  {
    std::vector<int> face3_nodes = get_face_nodes(elem.faces[3]);
    std::vector<int> face4_nodes = get_face_nodes(elem.faces[4]);

    if(!is_above(Vertices[face3_nodes[0]], Vertices[face3_nodes[1]], Vertices[face3_nodes[2]], Vertices[face4_nodes[0]]))
    {
      std::reverse(face3_nodes.begin(), face3_nodes.end());
    }
    elem.vertices = face3_nodes;

    for(unsigned int m=0; m<face3_nodes.size(); ++m)
    {
      for(unsigned int n=0; n<face4_nodes.size(); ++n)
        if( _edge_set.find(std::make_pair(face3_nodes[m], face4_nodes[n])) != _edge_set.end())
        {
          elem.vertices.push_back(face4_nodes[n]);
        }
    }
    assert(elem.vertices.size()==6);
  }

  void GRID::build_hex8_node(Element & elem)
  {
    std::vector<int> face4_nodes = get_face_nodes(elem.faces[4]);
    std::vector<int> face5_nodes = get_face_nodes(elem.faces[5]);

    if(!is_above(Vertices[face4_nodes[0]], Vertices[face4_nodes[1]], Vertices[face4_nodes[2]], Vertices[face5_nodes[0]]))
    {
      std::reverse(face4_nodes.begin(), face4_nodes.end());
    }
    elem.vertices = face4_nodes;

    for(unsigned int m=0; m<face4_nodes.size(); ++m)
    {
      for(unsigned int n=0; n<face5_nodes.size(); ++n)
        if( _edge_set.find(std::make_pair(face4_nodes[m], face5_nodes[n])) != _edge_set.end())
        {
          elem.vertices.push_back(face5_nodes[n]);
        }
    }
    assert(elem.vertices.size()==8);
  }


  void GRID::print( std::ostream & out ) const
  {
    out<<"Data {" << std::endl;
    out<< std::scientific << std::setprecision(15);

    // print CoordSystem
    out<<"  CoordSystem {" << std::endl;
    out<<"    translate = [  " << translate[0] << " " << translate[1] << " " << translate[2] <<"]" << std::endl;
    out<<"    transform = [  " << transform[0] << " " << transform[1] << " " << transform[2] << " "
    << transform[3] << " " << transform[4] << " " << transform[5] << " "
    << transform[6] << " " << transform[7] << " " << transform[8] <<"]" << std::endl;
    out<<"  }" << std::endl;
    out << std::endl;

    // print vertices
    out<<"  Vertices (" << Vertices.size() << ") {" <<std::endl;
    if(dimension==2)
    {
      for(unsigned int n=0; n<Vertices.size(); n++)
        out<< "    " << Vertices[n].coords[0] << " " << Vertices[n].coords[1] << std::endl;
    }
    if(dimension==3)
    {
      for(unsigned int n=0; n<Vertices.size(); n++)
        out<< "    " << Vertices[n].coords[0] << " " << Vertices[n].coords[1] << " " << Vertices[n].coords[2] << std::endl;
    }
    out<<"  }" << std::endl;
    out << std::endl;

    // print edges
    out<<"  Edges (" << Edges.size() << ") {" << std::endl;
    for(unsigned int n=0; n<Edges.size(); n++)
      out<< "    " << Edges[n].first << " " << Edges[n].second << std::endl;
    out<<"  }" << std::endl;
    out << std::endl;

    if(dimension==3)
    {
      out<<"  Faces (" << Faces.size() << ") {" << std::endl;
      for(unsigned int n=0; n<Faces.size(); n++)
      {
        out<< "    " << Faces[n].size();
        for(unsigned int e=0; e<Faces[n].size(); ++e)
          out<< " " << Faces[n][e] ;
        out << std::endl;
      }
      out<<"  }" << std::endl;
      out << std::endl;
    }

    // print location
    out<<"  Locations (" << Locations.size() << ") {" << std::endl;
    {
      int count=1;
      out<<"      ";
      for(unsigned int l=0; l<Locations.size(); ++l, ++count)
      {
        out<<Locations[l];
        if(count%10==0)
        {
          out<<std::endl;
          out<<"      ";
        }
      }
    }
    out<<"  }" << std::endl;
    out << std::endl;

    // print element
    out<<"  Elements (" << Elements.size() << ") {" << std::endl;
    for(unsigned int n=0; n<Elements.size(); n++)
    {
      const Element & elem = Elements[n];
      out<< "    " << elem.elem_code;
      for(unsigned int f=0; f<elem.faces.size(); ++f)
        out<< " " << elem.faces[f];
      out << std::endl;
    }
    out<<"  }" << std::endl;
    out << std::endl;

    // print region
    for(unsigned int r=0; r<regions.size(); ++r)
    {
      out<<"  Region (" << " \"" << regions[r]<< "\" " << ") {" << std::endl;
      out<<"    material = " << materials[r] << std::endl;
      const std::vector<unsigned int> & element = region_elements[r];
      out<<"    Elements (" << element.size() <<") {" << std::endl;
      int count=1;
      out<<"      ";
      for(unsigned int e=0; e<element.size(); ++e, ++count)
      {
        out << element[e];
        if(count%10==0)
        {
          out<<std::endl;
          out<<"      ";
        }
        else
          out<< " ";
      }
      out<<"    }"<< std::endl;
      out<<"  }"<< std::endl;
      out << std::endl;
    }
    out<<"}"<< std::endl;
  }
}

