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

#ifndef __dose_rate_h__
#define __dose_rate_h__

#include <vector>

#include "octree.h"

class SimulationSystem;
class MeshBase;

class DoseRate
{
  public:

    DoseRate(const SimulationSystem & system );

    ~DoseRate();

    /**
     * set minimal leaf size
     */
    void set_min_distance(double x) { min_leaf = x; }


    /**
     * refine the octree;
     */
    void refine();


    /**
     * calculate energy deposite of one track 
     */
    void energy_deposite(const Point &p1, const Point &p2, double e);

    /**
     * parallel sync
     */
    void sync_energy_deposite();

    /**
     * clear all energy deposite
     */
    void clear_energy_deposite();

    /**
     * get energy deposite density (J/kg) of a given point
     */
    double energy_deposite_density( const Point & ) const;

    /**
     * total energy on this octree
     */
    double total_energy() const;

    /**
     * debug only
     */
    void export_vtk(const std::string &file);



    /**
     * add a particle endpoint
     */
    void particle_endpoint(const Point &p);

    /**
     * FIXME no nothing at present
     * parallel sync
     */
    void sync_particle_endpoint() {}

    /**
     * clear particle endpoint list
     */
    void clear_particle_endpoint();

  protected:

    const SimulationSystem & _system;

    const MeshBase & _mesh;

    /**
     * density of each region
     */
    static std::vector<double> density;

    /**
     * build mesh constrain
     */
    static void build_region_density(const SimulationSystem & system) ;

    /**
     * mesh constrain for each region
     */
    static std::vector<double> constrain;

    /**
     * build mesh constrain
     */
    static void build_mesh_constrain(const MeshBase & mesh) ;

    /**
     * minimal leaf size
     */
    static double min_leaf;

    /**
     * compact elem
     */
    struct CElem
    {
      unsigned int id; //elem id
      Point center;    // elem center
      unsigned int subdomain; // elem subdomain
    };

    void mesh_elem(std::vector<CElem> &) const;


    class OcTreeDataDoseRate : public OcTreeDataBase
    {
      public:
        OcTreeDataDoseRate()
        : electron_energy(0.0), weighted_density(0.0)
        {}
        void calculate_density();

        std::vector<CElem> elems;
        std::vector<Point> electron_endpoint;
        double electron_energy;
        double weighted_density;

        std::vector<OcTreeDataBase *> subdivide(const OcTreeNode & node);
        bool refine(const OcTreeNode & node, const std::vector<OcTreeNode> &neighbor) const;
        double value(const std::string &) const;
    };

    /**
     * geometry octree
     */
    OcTree * _octree;
};

#endif

