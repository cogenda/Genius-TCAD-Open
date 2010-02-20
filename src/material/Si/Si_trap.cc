/*****************************************************************************/
/*                                                                           */
/*              8888888         88888888         88888888                    */
/*            8                8                8                            */
/*           8                 8                8                            */
/*           8                  88888888         88888888                    */
/*           8      8888                8                8                   */
/*            8       8                 8                8                   */
/*              888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS material database Version 0.4                                        */
/*  Last update: Feb 17, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Silicon

#include "PMI.h"

#include <algorithm>
#include <vector>
#include <map>
#include <iomanip>
#include <ctype.h>

class GSS_Si_Trap_Default : public PMIS_Trap
{
private:
  // {{{ class TrapSpec
  class TrapSpec
  {
    public:
      TrapChargeType charge_type;               // Neutral, Acceptor or Donor
      TrapType type;                            // Bulk or Interface


      /******** for bulk traps ***********************************************/
      // the doping profile name that specifies the spatial distribution of bulk traps
      std::string profile_name;
      // prefactor of trap density
      PetscScalar prefactor;


      /******** for interface traps ******************************************/
      // the interface with interface trap
      std::string interface_name;
      // interface trap density
      PetscScalar interface_density;

      PetscScalar E_t;      // trap energy with respect to midgap energy
      PetscScalar sigma_n;  // electron capture cross-section
      PetscScalar sigma_p;  // hole capture cross-section
      PetscScalar g_n;      // degeneracy factor with conduction band
      PetscScalar g_p;      // degeneracy factor with valance band

    public:
      TrapSpec()
      {
        profile_name = "";
        interface_name= "";
        prefactor = 0;
        interface_density = 0;
        charge_type = Neutral;
        type = Bulk;
        E_t  = 0;
        sigma_n = 0;
        sigma_p = 0;
        g_n = 1.0;
        g_p = 1.0;
      }
      TrapSpec(const std::string &name, const PetscScalar &pf, const PetscScalar &dit, const TrapChargeType ct, const TrapType t, const PetscScalar &Et, const PetscScalar &sn, const PetscScalar &sp, const PetscScalar &gn, const PetscScalar &gp)
      {
        if (t==Bulk)
        {
          // bulk trap
          profile_name = name;
          interface_name = "";
          prefactor = pf;
        }
        else
        {
          // interface trap
          profile_name = "";
          interface_name = name;
          interface_density = dit;
        }
        charge_type    = ct;
        type = t;
        E_t     = Et;
        sigma_n = sn;
        sigma_p = sp;
        g_n     = gn;
        g_p     = gp;
      }
  };
  // }}}

  std::vector<TrapSpec> TrapSpecs;

  // {{{ class Trap
  class Trap
  {
    public:
      int spec;
      PetscScalar N_tt;     // total trap density

      // density of trapped electron in the current iteration
      PetscScalar n_t;
      AutoDScalar n_t_AD;

      // density of trapped electron (convergenced value) at the last time step
      // We use BDF2 discretization here, two previous values is needed.
      PetscScalar n_t_last;
      PetscScalar n_t_last_last;

      // timestamp of the last two steps
      PetscScalar clock_last;
      PetscScalar clock_last_last;

    public:
      Trap() : n_t_AD(0)
      {
        spec = 0;
        N_tt = 0;
        n_t  = 0;
        n_t_last = 0;
        n_t_last_last = 0;
        clock_last = 0;
        clock_last_last = 0;
      }
      Trap(const unsigned int spc, const PetscScalar &Ntt, const PetscScalar &nt) : n_t_AD(nt)
      {
        spec = spc;
        N_tt = Ntt;
        n_t  = nt;
        n_t_last = nt;
        n_t_last_last = nt;
        clock_last = 0;
        clock_last_last = 0;
      }
  };
  // }}}

  typedef std::map<TrapLocation, std::vector<Trap>, TrapLocationComp> TrapStore_t;

  TrapStore_t TrapStore;
  Trap dummy_trap;

  void 	Trap_Init()
  {
  }

  struct ToLower
  {
    char operator() (char c) const  { return tolower(c); }
  };

public:

//----------------------------------------------------------------
// constructor and destructor
public:

  GSS_Si_Trap_Default(const PMIS_Environment &env):PMIS_Trap(env)
  {
    PMI_Info = "This is the Default carrier trapping model of Silicon"; 
    Trap_Init();
  }

  // {{{ int AddBulkTrapSpec(const std::string &profile_name, const PetscScalar &prefactor, const TrapChargeType ct, const PetscScalar &Et, const PetscScalar &sn, const PetscScalar &sp, const PetscScalar &gn, const PetscScalar &gp)
  /**
   * declare a class of bulk traps as specified in the PMI command
   * But the traps are not actually created until the PMIs are initilized with init_node().
   */
  int AddBulkTrapSpec(const std::string &profile_name, const PetscScalar &prefactor, const TrapChargeType ct, const PetscScalar &Et, const PetscScalar &sn, const PetscScalar &sp, const PetscScalar &gn, const PetscScalar &gp)
  {
    TrapSpec s(profile_name, prefactor, 0, ct, Bulk, Et,sn,sp,gn,gp);
    TrapSpecs.push_back(s);
    return (TrapSpecs.size()-1);
  }
  // }}}

  // {{{ int AddInterfaceTrapSpec(const std::string &if_name, const PetscScalar &interface_density, const TrapChargeType ct, const PetscScalar &Et, const PetscScalar &sn, const PetscScalar &sp, const PetscScalar &gn, const PetscScalar &gp)
  /**
   * declare a class of interface traps as specified in the PMI command
   * But the traps are not actually created until the PMIs are initilized with init_bc_node().
   */
  int AddInterfaceTrapSpec(const std::string &if_name, const PetscScalar &interface_density, const TrapChargeType ct, const PetscScalar &Et, const PetscScalar &sn, const PetscScalar &sp, const PetscScalar &gn, const PetscScalar &gp)
  {
    TrapSpec s(if_name, 0, interface_density, ct, Interface, Et,sn,sp,gn,gp);
    TrapSpecs.push_back(s);
    return (TrapSpecs.size()-1);
  }
  // }}}

  // {{{ void AddTrap (const Point &point, int &spec, const PetscScalar &Ntt)
  /**
   * Append a trap to TrapStore
   *
   * @Point  the coordinate of the node, which is used as the identifier of the trap
   * @spec   the spec index of the trap
   * @Ntt    the concentration of the trap
   */
  void AddTrap (const Point &point, const unsigned int &spec, const PetscScalar &Ntt)
  {
    genius_assert(spec<TrapSpecs.size());

    // trap location is used as the key in TrapStore
    TrapLocation tloc = TrapLocation(point,TrapSpecs[spec].type);

    PetscScalar n_t;
    if (TrapSpecs[spec].charge_type == Acceptor)
      n_t = 0;        // Acceptors are assumed to be initially empty
    else if (TrapSpecs[spec].charge_type == Donor)
      n_t = Ntt;      // Donors are assumed to be initially occupied
    else
      n_t = 0.5*Ntt;

    Trap trap=Trap(spec, Ntt, n_t);   // create the trap

    TrapStore_t::iterator it;
    it = TrapStore.find(tloc);
    if (it != TrapStore.end())
    {
      // we already have some traps defined at this location, append the new trap
      it->second.push_back(trap);
    }
    else
    {
      std::pair<TrapStore_t::iterator, bool> ret;
      std::vector<Trap> vec_trap;

      // this is the first trap at the location.
      // Create an entry in the map first, which is itself a vector
      // then insert an trap to the vector
      ret = TrapStore.insert(std::pair< TrapLocation,std::vector<Trap> >(tloc,vec_trap));
      ret.first->second.push_back(trap);
    }
  }
  // }}}

  // {{{ PetscScalar Charge(const bool flag_bulk)
  /**
   * returns the electric charge density due to trapped charge at this node
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  PetscScalar Charge(const bool flag_bulk)
  {
    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      // yes there are some traps at this node
      PetscScalar Q=0;
      std::vector<Trap> & traps = it->second;

      // loop over all traps at this location
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        if (TrapSpecs[ptrap->spec].charge_type == Acceptor)
          Q = Q - ptrap->n_t;                   // Acceptor: neutral when empty, negatively charged when occupied
        else if (TrapSpecs[ptrap->spec].charge_type == Donor)
          Q = Q + (ptrap->N_tt - ptrap->n_t);   // Donor: positively charged when empty, neutral when occupied
      }
      return Q*e;
    }
    else
      return 0;
  }
  // }}}

  // {{{ AutoDScalar ChargeAD(const bool flag_bulk)
  /**
   * returns the partial derivatives of electric charge density
   * w.r.t. the local V,n,p
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  AutoDScalar ChargeAD(const bool flag_bulk)
  {
    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      AutoDScalar Q=0;
      std::vector<Trap> & traps = it->second;
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        if (TrapSpecs[ptrap->spec].charge_type == Acceptor)
          Q = Q - ptrap->n_t_AD;
        else if (TrapSpecs[ptrap->spec].charge_type == Donor)
          Q = Q + (ptrap->N_tt - ptrap->n_t_AD);
      }
      return Q*e;
    }
    else
      return 0;
  }
  // }}}

  // {{{ PetscScalar ElectronTrapRate(const bool flag_bulk, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl)
  /**
   * Calculates the electron trapping rate
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  PetscScalar ElectronTrapRate(const bool flag_bulk, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl)
  {
    PetscScalar theta_n = 1.0e7*cm/s * sqrt(Tl/300/K);     // electron thermal velocity

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      // yes we have some traps at this location
      PetscScalar Cn=0, En=0;
      std::vector<Trap> & traps = it->second;

      // loop over all traps at this location
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_n,g_n,E_t;
        sigma_n = TrapSpecs[ptrap->spec].sigma_n;
        g_n     = TrapSpecs[ptrap->spec].g_n;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        Cn += sigma_n * theta_n * n * (ptrap->N_tt - ptrap->n_t);   // capture rate
        En += sigma_n * theta_n / g_n * ni * ptrap->n_t * exp(E_t/kb/Tl);  // emission rate
      }
      return Cn-En;
    }
    else
    {
      return 0.;
    }
  }
  // }}}

  // {{{ AutoDScalar ElectronTrapRate(const bool flag_bulk, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tl)
  /**
   * Calculates the partial derivatives of electron trapping rate
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  AutoDScalar ElectronTrapRate(const bool flag_bulk, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tl)
  {
    AutoDScalar theta_n = 1.0e7*cm/s * sqrt(Tl/300/K);

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      AutoDScalar Cn,En;
      Cn=0; En=0;
      std::vector<Trap> & traps = it->second;
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_n,g_n,E_t;
        sigma_n = TrapSpecs[ptrap->spec].sigma_n;
        g_n     = TrapSpecs[ptrap->spec].g_n;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        Cn += sigma_n * theta_n * n * (ptrap->N_tt - ptrap->n_t_AD);
        En += sigma_n * theta_n / g_n * ni * ptrap->n_t_AD * exp(E_t/kb/Tl);
      }
      return Cn-En;
    }
    else
    {
      return 0.;
    }
  }
  // }}}

  // {{{ PetscScalar HoleTrapRate(const bool flag_bulk, const PetscScalar &p, const PetscScalar &ni, const PetscScalar &Tl)
  /**
   * Calculates the hole trapping rate
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  PetscScalar HoleTrapRate(const bool flag_bulk, const PetscScalar &p, const PetscScalar &ni, const PetscScalar &Tl)
  {
    PetscScalar theta_p = 1.0e7*cm/s * sqrt(Tl/300/K);     // hole thermal velocity

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      // yes we have some traps at this location
      PetscScalar Cp=0, Ep=0;
      std::vector<Trap> & traps = it->second;

      // loop over all traps at this location
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_p,g_p,E_t;
        sigma_p = TrapSpecs[ptrap->spec].sigma_p;
        g_p     = TrapSpecs[ptrap->spec].g_p;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        Cp += sigma_p * theta_p * p * ptrap->n_t;                                     // capture rate
        Ep += sigma_p * theta_p / g_p * ni * (ptrap->N_tt - ptrap->n_t) * exp(-E_t/kb/Tl);   // emission rate
      }
      return Cp-Ep;
    }
    else
    {
      return 0.;
    }
  }
  // }}}

  // {{{ AutoDScalar HoleTrapRate(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &ni, const AutoDScalar &Tl)
  /**
   *  Calculates the partial derivatives of hole trapping rate
   * one should call Calculate() to calculate the electron occupancy before calling this function
   */
  AutoDScalar HoleTrapRate(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &ni, const AutoDScalar &Tl)
  {
    AutoDScalar theta_p = 1.0e7*cm/s * sqrt(Tl/300/K);

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      AutoDScalar Cp=0, Ep=0;
      std::vector<Trap> & traps = it->second;
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_p,g_p,E_t;
        sigma_p = TrapSpecs[ptrap->spec].sigma_p;
        g_p     = TrapSpecs[ptrap->spec].g_p;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        Cp += sigma_p * theta_p * p * ptrap->n_t_AD;
        Ep += sigma_p * theta_p / g_p * ni * (ptrap->N_tt - ptrap->n_t_AD) * exp(-E_t/kb/Tl);
      }
      return Cp-Ep;
    }
    else
    {
      return 0.;
    }
  }
  // }}}

  // {{{ PetscScalar TrapHeat(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tp, const PetscScalar &Tn, const PetscScalar &Tl, const PetscScalar &EcEi, const PetscScalar &EiEv)
  PetscScalar TrapHeat(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tp, const PetscScalar &Tn, const PetscScalar &Tl, const PetscScalar &EcEi, const PetscScalar &EiEv)
  {
    PetscScalar theta_n = 1.0e7*cm/s * sqrt(Tl/300/K);     // electron thermal velocity
    PetscScalar theta_p = 1.0e7*cm/s * sqrt(Tl/300/K);     // hole thermal velocity

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      // yes we have some traps at this location
      std::vector<Trap> & traps = it->second;

      PetscScalar H=0;
      // loop over all traps at this location
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_n,sigma_p,g_n,g_p,E_t;
        sigma_n = TrapSpecs[ptrap->spec].sigma_n;
        sigma_p = TrapSpecs[ptrap->spec].sigma_p;
        g_n     = TrapSpecs[ptrap->spec].g_n;
        g_p     = TrapSpecs[ptrap->spec].g_p;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        PetscScalar Cn = sigma_n * theta_n * n * (ptrap->N_tt - ptrap->n_t);   // capture rate
        PetscScalar En = sigma_n * theta_n / g_n * ni * ptrap->n_t * exp(E_t/kb/Tl);  // emission rate

        PetscScalar Cp = sigma_p * theta_p * p * ptrap->n_t;                                     // capture rate
        PetscScalar Ep = sigma_p * theta_p / g_p * ni * (ptrap->N_tt - ptrap->n_t) * exp(-E_t/kb/Tl);   // emission rate

        H += (Cn-En) * (1.5*kb*Tn + EcEi);
        H += (Cp-Ep) * (1.5*kb*Tp + EiEv);
      }
      return H;
    }
    else
    {
      return 0.;
    }
  }
  // }}}

  // {{{ AutoDScalar TrapHeat(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tp, const AutoDScalar &Tn, const AutoDScalar &Tl, const AutoDScalar &EcEi, const AutoDScalar &EiEv)
  AutoDScalar TrapHeat(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tp, const AutoDScalar &Tn, const AutoDScalar &Tl, const AutoDScalar &EcEi, const AutoDScalar &EiEv)
  {
    AutoDScalar theta_n = 1.0e7*cm/s * sqrt(Tl/300/K);     // electron thermal velocity
    AutoDScalar theta_p = 1.0e7*cm/s * sqrt(Tl/300/K);     // hole thermal velocity

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      // yes we have some traps at this location
      std::vector<Trap> & traps = it->second;

      AutoDScalar H=0;
      // loop over all traps at this location
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_n,sigma_p,g_n,g_p,E_t;
        sigma_n = TrapSpecs[ptrap->spec].sigma_n;
        sigma_p = TrapSpecs[ptrap->spec].sigma_p;
        g_n     = TrapSpecs[ptrap->spec].g_n;
        g_p     = TrapSpecs[ptrap->spec].g_p;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        AutoDScalar Cn = sigma_n * theta_n * n * (ptrap->N_tt - ptrap->n_t);   // capture rate
        AutoDScalar En = sigma_n * theta_n / g_n * ni * ptrap->n_t * exp(E_t/kb/Tl);  // emission rate

        AutoDScalar Cp = sigma_p * theta_p * p * ptrap->n_t;                                     // capture rate
        AutoDScalar Ep = sigma_p * theta_p / g_p * ni * (ptrap->N_tt - ptrap->n_t) * exp(-E_t/kb/Tl);   // emission rate

        H += (Cn-En) * (1.5*kb*Tn + EcEi);
        H += (Cp-Ep) * (1.5*kb*Tp + EiEv);
      }
      return H;
    }
    else
    {
      return 0.;
    }
  }
  // }}}

  // {{{ AutoDScalar ElectronTrapHeat(const bool flag_bulk, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tn, const AutoDScalar &Tl, const AutoDScalar &EcEi)
  AutoDScalar ElectronTrapHeat(const bool flag_bulk, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tn, const AutoDScalar &Tl, const AutoDScalar &EcEi)
  {
    AutoDScalar theta_n = 1.0e7*cm/s * sqrt(Tl/300/K);     // electron thermal velocity

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      // yes we have some traps at this location
      std::vector<Trap> & traps = it->second;

      AutoDScalar H=0;
      // loop over all traps at this location
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_n,g_n,E_t;
        sigma_n = TrapSpecs[ptrap->spec].sigma_n;
        g_n     = TrapSpecs[ptrap->spec].g_n;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        AutoDScalar Cn = sigma_n * theta_n * n * (ptrap->N_tt - ptrap->n_t);   // capture rate
        AutoDScalar En = sigma_n * theta_n / g_n * ni * ptrap->n_t * exp(E_t/kb/Tl);  // emission rate

        H += (Cn-En) * (1.5*kb*Tn + EcEi);
      }
      return H;
    }
    else
    {
      return 0.;
    }
  }
  // }}}

  // {{{ void init_node()
  /**
   * Actually create traps as specified in PMI commands
   * Called when PMIs are initialized.
   */
  void init_node()
  {
    for (unsigned int i=0; i<TrapSpecs.size(); i++)
    {
      if (TrapSpecs[i].type!=Bulk) continue;  // we only process bulk traps here

      PetscScalar conc = ReadUserScalarValue(TrapSpecs[i].profile_name); // read concentration from profile
      conc=conc*TrapSpecs[i].prefactor;     // concentration is scaled by the prefactor
      if (conc>0)
        AddTrap(*(*pp_point),i,conc);
    }
  }
  // }}}

  // {{{ void init_bc_node(const std::string & bc_label)
  /**
   * Actually create interface traps as specified in PMI commands
   * Called when PMIs are initialized.
   */
  void init_bc_node(const std::string & bc_label)
  {
    for (unsigned int i=0; i<TrapSpecs.size(); i++)
    {
      // we only process interface traps here
      if (TrapSpecs[i].type!=Interface) continue;

      // does the interface ID match with that specified in PMI command?
      if (TrapSpecs[i].interface_name != bc_label) continue;

      PetscScalar conc = TrapSpecs[i].interface_density;
      if (conc>0)
        AddTrap(*(*pp_point),i,conc);
    }
  }
  // }}}

  // {{{void Calculate(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl)
  /**
   * Calculate trap occupancy. The time-derivative term is included with BDF1 discretization
   */
  void Calculate(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl)
  {

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it!=TrapStore.end())
    {
      PetscScalar theta_n = 1.0e7*cm/s * sqrt(Tl/300/K); // electron thermal velocity
      PetscScalar theta_p = 1.0e7*cm/s * sqrt(Tl/300/K); // hole thermal velocity
      std::vector<Trap> & traps = it->second;
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_n,sigma_p,g_n,g_p,E_t;
        sigma_n = TrapSpecs[ptrap->spec].sigma_n;
        sigma_p = TrapSpecs[ptrap->spec].sigma_p;
        g_n     = TrapSpecs[ptrap->spec].g_n;
        g_p     = TrapSpecs[ptrap->spec].g_p;
        E_t     = TrapSpecs[ptrap->spec].E_t;


        PetscScalar A,B;
        A = sigma_n * theta_n * (n + ni/g_n*exp( E_t/kb/Tl)) + sigma_p * theta_p * (p + ni/g_p*exp(-E_t/kb/Tl));
        B = ptrap->N_tt * ( sigma_n * theta_n * n + sigma_p * theta_p * ni / g_p * exp(-E_t/kb/Tl) );

        if (ReadTime() > ptrap->clock_last)
        {
          // we should consider the time derivative of trap occupancy
          if (ptrap->clock_last > ptrap->clock_last_last)
          {
            // We have two previous time step available, use BDF2 discretization
            PetscScalar d =  1.0 / (ReadTime() - ptrap->clock_last_last);
            PetscScalar r = (ptrap->clock_last - ptrap->clock_last_last) * d;

            A = A + d * (2-r)/(1-r);
            B = B + d * ( 1/r/(1-r) * ptrap->n_t_last - (1-r)/r * ptrap->n_t_last_last );
          }
          else
          {
            // We have two previous time step available, use BDF1 discretization
            PetscScalar d =  1.0 / (ReadTime() - ptrap->clock_last);
            A = A + d;
            B = B + ptrap->n_t_last * d;

          }
        }

        ptrap->n_t = B/A;
      }
    }
  }
  // }}}

  // {{{ void Calculate(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tl)
  /**
   * partial derivatives of trap occupancy w.r.t. local V,n,p,
   * The time-derivative term is included with BDF1 discretization
   */
  void Calculate(const bool flag_bulk, const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &ni, const AutoDScalar &Tl)
  {

    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it!=TrapStore.end())
    {
      AutoDScalar theta_n = 1.0e7*cm/s * sqrt(Tl/300/K);
      AutoDScalar theta_p = 1.0e7*cm/s * sqrt(Tl/300/K);
      std::vector<Trap> & traps = it->second;
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        PetscScalar sigma_n,sigma_p,g_n,g_p,E_t;
        sigma_n = TrapSpecs[ptrap->spec].sigma_n;
        sigma_p = TrapSpecs[ptrap->spec].sigma_p;
        g_n     = TrapSpecs[ptrap->spec].g_n;
        g_p     = TrapSpecs[ptrap->spec].g_p;
        E_t     = TrapSpecs[ptrap->spec].E_t;

        AutoDScalar A,B;
        A = sigma_n * theta_n * (n + ni/g_n*exp( E_t/kb/Tl)) + sigma_p * theta_p * (p + ni/g_p*exp(-E_t/kb/Tl));
        B = ptrap->N_tt * ( sigma_n * theta_n * n + sigma_p * theta_p * ni / g_p * exp(-E_t/kb/Tl) );

        if (ReadTime() > ptrap->clock_last)
        {
          // we should consider the time derivative of trap occupancy
          if (ptrap->clock_last > ptrap->clock_last_last)
          {
            // We have two previous time step available, use BDF2 discretization
            PetscScalar d =  1.0 / (ReadTime() - ptrap->clock_last_last);
            PetscScalar r = (ptrap->clock_last - ptrap->clock_last_last) * d;

            A = A + d * (2-r)/(1-r);
            B = B + d * ( 1/r/(1-r) * ptrap->n_t_last - (1-r)/r * ptrap->n_t_last_last );
          }
          else
          {
            // We have two previous time step available, use BDF1 discretization
            PetscScalar d =  1.0 / (ReadTime() - ptrap->clock_last);
            A = A + d;
            B = B + ptrap->n_t_last * d;
          }
        }

        ptrap->n_t_AD = B/A;
      }
    }
  }
  // }}}

  // {{{ void Update(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl)
  /**
   * Converged results are obtained, so we should update the solution
   */
  void Update(const bool flag_bulk, const PetscScalar &p, const PetscScalar &n, const PetscScalar &ni, const PetscScalar &Tl)
  {
    TrapLocation tloc = TrapLocation(*(*pp_point),flag_bulk?Bulk:Interface);

    TrapStore_t::iterator it = TrapStore.find(tloc);

    if (it != TrapStore.end())
    {
      // Trap::n_t should already contain the correct solution, however let's play safe and calculate again
      Calculate(flag_bulk, p, n, ni, Tl);

      std::vector<Trap> & traps = it->second;
      for (std::vector<Trap>::iterator ptrap=traps.begin(); ptrap!=traps.end(); ptrap++)
      {
        // update trap occupancy
        ptrap->n_t_last_last = ptrap->n_t_last;
        ptrap->n_t_last = ptrap->n_t;

        // update timestamp
        ptrap->clock_last_last = ptrap->clock_last;
        ptrap->clock_last = ReadTime();
      }
    }
  }
  // }}}

  // {{{ int calibrate(std::map<std::string, double> & pmi_calibrate_numeric, std::map<std::string, std::string> & pmi_calibrate_string)
  /**
   * Read the parameters of PMI command, and create TrapSpec accordingly
   */
  int calibrate(const std::vector<Parser::Parameter> & pmi_parameter)
  {

    std::string profile_name, interface_name;

    // default values for parameters
    PetscScalar prefactor = 1.0;
    TrapType type = Bulk;
    TrapChargeType charge_type = Neutral;
    PetscScalar Et = 0.0;
    PetscScalar sigman = 4e-16 * cm * cm;
    PetscScalar sigmap = 4e-16 * cm * cm;
    PetscScalar gn=1.0, gp=1.0;
    PetscScalar Dit=0;
    bool has_profile_name=false, has_interface_name=false,
         has_charge_type=false, has_interface_density=false;
    bool has_some_param=false;

    std::vector<Parser::Parameter>::const_iterator it;
    // string parameters first
    for (it=pmi_parameter.begin(); it!=pmi_parameter.end(); it++)
    {
      if(it->type()!=Parser::STRING) continue;

      std::string name,val;
      name = it->name();
      val  = it->get_string();
      std::transform (name.begin(), name.end(), name.begin(), ToLower());
      std::transform (val.begin(), val.end(), val.begin(), ToLower());
      if (name == "profile")
      {
        has_profile_name = true; has_some_param=true;
        profile_name = "Profile_" + it->get_string();
      }
      if (name == "interface")
      {
        has_interface_name = true; has_some_param=true;
        interface_name = it->get_string();
      }
      if (name == "chargetype")
      {
        has_charge_type = true; has_some_param=true;
        if (val=="donor")
          charge_type=Donor;
        if (val=="acceptor")
          charge_type=Acceptor;
      }
      if (name == "type")
      { 
        has_some_param=true;
        if (val=="interface")
          type = Interface;
      }
    }

    if (!has_some_param)
    {
      // no meaningful parameter is specified
      // maybe the user just want to print
      return 0;
    }

    // for interface trap, we must have interface name
    // for bulk trap, we must have profile name
    genius_assert( (type == Interface && has_interface_name) || (type == Bulk && has_profile_name) );
    genius_assert(has_charge_type);

    // then numerical parameters
    for (it=pmi_parameter.begin(); it!=pmi_parameter.end(); it++)
    {
      if(it->type()!=Parser::REAL) continue;

      std::string name;
      PetscScalar val;
      name = it->name();
      val  = it->get_real();
      std::transform (name.begin(), name.end(), name.begin(), ToLower());
      if (name == "prefactor")
        prefactor = val;
      if (name == "energy")
        Et=val * eV;
      if (name == "sigman")
        sigman= val * cm*cm;
      if (name == "sigmap")
        sigmap= val * cm*cm;
      if (name == "gn")
        gn = val;
      if (name == "gp")
        gp = val;
      if (name == "if.density")
      {
        Dit= val / (cm*cm);
        has_interface_density = true;
      }
    }

    if (type==Interface)
    {
      // for interface traps, we must have interface density
      genius_assert( has_interface_density );

      // append a new class of interface traps,
      // actual traps are created later on
      AddInterfaceTrapSpec(interface_name,Dit,charge_type,Et,sigman,sigmap,gn,gp);
    }
    else
    {
      // append a new class of bulk traps,
      // actual traps are created later on
      AddBulkTrapSpec(profile_name,prefactor,charge_type,Et,sigman,sigmap,gn,gp);
    }
    TrapStore.clear();
    return 0;
  }
  // }}}

  // {{{ std::string get_parameter_string(const int verbosity)
  std::string get_parameter_string(const int verbosity)
  {
    std::stringstream output;
   
    // print out bulk traps first
    bool header_printed = false; 
    for ( std::vector<TrapSpec>::const_iterator it = TrapSpecs.begin();
          it != TrapSpecs.end(); it++ )
    {
      if (it->type != Bulk)
        continue;
      
      int wd_name=15, wd_wide=20, wd_narrow=12;
      if ( !header_printed) 
      {
        output << "Bulk traps:" << std::endl;
        output << std::setw(wd_name)   << "Profile Name"
               << std::setw(wd_narrow) << "Prefactor"
               << std::setw(wd_narrow) << "Et (eV)"
               << std::setw(wd_wide)   << "sigma_n (cm^-2)"
               << std::setw(wd_wide)   << "sigma_p (cm^-2)"
               << std::setw(wd_narrow) << "g_n"
               << std::setw(wd_narrow) << "g_p"
               << std::endl;
        header_printed = true;
      }

      output   << std::setw(wd_name)   << it->profile_name
               << std::setw(wd_narrow) << it->prefactor
               << std::setw(wd_narrow) << it->E_t
               << std::setw(wd_wide)   << it->sigma_n
               << std::setw(wd_wide)   << it->sigma_p
               << std::setw(wd_narrow) << it->g_n
               << std::setw(wd_narrow) << it->g_p
               << std::endl;
    }
    if (header_printed)
      output   << std::endl;

    header_printed = false;
    for ( std::vector<TrapSpec>::const_iterator it = TrapSpecs.begin();
          it != TrapSpecs.end(); it++ )
    {
      if (it->type != Interface)
        continue;
      
      int wd_name=15, wd_wide=20, wd_narrow=12;
      if ( !header_printed) 
      {
        output << "Interface traps:" << std::endl;
        output << std::setw(wd_name)   << "Interface Name"
               << std::setw(wd_wide)   << "conc. (cm^-2)"
               << std::setw(wd_narrow) << "Et (eV)"
               << std::setw(wd_wide)   << "sigma_n (cm^-2)"
               << std::setw(wd_wide)   << "sigma_p (cm^-2)"
               << std::setw(wd_narrow) << "g_n"
               << std::setw(wd_narrow) << "g_p"
               << std::endl;
        header_printed = true;
      }

      output   << std::setw(wd_name)   << it->interface_name
               << std::setw(wd_narrow) << it->interface_density
               << std::setw(wd_narrow) << it->E_t
               << std::setw(wd_wide)   << it->sigma_n
               << std::setw(wd_wide)   << it->sigma_p
               << std::setw(wd_narrow) << it->g_n
               << std::setw(wd_narrow) << it->g_p
               << std::endl;
    }
    if (header_printed)
      output   << std::endl;

    return output.str();
  }
  // }}}

  ~GSS_Si_Trap_Default()
  {
  }

}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Trap* PMIS_Si_Trap_Default (const PMIS_Environment& env)
  {
    return new GSS_Si_Trap_Default(env);
  }
}

