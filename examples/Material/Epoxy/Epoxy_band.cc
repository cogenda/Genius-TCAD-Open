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
// Material Type: Epoxy


#include "PMI.h"

class GSS_Epoxy_BandStructure : public PMII_BandStructure
{
private:

  PetscScalar FN_A;      // Coefficient of the pre-exponential term for the Fowler-Nordheim  tunneling model
  PetscScalar FN_B;      // Coefficient of the exponential term for the Fowler-Nordheim tunneling model.

  PetscScalar JE_ALPHA;
  PetscScalar JE_BETA;

  void   Band_Init()
  {
    FN_A      = 6.320000E-07*A/(V*V);
    FN_B      = 2.210000E+08*V/cm;

    JE_ALPHA  = 1e-17*A/(m*m);
    JE_BETA   = 1.0;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("JE_ALPHA",    PARA("JE_ALPHA", "The parameter of emit current density", "A/m^2", A/(m*m), &JE_ALPHA)) );
    parameter_map.insert(para_item("JE_BETA",     PARA("JE_BETA",  "The power parameter of emit current relative to electrical field", "-", 1.0, &JE_BETA)) );
#endif

  }

public:

  PetscScalar Eg            (const PetscScalar &Tl) const
  { return 10.0*eV;  }

  PetscScalar EffecElecMass (const PetscScalar &Tl) const { return me; }

  PetscScalar EffecHoleMass  (const PetscScalar &Tl) const { return me; }

  PetscScalar ARichardson() const { return 1.100000e+02*A/(K*cm)/(K*cm); }

  PetscScalar ElecInject(const PetscScalar &W, const PetscScalar &E, const PetscScalar &Tl) const
  {
    PetscScalar E0 = 1*V/m;
    return JE_ALPHA*std::pow(E/E0, JE_BETA);
  }

  AutoDScalar ElecInject(const AutoDScalar &W, const AutoDScalar &E, const AutoDScalar &Tl) const
  {
    PetscScalar E0 = 1*V/m;
    return JE_ALPHA*adtl::pow(E/E0, JE_BETA);
  }

  /*
  PetscScalar ElecInject(const PetscScalar &W, const PetscScalar &E, const PetscScalar &Tl)
  {
    PetscScalar FN_A = 6.320000E-07*A/(V*V);  // Coefficient of the pre-exponential term for the Fowler-Nordheim  tunneling model
    PetscScalar FN_B = 2.210000E+08*V/cm;     // Coefficient of the exponential term for the Fowler-Nordheim tunneling model.
    return FN_A*E*E*exp(-FN_B/E);
  }

  AutoDScalar ElecInject(const AutoDScalar &W, const AutoDScalar &E, const AutoDScalar &Tl)
  {
    PetscScalar FN_A = 6.320000E-07*A/(V*V);  // Coefficient of the pre-exponential term for the Fowler-Nordheim  tunneling model
    PetscScalar FN_B = 2.210000E+08*V/cm;     // Coefficient of the exponential term for the Fowler-Nordheim tunneling model.
    return FN_A*E*E*exp(-FN_B/E);
  }
  */

  PetscScalar HCI_Barrier_n(const PetscScalar &affinity_semi, const PetscScalar &, const PetscScalar &affinity_ins,
                            const PetscScalar &, const PetscScalar &) const
  { return affinity_semi - 0.0; }

  PetscScalar HCI_Barrier_p(const PetscScalar &affinity_semi, const PetscScalar & Eg_semi, const PetscScalar &affinity_ins,
                            const PetscScalar &, const PetscScalar &) const
  { return 0.0 - affinity_semi -  Eg_semi; }

  PetscScalar HCI_Probability_Insulator_n(const PetscScalar &, const PetscScalar &) const
  { return 0.0;  }

  PetscScalar HCI_Probability_Insulator_p(const PetscScalar &, const PetscScalar &) const
  { return 0.0;  }

  PetscScalar J_FN_Tunneling(const PetscScalar &E_ins, const PetscScalar &alpha) const
  {
    PetscScalar E = fabs(E_ins) + 1*V/cm;
    if(alpha == 1.0)
    {
      if( FN_B/E > 70.0 ) return 0.0;
      return FN_A*E_ins*E_ins*exp( - FN_B/E );
    }
    else
    {

      if( FN_B/E*(1-std::pow(1.0-alpha, 1.5)) > 70.0 ) return 0.0;
      return FN_A*E_ins*E_ins/std::pow(1.0-sqrt(1.0-alpha), 2.0)*exp( - FN_B/E*(1-std::pow(1.0-alpha, 1.5)) );
    }
  }

  PetscScalar J_CBET_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                               const PetscScalar &Efn1, const PetscScalar &Efn2,
                               const PetscScalar &Ec1,  const PetscScalar &Ec2,
                               const PetscScalar &B1,   const PetscScalar &B2,
                               const PetscScalar &t) const
  { return 0.0; }

  AutoDScalar J_CBET_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                               const AutoDScalar &Efn1, const AutoDScalar &Efn2,
                               const AutoDScalar &Ec1,  const AutoDScalar &Ec2,
                               const AutoDScalar &B1,   const AutoDScalar &B2,
                               const PetscScalar &t) const
  { return 0.0; }


  PetscScalar J_VBHT_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                               const PetscScalar &Efp1, const PetscScalar &Efp2,
                               const PetscScalar &Ev1,  const PetscScalar &Ev2,
                               const PetscScalar &B1,   const PetscScalar &B2,
                               const PetscScalar &t) const
  { return 0.0; }

  AutoDScalar J_VBHT_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                               const AutoDScalar &Efp1, const AutoDScalar &Efp2,
                               const AutoDScalar &Ev1,  const AutoDScalar &Ev2,
                               const AutoDScalar &B1,   const AutoDScalar &B2,
                               const PetscScalar &t) const
  { return 0.0; }

  PetscScalar J_VBET_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                               const PetscScalar &Efn1, const PetscScalar &Efn2,
                               const PetscScalar &Ec1,  const PetscScalar &Ec2,
                               const PetscScalar &Ev1,  const PetscScalar &Ev2,
                               const PetscScalar &B1,   const PetscScalar &B2,
                               const PetscScalar &t) const
  { return 0.0; }

  AutoDScalar J_VBET_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                               const AutoDScalar &Efn1, const AutoDScalar &Efn2,
                               const AutoDScalar &Ec1,  const AutoDScalar &Ec2,
                               const AutoDScalar &Ev1,  const AutoDScalar &Ev2,
                               const AutoDScalar &B1,   const AutoDScalar &B2,
                               const PetscScalar &t) const
  { return 0.0; }



  GSS_Epoxy_BandStructure(const PMII_Environment &env):PMII_BandStructure(env)
  {
    Band_Init();
  }

  ~GSS_Epoxy_BandStructure()
  {}
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_BandStructure* PMII_Epoxy_BandStructure_Default (const PMII_Environment& env)
  {
    return new GSS_Epoxy_BandStructure(env);
  }
}
