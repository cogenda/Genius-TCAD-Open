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
// Material Type: Air


#include "PMI.h"

class GSS_Air_BandStructure : public PMII_BandStructure
{
private:

  void   Band_Init()
  {}

public:

  PetscScalar Eg            (const PetscScalar &Tl) const
  { return 0.0;  }

  PetscScalar HCI_Barrier_n(const PetscScalar &affinity_semi, const PetscScalar &,
                            const PetscScalar &, const PetscScalar &) const
  { return affinity_semi - 0.0; }

  PetscScalar HCI_Barrier_p(const PetscScalar &affinity_semi, const PetscScalar & Eg_semi,
                            const PetscScalar &, const PetscScalar &) const
  { return 0.0 - affinity_semi -  Eg_semi; }

  PetscScalar HCI_Probability_Insulator_n(const PetscScalar &, const PetscScalar &) const
  { return 0.0;  }

  PetscScalar HCI_Probability_Insulator_p(const PetscScalar &, const PetscScalar &) const
  { return 0.0;  }

  PetscScalar J_FN_Tunneling(const PetscScalar &E_ins, const PetscScalar &alpha) const
  { return 0.0; }

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



  GSS_Air_BandStructure(const PMII_Environment &env):PMII_BandStructure(env)
  {
    Band_Init();
  }

  ~GSS_Air_BandStructure()
  {}
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_BandStructure* PMII_Air_BandStructure_Default (const PMII_Environment& env)
  {
    return new GSS_Air_BandStructure(env);
  }
}
