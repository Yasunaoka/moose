//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Reference book: 《计算动力学》，张雄，王天舒，清华大学出版社
#include "ComputeStressBase.h"
#include "Material.h"
#include "PropertyReadFile.h"


/**
 * JCEOSnl implements rate-dependent associative J2 plasticity
 * with isotropic hardening in the finite-strain framework.
 * Yield function = sqrt(3*s_ij*s_ij/2) - K(JC)
 * where s_ij = stress_ij - delta_ij*trace(stress)/3 is the deviatoric stress
 * and K is the yield stress.
 * the hydrostatic stress part updates by Gruneisen EOS.
 * no iteration oweing to small step in explicite dynamics
 */
class JCEOSnl : public ComputeStressBase
{
public:
  static InputParameters validParams();

  JCEOSnl(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpStress();

  //JCpart
  Real _A;
  Real _B;
  Real _n;
  Real _C;
  Real _reference_strain_rate;
  Real _m;
  Real _T_m;
  Real _T_r;
  Real _T_i;
  Real _SpecHeat;


  MaterialProperty<Real>&_yield_factor;
  MaterialProperty<Real> & _T;
  const MaterialProperty<Real>& _T_old;
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  const MaterialProperty<RankTwoTensor> & _strain_rate;
  MaterialProperty<Real> & _eqv_plastic_strain;
  const MaterialProperty<Real> & _eqv_plastic_strain_old;

  //eos part
  Real _D_i;
  Real _ShearModule;
  Real _C0;
  Real _S1;
  Real _gamma0;
  Real _a;

  const MaterialProperty<Real>& _density;
  const Real & _current_elem_volume;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  MaterialProperty<Real> & _element_vol;
  const MaterialProperty<Real> & _element_vol_old;
  MaterialProperty<Real> & _E;
  const MaterialProperty<Real> & _E_old;

  //phase field dgration
  const VariableValue & _c;
  MaterialProperty<Real> & _shearband_driving_force;
  const MaterialProperty<Real> & _shearband_driving_force_old;
  const MaterialProperty<Real> & _D;


  const PropertyReadFile * const _read_prop_user_object;



  Real getYieldStress(const Real equivalent_plastic_strain, 
                      const Real A, 
                      const Real B, 
                      const Real n, 
                      const Real C, 
                      Real effstrain_rate, 
                      const Real reference_strain_rate,
                      const Real T,
                      const Real m,
                      const Real T_m,
                      const Real T_r,
                      const Real yield_factor);


  Real getdYieldStressdPlasticStrain(const Real equivalent_plastic_strain, 
                                      const Real B, 
                                      const Real n, 
                                      const Real C, 
                                      Real effstrain_rate, 
                                      const Real reference_strain_rate,
                                      const Real T,
                                      const Real m,
                                      const Real T_m,
                                      const Real T_r,
                                      const Real yield_factor);


  usingTensorIndices(i_, j_, k_, l_);
};
