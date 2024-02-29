//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "JCEOSnl.h"
#include "libmesh/utility.h"
#include "MathUtils.h"
#include "Assembly.h"

registerMooseObject("TensorMechanicsApp", JCEOSnl);

InputParameters
JCEOSnl::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Associative J2 plasticity with jc and gruneisen.");

  //JCpart
  params.addRequiredParam<Real>("A", "A+B(epsilon_p)^n");
  params.addRequiredParam<Real>("B", "A+B(epsilon_p)^n");
  params.addRequiredParam<Real>("n", "A+B(epsilon_p)^n");
  params.addRequiredParam<Real>("C", "1+Cln(epsilon_p_dot/reference epsilon_dot)");
  params.addRequiredParam<Real>("reference_strain_rate", "1+Cln(epsilon_p_dot/reference epsilon_dot)");
  params.addRequiredParam<Real>("m", "1-((T-T_r)/(T_m-T_r))^m");
  params.addRequiredParam<Real>("T_m", "1-((T-T_r)/(T_m-T_r))^m,melting_T");
  params.addRequiredParam<Real>("T_r", "1-((T-T_r)/(T_m-T_r))^m,reference_T");
  params.addRequiredParam<Real>("T_i", "initial temp");
  params.addRequiredParam<Real>("SpecHeat", "Name of Material Property that provides the SpecHeat");

  //eos part
  params.addRequiredParam<Real>("D_i", "initial density");
  params.addRequiredParam<Real>("ShearModule", "for the deviatoric stress update");
  params.addRequiredParam<Real>("C0", "us = C0 + S1 * up");
  params.addRequiredParam<Real>("S1", "us = C0 + S1 * up");
  params.addRequiredParam<Real>("gamma0", "Gruneusen parameter");
  params.addRequiredParam<Real>("a", "gamma = a + (gamma0 - a) * V/V0");

  //phase field dgration
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addParam<MaterialPropertyName>(
      "shearband_driving_force_name",
      "local_energy_dev",
      "Name of material property for local mechanical energy function.");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");


  params.addParam<UserObjectName>("read_prop_user_object",
                                  "The ElementReadPropertyFile "
                                  "GeneralUserObject to read element "
                                  "specific property values from file");

  return params;
}

JCEOSnl::JCEOSnl(const InputParameters & parameters)
  : ComputeStressBase(parameters),

  //JCpart
  _A(getParam<Real>("A")),
  _B(getParam<Real>("B")),
  _n(getParam<Real>("n")),
  _C(getParam<Real>("C")),
  _reference_strain_rate(getParam<Real>("reference_strain_rate")),
  _m(getParam<Real>("m")),
  _T_m(getParam<Real>("T_m")),
  _T_r(getParam<Real>("T_r")),
  _T_i(getParam<Real>("T_i")),
  _SpecHeat(getParam<Real>("SpecHeat")),


  _yield_factor(declareProperty<Real>(_base_name +"yield_factor")),
  _T(declareProperty<Real>(_base_name + "T")),//stateful
  _T_old(getMaterialPropertyOld<Real>(_base_name + "T")),
  _plastic_strain(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),//stateful
  _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain")),
  _strain_rate(getMaterialProperty<RankTwoTensor>(_base_name + "strain_rate")),
  _eqv_plastic_strain(declareProperty<Real>(_base_name + "eqv_plastic_strain")),//stateful
  _eqv_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "eqv_plastic_strain")),



  //eos part
  _D_i(getParam<Real>("D_i")),
  _ShearModule(getParam<Real>("ShearModule")),
  _C0(getParam<Real>("C0")),
  _S1(getParam<Real>("S1")),
  _gamma0(getParam<Real>("gamma0")),
  _a(getParam<Real>("a")),

  _density(getGenericMaterialProperty<Real,false>(_base_name + "density")),
  _current_elem_volume(_assembly.elemVolume()),
  _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),//stateful
  _strain_increment(getMaterialProperty<RankTwoTensor>(_base_name + "strain_increment")),
  _rotation_increment(getMaterialProperty<RankTwoTensor>(_base_name + "rotation_increment")),
  _elasticity_tensor_name(_base_name + "elasticity_tensor"),
  _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
  _element_vol(declareProperty<Real>(_base_name + "element_volume")),
  _element_vol_old(getMaterialPropertyOld<Real>(_base_name + "element_volume")),
  _E(declareProperty<Real>(_base_name + "E")),
  _E_old(getMaterialPropertyOld<Real>(_base_name + "E")),

  //phase field dgration
  _c(coupledValue("c")),
  _shearband_driving_force(declareProperty<Real>(getParam<MaterialPropertyName>("shearband_driving_force_name"))),
  _shearband_driving_force_old(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("shearband_driving_force_name"))),
  _D(getMaterialProperty<Real>("D_name")),    


  _read_prop_user_object(isParamValid("read_prop_user_object")
                            ? &getUserObject<PropertyReadFile>("read_prop_user_object")
                            : nullptr)
{

}

void
JCEOSnl::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _plastic_strain[_qp].zero();
  _eqv_plastic_strain[_qp] = 0.0;
  _T[_qp] = _T_i;
  _E[_qp] = 0.0;
  _shearband_driving_force[_qp] = 0.0;
}

void
JCEOSnl::computeQpStress()
{
  _element_vol[_qp] = _current_elem_volume;
    if (_read_prop_user_object)
  {
    _yield_factor[_qp] = _read_prop_user_object->getData(_current_elem, 0);
  }
  else
  {
    _yield_factor[_qp] = 1;
  }
  // calculate the deviatoric part
  RankTwoTensor stress_dev_old = _stress_old[_qp].deviatoric();
  RankFourTensor dev_E_ijkl = _elasticity_tensor[_qp];
  dev_E_ijkl(0, 0, 0, 0) = 2 * _ShearModule;
  dev_E_ijkl(1, 1, 1, 1) = 2 * _ShearModule;
  dev_E_ijkl(2, 2, 2, 2) = 2 * _ShearModule;
  dev_E_ijkl(0, 0, 1, 1) = 0.0;
  dev_E_ijkl(1, 1, 0, 0) = 0.0;
  dev_E_ijkl(0, 0, 2, 2) = 0.0;
  dev_E_ijkl(2, 2, 0, 0) = 0.0;
  dev_E_ijkl(1, 1, 2, 2) = 0.0;
  dev_E_ijkl(2, 2, 1, 1) = 0.0;
  // RankTwoTensor undegration_dstress_dev = dev_E_ijkl * _strain_increment[_qp].deviatoric();
  RankTwoTensor dstress_dev = _D[_qp] * dev_E_ijkl * _strain_increment[_qp].deviatoric();
  RankTwoTensor stress_dev = stress_dev_old + dstress_dev;//trial stress
  RankTwoTensor trial_stress_dev = stress_dev;
  Real effstrainrate = std::sqrt(2.0 / 3.0 * _strain_rate[_qp].doubleContraction(_strain_rate[_qp]));
  Real yield_stress = getYieldStress(_eqv_plastic_strain_old[_qp], 
                                      _A, 
                                      _B, 
                                      _n,  
                                      _C, 
                                      effstrainrate, 
                                      _reference_strain_rate,
                                      _T_old[_qp],
                                      _m,
                                      _T_m,
                                      _T_r,
                                      _yield_factor[_qp]); // yield stress at old step
  Real Ep = getdYieldStressdPlasticStrain(_eqv_plastic_strain_old[_qp], 
                                          _B, 
                                          _n,  
                                          _C, 
                                          effstrainrate, 
                                          _reference_strain_rate,
                                          _T_old[_qp],
                                          _m,
                                          _T_m,
                                          _T_r,
                                          _yield_factor[_qp]);

  Real trial_eqv_stress = std::sqrt(1.5 * trial_stress_dev.doubleContraction(trial_stress_dev));
  if (yield_stress < 0)
  {
    mooseError("yield_stress < 0 !",_current_elem->id());
  }
  // plastic loading occurs (maybe)
  if (trial_eqv_stress > _D[_qp] * yield_stress)
  {
    mooseInfo("yield?");
    Real deqv_plastic_strain = (trial_eqv_stress - _D[_qp] * yield_stress)/(3 * _ShearModule + Ep);
    Real dshearband_driving_force = yield_stress * deqv_plastic_strain;
    _eqv_plastic_strain[_qp] = _eqv_plastic_strain_old[_qp] + deqv_plastic_strain;
    // yield_stress += Ep * deqv_plastic_strain;//yield stress at this step
    yield_stress = getYieldStress(_eqv_plastic_strain[_qp], 
                                  _A, 
                                  _B, 
                                  _n,  
                                  _C, 
                                  effstrainrate, 
                                  _reference_strain_rate,
                                  _T_old[_qp],
                                  _m,
                                  _T_m,
                                  _T_r,
                                  _yield_factor[_qp]); // yield stress at this step
    if (trial_eqv_stress < _D[_qp] * yield_stress)//indicates no plastic in fact
    {
      mooseInfo("not yet...");
      _eqv_plastic_strain[_qp] = _eqv_plastic_strain_old[_qp];
      dshearband_driving_force = 0;
    }
    else//yield occuers
    {
      mooseInfo("indeed!");
      Real ratio = _D[_qp] * yield_stress/trial_eqv_stress;
      RankTwoTensor dplastic_strain = dev_E_ijkl.invSymm() * ((1-ratio) * stress_dev);
      _plastic_strain[_qp] = _plastic_strain_old[_qp] + dplastic_strain;
      stress_dev *= ratio;
      Real eqv_stress = std::sqrt(1.5 * stress_dev.doubleContraction(stress_dev));
      Real dT = (eqv_stress * deqv_plastic_strain)/(_density[_qp] * _SpecHeat);
      _T[_qp] = _T_old[_qp] + dT;
      _shearband_driving_force[_qp] = _shearband_driving_force_old[_qp] + dshearband_driving_force;//this term won't decrase,so the hist is not needed
    }
  }
  RankTwoTensor stress_vol_old = _stress_old[_qp] - stress_dev_old;
  _element_vol[_qp] = _current_elem_volume;//restoring for next step
  Real dE = _current_elem_volume * 0.5 * (stress_dev + stress_dev_old).doubleContraction(_strain_increment[_qp]);
  dE = dE + 0.5 * (_current_elem_volume - _element_vol_old[_qp]) * stress_vol_old.trace()/3;// waiting for q updating......
  _E[_qp] = _E_old[_qp] + dE;
  Real miu = (_density[_qp]/_D_i) - 1;
  RankTwoTensor stress_vol;
  stress_vol.zero();
  if (miu >= 0)
  {
    stress_vol.addIa(-(_D_i * _C0 * _C0 * miu * (1 + (1 - 0.5 * _gamma0) * miu - 0.5 * _a * miu * miu))/std::pow(1 - (_S1 - 1) * miu, 2));
  }
  else
  {
    stress_vol.addIa(-_D_i * _C0 * _C0 * miu);
  }
  Real V0 = _current_elem_volume * _density[_qp]/_D_i;
  stress_vol.addIa(-(_gamma0 + _a * miu) * _E[_qp]/V0);
  if (miu < 0)//degrad pressure when expand
  {
    stress_vol *= _D[_qp];
  }
  
  _E[_qp] = _E[_qp] + 0.5 * (_current_elem_volume - _element_vol_old[_qp]) * stress_vol.trace()/3;
  _stress[_qp] = stress_vol + stress_dev;
  // Rotate the stress tensor to the current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();

  // Rotate plastic strain tensor to the current configuration
  _plastic_strain[_qp] =_rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();

  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];
  
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];

}

// Obtain yield stress
Real
JCEOSnl::getYieldStress(const Real equivalent_plastic_strain, 
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
                        const Real yield_factor)
  {
  if (effstrain_rate < 0.001)
  {
    effstrain_rate = 0.001;
  }

  Real T_dim = (T - T_r)/(T_m - T_r);

  if (T_dim <= 0)
  {
    T_dim = 0.00001;
  }
  if (T_dim >= 1)
  {
    T_dim = 0.99999;
  }
  
  
  
  return yield_factor * (
                          (A + B * std::pow(equivalent_plastic_strain, n)) *
                          (1 + C * std::log(effstrain_rate/ reference_strain_rate)) *
                          (1 - std::pow(T_dim, m))
                          );
}

Real
JCEOSnl::getdYieldStressdPlasticStrain(const Real equivalent_plastic_strain, 
                                      const Real B, 
                                      const Real n, 
                                      const Real C, 
                                      Real effstrain_rate, 
                                      const Real reference_strain_rate,
                                      const Real T,
                                      const Real m,
                                      const Real T_m,
                                      const Real T_r,
                                      const Real yield_factor)
{
    if (effstrain_rate < 0.001)
  {
    effstrain_rate = 0.001;//to aviod yield_stress < 0
  }

  Real T_dim = (T - T_r)/(T_m - T_r);

  if (T <= T_r)
  {
    T_dim = 0.00001;
  }
  if (T >= T_r)
  {
    T_dim = 0.99999;
  }


  return yield_factor * (
                          (B * n * std::pow(equivalent_plastic_strain+0.0000000001, n-1)) *
                          (1 + C * std::log(effstrain_rate / reference_strain_rate)) *
                          (1 - std::pow(T_dim, m))
                          );//to aviod equivalent_plastic_strain=0 which leads equivalent_plastic_strain^-0.64 = infinite
}