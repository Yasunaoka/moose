//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NEML2ModelInterface.h"

/**
 * Interface class to provide common input parameters, members, and methods for MOOSEObjects that
 * use NEML2 solid mechanics models.
 */
template <class T>
class NEML2SolidMechanicsInterface : public NEML2ModelInterface<T>
{
public:
  static InputParameters validParams();

  template <typename... P>
  NEML2SolidMechanicsInterface(const InputParameters & params, P &&... args);

#ifdef NEML2_ENABLED
protected:
  using NEML2ModelInterface<T>::getLabeledAxisAccessor;
  using NEML2ModelInterface<T>::model;

  virtual void validateModel() const override;

  // @{ Getters for NEML2 variable accessors
  const neml2::LabeledAxisAccessor & stress() const { return _neml2_stress; }
  const neml2::LabeledAxisAccessor & strain() const { return _neml2_strain; }
  const neml2::LabeledAxisAccessor & temperature() const { return _neml2_temperature; }
  const neml2::LabeledAxisAccessor & time() const { return _neml2_time; }

  const neml2::LabeledAxisAccessor & stressOld() const { return _neml2_stress_n; }
  const neml2::LabeledAxisAccessor & strainOld() const { return _neml2_strain_n; }
  const neml2::LabeledAxisAccessor & temperatureOld() const { return _neml2_temperature_n; }
  const neml2::LabeledAxisAccessor & timeOld() const { return _neml2_time_n; }
  // @}

private:
  // @{ Variable accessors for the NEML2 material model
  const neml2::LabeledAxisAccessor _neml2_stress;
  const neml2::LabeledAxisAccessor _neml2_strain;
  const neml2::LabeledAxisAccessor _neml2_temperature;
  const neml2::LabeledAxisAccessor _neml2_time;

  const neml2::LabeledAxisAccessor _neml2_stress_n;
  const neml2::LabeledAxisAccessor _neml2_strain_n;
  const neml2::LabeledAxisAccessor _neml2_temperature_n;
  const neml2::LabeledAxisAccessor _neml2_time_n;
  // @}
#endif // NEML2_ENABLED
};

template <class T>
InputParameters
NEML2SolidMechanicsInterface<T>::validParams()
{
  InputParameters params = NEML2ModelInterface<T>::validParams();

  params.addParam<std::string>("neml2_stress", "S", "NEML2 variable accessor for stress");
  params.addParam<std::string>("neml2_strain", "E", "NEML2 variable accessor for strain");
  params.addParam<std::string>("neml2_temperature", "T", "NEML2 variable accessor for temperature");
  params.addParam<std::string>("neml2_time", "t", "NEML2 variable accessor for time");

  return params;
}

template <class T>
template <typename... P>
NEML2SolidMechanicsInterface<T>::NEML2SolidMechanicsInterface(const InputParameters & params,
                                                              P &&... args)
  : NEML2ModelInterface<T>(params, args...)
#ifdef NEML2_ENABLED
    ,
    // The LabeledAxisAccessors for the NEML2 variables:
    _neml2_stress(getLabeledAxisAccessor(params.get<std::string>("neml2_stress")).on("state")),
    _neml2_strain(getLabeledAxisAccessor(params.get<std::string>("neml2_strain")).on("forces")),
    _neml2_temperature(
        getLabeledAxisAccessor(params.get<std::string>("neml2_temperature")).on("forces")),
    _neml2_time(getLabeledAxisAccessor(params.get<std::string>("neml2_time")).on("forces")),
    // and their corresponding old state/forces:
    _neml2_stress_n(_neml2_stress.slice(1).on("old_state")),
    _neml2_strain_n(_neml2_strain.slice(1).on("old_forces")),
    _neml2_temperature_n(_neml2_temperature.slice(1).on("old_forces")),
    _neml2_time_n(_neml2_time.slice(1).on("old_forces"))
#endif // NEML2_ENABLED
{
}

#ifdef NEML2_ENABLED
template <class T>
void
NEML2SolidMechanicsInterface<T>::validateModel() const
{
  NEML2ModelInterface<T>::validateModel();

  // The NEML2 model should have strain as an input
  if (!model().input().has_variable(_neml2_strain))
    mooseError("The NEML2 model does not have ", _neml2_strain, " as an input.");

  // The NEML2 model should have stress as an output
  if (!model().output().has_variable(_neml2_stress))
    mooseError("The NEML2 model does not have ", _neml2_stress, " as an output.");
}
#endif
