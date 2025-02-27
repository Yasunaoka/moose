//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FsiApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

#include "NavierStokesApp.h"
#include "TensorMechanicsApp.h"

InputParameters
FsiApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  return params;
}

FsiApp::FsiApp(InputParameters parameters) : MooseApp(parameters)
{
  FsiApp::registerAll(_factory, _action_factory, _syntax);
}

FsiApp::~FsiApp() {}

void
FsiApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"FsiApp"});
  Registry::registerActionsTo(af, {"FsiApp"});

  NavierStokesApp::registerAll(f, af, s);
  TensorMechanicsApp::registerAll(f, af, s);
}

void
FsiApp::registerApps()
{
  registerApp(FsiApp);

  NavierStokesApp::registerApps();
  TensorMechanicsApp::registerApps();
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
FsiApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  FsiApp::registerAll(f, af, s);
}
extern "C" void
FsiApp__registerApps()
{
  FsiApp::registerApps();
}
