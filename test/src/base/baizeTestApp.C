//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "baizeTestApp.h"
#include "baizeApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
baizeTestApp::validParams()
{
  InputParameters params = baizeApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

baizeTestApp::baizeTestApp(InputParameters parameters) : MooseApp(parameters)
{
  baizeTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

baizeTestApp::~baizeTestApp() {}

void
baizeTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  baizeApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"baizeTestApp"});
    Registry::registerActionsTo(af, {"baizeTestApp"});
  }
}

void
baizeTestApp::registerApps()
{
  registerApp(baizeApp);
  registerApp(baizeTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
baizeTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  baizeTestApp::registerAll(f, af, s);
}
extern "C" void
baizeTestApp__registerApps()
{
  baizeTestApp::registerApps();
}
