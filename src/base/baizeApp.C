#include "baizeApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
baizeApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

baizeApp::baizeApp(InputParameters parameters) : MooseApp(parameters)
{
  baizeApp::registerAll(_factory, _action_factory, _syntax);
}

baizeApp::~baizeApp() {}

void 
baizeApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<baizeApp>(f, af, s);
  Registry::registerObjectsTo(f, {"baizeApp"});
  Registry::registerActionsTo(af, {"baizeApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
baizeApp::registerApps()
{
  registerApp(baizeApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
baizeApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  baizeApp::registerAll(f, af, s);
}
extern "C" void
baizeApp__registerApps()
{
  baizeApp::registerApps();
}
