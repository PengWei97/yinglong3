#pragma once
#include "FeatureDataVectorPostprocessor.h" 

// FeatureMatPropCplVectorPostprocessor

enum class state_variable {slip_resistance, backstress};

class FeatureMatPropCplVectorPostprocessor : public FeatureDataVectorPostprocessor
{
public:
  static InputParameters validParams();

  virtual void execute() override;
  virtual void finalize() override;

  FeatureMatPropCplVectorPostprocessor(const InputParameters & parameters);

  // _slip_resistance, _backstress
  std::vector<Real> getStateVariable(unsigned int feature_id, const state_variable & state_variable_name  = state_variable::slip_resistance);

protected:
  /// Add volume contributions to one or entries in the feature volume vector
  void accumulateStateVariables(const Elem * elem,
                                 const std::vector<unsigned int> & var_to_features,
                                 std::size_t num_features);

  std::vector<Real> computeStateVarivalesIntegral(std::size_t var_index, const MaterialProperty<std::vector<Real>> & state_variable_copy) const;

  void sum_state_variables(std::vector<std::vector<Real>> & stat_variables, const unsigned int & num_features);

  const std::string _base_name;
  // Get from material class CPKalindindiCplUpdate
  const MaterialProperty<std::vector<Real>> & _slip_resistance_copy;
  const MaterialProperty<std::vector<Real>> & _backstress_copy;

  const unsigned int _number_slip_systems;

  std::vector<std::vector<Real>> & _slip_resistances;
  std::vector<std::vector<Real>> & _backstresses;  
};