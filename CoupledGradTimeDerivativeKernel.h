//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"
// Forward Declaration

/*** This calculates a modified coupled time derivative that multiplies the time derivative of the grandient of a coupled variable by a function of the variables
 */
class CoupledGradTimeDerivativeKernel
  : public DerivativeMaterialInterface<Kernel>
 /// : public DerivativeMaterialInterface<JvarMapKernelInterface<CoupledTimeDerivative>>
{
public:
  static InputParameters validParams();

  CoupledGradTimeDerivativeKernel(const InputParameters & parameters);
  virtual void initialSetup();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  
  const VariableGradient & _grad_v_dot;
  const VariableGradient & _grad_v;

  /// The function multiplied by the coupled time derivative
  const MaterialProperty<Real> & _F;

  /// function derivative w.r.t. the kernel variable
  const MaterialProperty<Real> & _dFdu;

};
