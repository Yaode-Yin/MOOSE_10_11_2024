//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledGradTimeDerivative.h"

registerMooseObject("MooseApp", CoupledGradTimeDerivative);

InputParameters
CoupledGradTimeDerivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Time derivative Kernel that acts on a coupled variable. Weak form: "
                             "$(\\psi_i, \\frac{\\partial v_h}{\\partial t})$.");
  params.addRequiredCoupledVar("v", "Coupled variable");
  return params;
}

CoupledGradTimeDerivative::CoupledGradTimeDerivative(const InputParameters & parameters)
  : ArrayKernel(parameters), _grad_v_dot(coupledGradientDot("v"))
{
}

Real
CoupledGradTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * _grad_v_dot[_qp];
}

Real
CoupledGradTimeDerivative::computeQpJacobian()
{
  return 0.0;
}

//Real
//CoupledTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
//{
//  if (jvar == _v_var)
//    return _test[_i][_qp] * _phi[_j][_qp] * _dv_dot[_qp];
//
//  return 0.0;
//}
