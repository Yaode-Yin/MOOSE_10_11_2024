//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledGradTimeDerivativeKernel.h"

registerMooseObject("PhaseFieldApp", CoupledSusceptibilityTimeDerivative);

InputParameters
CoupledGradTimeDerivativeKernel::validParams()
{
  //InputParameters params = JvarMapKernelInterface<CoupledTimeDerivative>::validParams();
  InputParameters params = Kernel::validParams();
  params.addClassDescription("A modified coupled time derivative Kernel that multiplies the time "
                             "derivative of the gradient of a coupled variable");
  params.addRequiredCoupledVar("v", "Coupled variable");
  params.addRequiredParam<MaterialPropertyName>(
      "f_name", "Susceptibility function F defined in a FunctionMaterial");
  return params;
}

CoupledGradTimeDerivativeKernel::CoupledGradTimeDerivativeKernel(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _F(getMaterialProperty<Real>("f_name")),
    _dFdu(getMaterialPropertyDerivative<Real>("f_name", _var.name())),
    _grad_v_dot(coupledGradientDot("v")), 
    _grad_v(coupledGradient("v"))
  {
  }

void
CoupledGradTimeDerivativeKernel::initialSetup()
{
  validateNonlinearCoupling<Real>("f_name");
}

Real
CoupledGradTimeDerivativeKernel::computeQpResidual()
{
  return _test[_i][_qp] * _grad_v_dot[_qp] * _grad_v[_qp] * _u[_qp] * _F[_qp];
}

Real
CoupledGradTimeDerivativeKernel::computeQpJacobian()
{
  return _test[_i][_qp] * _grad_v_dot[_qp] * _grad_v[_qp] * _phi[_j][_qp] *_F[_qp]  + _test[_i][_qp] * _grad_v_dot[_qp] * _grad_v[_qp] * _u[_qp] * _dFdu[_qp];
}

// Real
// CoupledSusceptibilityTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
//{
  // get the coupled variable jvar is referring to
 // const unsigned int cvar = mapJvarToCvar(jvar);

 // if (jvar == _v_var)
  //  return CoupledTimeDerivative::computeQpOffDiagJacobian(jvar) * _F[_qp] +
 //          CoupledTimeDerivative::computeQpResidual() * _phi[_j][_qp] * (*_dFdarg[cvar])[_qp];

//  return CoupledTimeDerivative::computeQpResidual() * _phi[_j][_qp] * (*_dFdarg[cvar])[_qp];
// }
