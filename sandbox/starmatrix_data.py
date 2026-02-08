import starmatrix
import starmatrix.settings as smsettings
from starmatrix.model import Model
import starmatrix.matrix as matrix


def config(sm_params):
  custom_params = smsettings.default_settings()
  custom_params['binary_fraction'] = sm_params['binary_fraction']
  custom_params['sol_ab'] = sm_params['sol_ab']
  custom_params['dtd_sn'] = sm_params['dtd_sn']

  context = smsettings.validate(custom_params)

  return context

def sn_matrix(stellar_mass):
  contribution_matrix = matrix.q_sn(stellar_mass)
  return contribution_matrix
