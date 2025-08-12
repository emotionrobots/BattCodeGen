#============================================================================================
#
#  batt_codegen.py
#
#  pyBAMM Battery C code generation capability
#
#============================================================================================
import os 
import pybamm
import casadi 
import subprocess
import numpy as np
from numpy import testing


#-------------------------------------------------------------------------------------------
#  Create DFN model with default configurations 
#-------------------------------------------------------------------------------------------
model = pybamm.lithium_ion.DFN()
geometry = model.default_geometry
params = model.default_parameter_values

params.process_geometry(geometry)
params.process_model(model)

#-------------------------------------------------------------------------------------------
#  Must discretize the model first 
#-------------------------------------------------------------------------------------------
var = pybamm.standard_spatial_vars
var_pts = {var.x_n: 30, var.x_s: 30, var.x_p: 30, var.r_n: 10, var.r_p: 10}
mesh = pybamm.Mesh(geometry, model.default_submesh_types, var_pts)
disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
disc.process_model(model)


#-------------------------------------------------------------------------------------------
#  Generate C code using Casadi
#
#  API: 
#     def generate(self, filename, variable_names, input_parameter_order=None, cg_options=None)
#
#        Parameters
#        ----------
#        filename : str
#            Name of the file to which to save the code
#        variable_names : list
#            Variables to be exported alongside the model structure
#        input_parameter_order : list, optional
#            Order in which the input parameters should be stacked.
#            If input_parameter_order=None and len(self.input_parameters) > 1, a
#            ValueError is raised (this helps to avoid accidentally using the wrong
#            order)
#        cg_options : dict
#            Options to pass to the code generator.
#            See https://web.casadi.org/docs/#generating-c-code
#       
#-------------------------------------------------------------------------------------------
opts = dict(main=True, with_header=True, indent=3, verbose=True) 
model.generate("test.c", ["Battery voltage [V]"], input_parameter_order=None, cg_options=opts)



#-------------------------------------------------------------------------------------------
# Compile
#-------------------------------------------------------------------------------------------
subprocess.run(["gcc", "-fPIC", "-shared", "-o", "test.so", "test.c"])  # nosec

'''
# Read the generated functions
x0_fn = casadi.external("x0", "./test.so")
z0_fn = casadi.external("z0", "./test.so")
rhs_fn = casadi.external("rhs_", "./test.so")
alg_fn = casadi.external("alg_", "./test.so")
jac_rhs_fn = casadi.external("jac_rhs", "./test.so")
jac_alg_fn = casadi.external("jac_alg", "./test.so")
var_fn = casadi.external("variables", "./test.so")

# Test that function values are as expected

#self.assertEqual(x0_fn([2, 5]), 5)
#self.assertEqual(z0_fn([0, 0]), 1)
#self.assertEqual(rhs_fn(0, 3, 2, [7, 2]), -21)
#self.assertEqual(alg_fn(0, 3, 2, [7, 2]), 1)

np.testing.assert_array_equal(np.array(jac_rhs_fn(5, 6, 7, [8, 9])), [[-8, 0]])
np.testing.assert_array_equal(np.array(jac_alg_fn(5, 6, 7, [8, 9])), [[1, -1]])
#self.assertEqual(var_fn(6, 3, 2, [7, 2]), -1)

'''
# Remove generated files.
#os.remove("test.c")
os.remove("test.so")
