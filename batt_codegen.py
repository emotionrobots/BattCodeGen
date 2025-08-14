#============================================================================================
#
#  batt_codegen.py
#
#  pyBAMM Battery C code generation
#
#  The script outputs batt_model.c and batt_model.h
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
#  Setup the spatial variable to define Mesh
#-------------------------------------------------------------------------------------------
var = pybamm.standard_spatial_vars
var_pts = {var.x_n: 30, var.x_s: 30, var.x_p: 30, var.r_n: 10, var.r_p: 10}
mesh = pybamm.Mesh(geometry, model.default_submesh_types, var_pts)


#-------------------------------------------------------------------------------------------
#  CasADi code-generation must discretize the model first 
#-------------------------------------------------------------------------------------------
disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
disc.process_model(model)


#-------------------------------------------------------------------------------------------
#
#  Generate C code using CasADi
#
#  API: 
#     def generate(self, filename, variable_names, input_parameter_order=None, cg_options=None)
#
#        Parameters
#        ----------
#        filename : str
#            Name of the file to which to save the code
#
#        variable_names : list
#            Variables to be exported alongside the model structure
#
#        input_parameter_order : list (optional)
#            Order in which the input parameters should be stacked.
#            If input_parameter_order=None and len(self.input_parameters) > 1, a
#            ValueError is raised (this helps to avoid accidentally using the wrong
#            order)
#
#        cg_options : dict (optional)
#            Options to pass to the code generator.
#            See https://web.casadi.org/docs/#generating-c-code
#       
#-------------------------------------------------------------------------------------------
opts = dict(main=False, with_header=True, indent=3, verbose=True) 
#var_list = ["Terminal voltage [V]"]
var_list = list(model.variables.keys()) 
model.generate("batt_model.c", 
               var_list, 
               input_parameter_order=[], 
               cg_options=opts,
               write_names=True)




