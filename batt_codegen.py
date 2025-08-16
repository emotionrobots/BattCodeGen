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
#  Create DFN model with required configurations 
#-------------------------------------------------------------------------------------------
model = pybamm.lithium_ion.DFN(options={
            "thermal": "lumped",                         # lumped thermal ODE
            "SEI": "ec reaction limited",                # SEI -> contributes to LLI
            "loss of active material": "reaction-driven",  # LAM mechanism
            #"loss of active material": "stress-driven",  # LAM mechanism
            # (You can switch to "reaction-driven" or tuple per electrode if desired)
        })
'''
model = pybamm.lithium_ion.DFN()
'''

#-------------------------------------------------------------------------------------------
# Add explicit anode/cathode potentials as variables (convenient for codegen)
#-------------------------------------------------------------------------------------------
model.variables["Anode potential [V]"] = model.variables[
    "Negative electrode surface potential difference at separator interface [V]"
]

model.variables["Cathode potential [V]"] = model.variables[
    "Positive electrode surface potential difference at separator interface [V]"
]


#-------------------------------------------------------------------------------------------
#  Use default geometry
#-------------------------------------------------------------------------------------------
geometry = model.default_geometry


#-------------------------------------------------------------------------------------------
#  Setup default parameters but tweak capacity and cut-off voltages
#-------------------------------------------------------------------------------------------
params = model.default_parameter_values
params.update({
        "Nominal cell capacity [A.h]": 4.0,
        "Upper voltage cut-off [V]": 4.25,
        "Lower voltage cut-off [V]": 3.0, 
    })


#-------------------------------------------------------------------------------------------
#  Mark current and ambient temperature as model simulation inputs  
#-------------------------------------------------------------------------------------------
params.update({
    "Current function [A]": pybamm.InputParameter("Current function [A]"),
    "Ambient temperature [K]": pybamm.InputParameter("Ambient temperature [K]"),
    },
    check_already_exists=False)


#-------------------------------------------------------------------------------------------
#  Apply parameters to geoemetry and model
#-------------------------------------------------------------------------------------------
params.process_geometry(geometry)
params.process_model(model)


#-------------------------------------------------------------------------------------------
#  Setup the spatial variable to define Mesh for discretisation
#-------------------------------------------------------------------------------------------
var = pybamm.standard_spatial_vars
var_pts = {var.x_n: 30, var.x_s: 30, var.x_p: 30, var.r_n: 10, var.r_p: 10}
mesh = pybamm.Mesh(geometry, model.default_submesh_types, var_pts)

disc = pybamm.Discretisation(mesh, model.default_spatial_methods)
disc.process_model(model)


#-------------------------------------------------------------------------------------------
#  Build a concise list of variables to export from C 
#-------------------------------------------------------------------------------------------
export_vars = [
    "Voltage [V]",                                          # terminal voltage
    "Current [A]",                                          # sign convention (out of cell)
    "Discharge capacity [A.h]",                             # for SOC calc in C
    "Total lithium capacity [A.h]",                         # for SOC calc in C
    "X-averaged cell temperature [K]",                      # lumped temperature
    "Anode potential [V]",                                  # anode potential
    "Cathode potential [V]",                                # cathode poential
    "Loss of lithium inventory [%]",                        # LLI
    "Loss of active material in negative electrode [%]",    # anaode LAM 
    "Loss of active material in positive electrode [%]",    # cathode LAM
]


#-------------------------------------------------------------------------------------------
# Order of runtime inputs given to generated C entry points:
#-------------------------------------------------------------------------------------------
inputs = ["Current function [A]", "Ambient temperature [K]"]



#-------------------------------------------------------------------------------------------
#
#  Generate C code using CasADi
#
#  API: 
#     generate(filename, variable_names, input_parameter_order=None, cg_options=None)
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
opts = {
        "main": False, 
        "with_header": True, 
        "with_mem": False,
        "casadi_real": "double",
        "casadi_int" : "long long int",
        "indent": 3, 
        "verbose": True 
}

model.generate(
        filename="batt_model.c", 
        variable_names=export_vars, 
        input_parameter_order=inputs, 
        cg_options=opts)




