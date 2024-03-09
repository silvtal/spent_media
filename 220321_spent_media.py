#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:19:35 2021

@author: shiru
"""
from carveme.reconstruction.utils import medium_to_constraints
from reframed.core.transformation import disconnected_metabolites
from reframed.solvers import solver_instance
from reframed.solvers.solver import VarType
from reframed.solvers.solution import Status
from reframed import load_cbmodel, FBA, CAFBA, Environment, plot_flux_envelope
from reframed import *
from carveme.reconstruction.utils import load_media_db
import  os, cplex


# =============================================================================
# Creo dos modelos y dos soluciones
# =============================================================================
os.chdir("/home/shiru/AAA/cit_check_TRUE/Node35562")
flavor= "fbc2"
max_uptake=10
medium="M9[cit]"
media=list(medium)
media_db = load_media_db("../../my_media.tsv", compound_col="compound")
init_env = Environment.from_compounds(media_db[medium])


# Entero
f="691423.xml"

merged_model = load_cbmodel(f,flavor="fbc2")

# Pseudo
spent="../Node28866/4419276.xml"
spent_model = load_cbmodel(spent, flavor=flavor)
solver0 = solver_instance(spent_model)

# FRom gapfilling.py:
for medium_name in media:
    if medium_name in media_db:
        compounds = set(media_db[medium_name])
        constraints = medium_to_constraints(merged_model, compounds, max_uptake=max_uptake, inplace=False, verbose=False)

        if spent_model:
            constraints0 = medium_to_constraints(spent_model, compounds, max_uptake=max_uptake, inplace=False, verbose=False)
            print(constraints0)
            # for r_id in spent_model.get_exchange_reactions(): # ESTA es la clave
            #     if r_id in constraints:
            #         sol = FBA(spent_model, objective={r_id: 1}, constraints=constraints0, solver=solver0, get_values=False)
            #         if sol.fobj > ABSTOL:
            #             constraints[r_id] = (-max_uptake, None)
            #             print("added", r_id[5:-2], "to", medium_name)
                        
spent_model.get_exchange_reactions()
