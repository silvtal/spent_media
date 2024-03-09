#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 13:21:21 2021

@author: shiru
"""

import os
from reframed import load_cbmodel, FBA, CAFBA, Environment, plot_flux_envelope
from reframed import *
from carveme.reconstruction.utils import load_media_db
import pandas as pd
import itertools
import stats
import copy

# input
wd="/home/silvia/AAA/2021-04-20_CAFBA_M9[caca]"
mediadb = "/home/silvia/AAA/my_media.tsv"
medium="M9[glc]"
to_remove="glc__D"
# node 1: we check their growth;; node 2: spent medium
reverse=False  #(reverse true si cambio node1 por node2 xd)
#node2 = "/home/silvia/AAA/glc_check_TRUE/models/Node35562" #entero
#node1 = "/home/silvia/AAA/glc_check_TRUE/models/Node27828" #pseudo
node2="/home/silvia/temp_NodeE_E2"
node1="/home/silvia/temp_NodePs_E2"
outputname = "E2_results_CAFBA_for_all_M9[sin_glc]_con_caca_MEDIA_DE_RATIOS_con_caca_DE_ENTERO" #edited by hand later to include everything
cores=12

# definitions
def load_xml(folder,wd):
    models = []
    for sbml in os.listdir(folder):
      if ".xml" in sbml:
        models.append(folder+"/"+sbml)
    return models
        
def create_spare_medium(f,medium,media_db,init_env=None, outputdir=".", outputname="tempmedium.tsv", include_original=True, to_remove=None):
    """ Medium creation.
    
    =============================================================================
    There are three different types of pre-defined boundary reactions: exchange, 
    demand, and sink reactions. All of them are unbalanced pseudo reactions, that 
    means they fulfill a function for modeling by adding to or removing metabolites
    from the model system but are not based on real biology. An exchange reaction 
    is a reversible reaction that adds to or removes an extracellular metabolite 
    from the extracellular compartment. A demand reaction is an irreversible 
    reaction that consumes an intracellular metabolite. A sink is similar to an 
    exchange but specifically for intracellular metabolites, i.e., a reversible 
    reaction that adds or removes an intracellular metabolite.
    =============================================================================
    
    =============================================================================
    --> i need all exchange ones except sinks
    --> i only care about the external metabolites that are being produced/dumped
    --> i need those external metabs that
          1) have a exchange reaction in my studied model
          2) and whose exchange reaction has a positive flux (?)
    =============================================================================
    
    """
    # load and solve the spent_model
    s_model = load_cbmodel(f,flavor="fbc2")   
    if init_env==None:
        init_env = Environment.from_compounds(media_db[medium])
    
    init_env.apply(s_model)
    solution = CAFBA(s_model)  
    
    # select all positive flux reactions (?)
    df=solution.to_dataframe()
    positive_fluxes= df[df["value"]>0]
    
    # pick all exchange reactions except sinks
    my_r= [r for r in s_model.get_exchange_reactions() if "sink" not in r]
    
    spare=[r for r in list(positive_fluxes.index) if r in my_r]
    
    # spare
    # Out[81]: 
    # ['R_EX_h2o_e',
    #  'R_EX_succ_e',
    #  'R_EX_co2_e',
    #  'R_EX_fe2_e',
    #  'R_EX_glyclt_e',
    #  'R_EX_urea_e']
    
    # get metabolites from those reactions
    my_m=["".join(r.split("_")[2:-1]) for r in spare]

    # create file
    if include_original:
        if to_remove==None:
            my_m= media_db[medium]+my_m
        else:
            medcopy=copy.copy(media_db[medium])
            medcopy.remove(to_remove)
            my_m= medcopy+my_m
    
    temp_media=pd.DataFrame({
        "medium" : ["temp"]*len(my_m),
        "description" : ["spent medium without original"]*len(my_m),
        "compound" : my_m
        })
    
    temp_media.to_csv( outputdir+"/"+outputname,sep="\t")
    
def ratio_calculator(models, s_models, medium, media_db, outputdir=".",outputname=outputname, include_original=False, to_remove=None, reverse=False):
    """
    Calculates all the ratios and the ratio mean
    
    models   -> the ones which growth we are checking.
    s_models -> the ones which spend media we are using
    """
    
    pairs = list(itertools.product(models1, s_models))

    if os.path.isfile(outputdir+"/"+outputname): # if exists:
        exists=True
        try: 
            already_done=pd.read_csv(outputdir+"/"+outputname)
        except:
            print("failed")
    else: # if doesnt exist yet
        exists=False
        header=pd.DataFrame([["Medio","Modelo","spent_model","Ratio"]])
        header.to_csv(outputdir+"/"+outputname,index=False,header=False)
        
    # This first loop leaves out already done pairs.
    todo_pairs = []
    for f,s in pairs:  
        if exists and (f.split("/")[-1]+s.split("/")[-1] in list(already_done.iloc[:,1]+already_done.iloc[:,2])):
            continue
        else:
            todo_pairs.append((f,s))
    pairs=todo_pairs
    
    for f,s in pairs:  
        # Compute FIRST growth rate (the one that produces the spent medium)
        model = load_cbmodel(s,flavor="fbc2")
        original_medium = load_media_db(mediadb, compound_col="compound")
 
        init_env = Environment.from_compounds(original_medium[medium])
        
        init_env.apply(model)
        solution = CAFBA(model,   objective="Growth")    
        
        if solution.fobj == 0:
            print("The first species ("+s+") does not grow or produce any extracellular compounds. Moving to the next pair...")
            continue
        first_gr = copy.deepcopy(solution.fobj)
        
        # Compute SECOND growth rate (the one that grows in spent medium)
        model = load_cbmodel(f,flavor="fbc2")
        create_spare_medium(s, medium, media_db,outputdir=outputdir,include_original=include_original,to_remove=to_remove) # !!
        temp_medium = load_media_db(outputdir+"/tempmedium.tsv", compound_col="compound")
        init_env = Environment.from_compounds(temp_medium["temp"])
        
        init_env.apply(model)
        solution = CAFBA(model,   objective="Growth")    
        
        second_gr = copy.deepcopy(solution.fobj)
        
        if reverse:
            if second_gr==0:
                continue
            else:
                ratio_of_pair = first_gr/second_gr
        else:
            ratio_of_pair = second_gr/first_gr

        results=[medium,f.split("/")[-1],s.split("/")[-1],ratio_of_pair]
        all_data=pd.DataFrame([results])
        all_data.to_csv(outputdir+"/"+outputname,index=False,mode="a",header=False)

# run
os.chdir(wd)

models1 = load_xml(node1,wd)
models2 = load_xml(node2,wd)

media_db = load_media_db(mediadb, compound_col="compound")

init_env = Environment.from_compounds(media_db[medium])

ratio_calculator(models1,models2, medium, media_db, outputdir=wd,include_original=True, to_remove=to_remove, reverse=reverse)
