#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 16:13:18 2021

@author: Silvia
"""

# =============================================================================
# Cada nodo recibirá una lista de parejas de una longitud de 100. Así
# ahorramos en tiempo de copipega: copipegamos solo cada 100 veces.
# Cada vez que TODOS los nodos hayan acabado con sus 100, se guardará a disco,
# ya que mando trabajos en cantidades de 1200 pares y uso 12 nodos.
#
# Si no interrumpo manualmente, el script seguirá: cada nodo cogerá otros 100.
# 
# Si interrumpo, en la siguiente ejecución verá primero qué pares ya han sido
# procesados y los quitará de la lista.
# 
# bucle: divide pairs en chunks_of_pairs de 1200. 
# Parallel: divide chunks_of_pairs en pairlist de 100. Cada nodo de 12 coge 1.
# process: procesa cada pair en pairlist. Devuelve el resultado junto
# =============================================================================


import os
from reframed import load_cbmodel, FBA, CAFBA, Environment, plot_flux_envelope
from carveme.reconstruction.utils import load_media_db
import pandas as pd
import itertools
import stats
import copy
import numpy as np
import random
from joblib import Parallel, delayed

# =============================================================================
# input
# =============================================================================
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
    
    #add a random number to the tempmedia file to make it unique and allow paralellization.
    outputname="n"+str(int(random.random()*3456789))+outputname
    temp_media.to_csv( outputdir+"/"+outputname,sep="\t")
    
    return(outputname)
    
def ratio_calculator_pairs(models, s_models, medium, media_db, outputdir=".",outputname="ratiocalcresults", include_original=False, to_remove=None, reverse=False):
    
    pairs = list(itertools.product(models, s_models))

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
        already_done=[0]
        
    # # This first loop leaves out already done pairs.
    # todo_pairs = []
    # for f,s in pairs:  
    #     if exists and (f.split("/")[-1]+s.split("/")[-1] in list(already_done.iloc[:,1]+already_done.iloc[:,2])):
    #         print("está")
    #         pass
    #     else:
    #         print("no está, me lo apunto")
    #         todo_pairs.append((f,s))
    return(pairs[(len(already_done)-1):-1])

def ratio_calculator(pair, medium, media_db, 
            outputdir=".",outputname="ratiocalcresults",
            include_original=False, to_remove=None, reverse=False):
    """
    Calculates all the ratios and the ratio mean
    
    models   -> the ones which growth we are checking.
    s_models -> the ones which spend media we are using
    """
    f=pair[0]
    s=pair[1]

    # Compute FIRST growth rate (the one that produces the spent medium)
    model = load_cbmodel(s,flavor="fbc2")
    original_medium = load_media_db(mediadb, compound_col="compound")
 
    init_env = Environment.from_compounds(original_medium[medium])
    
    init_env.apply(model)
    solution = CAFBA(model,   objective="Growth")    
    
    if solution.fobj != 0:
        first_gr = copy.deepcopy(solution.fobj)
    
        # Compute SECOND growth rate (the one that grows in spent medium)
        model = load_cbmodel(f,flavor="fbc2")
        tempmediumname = create_spare_medium(s, medium, media_db,outputdir=outputdir,include_original=include_original,to_remove=to_remove) # !!
        temp_medium = load_media_db(outputdir+"/"+tempmediumname, compound_col="compound")
        os.system("rm "+outputdir+"/"+tempmediumname)
        init_env = Environment.from_compounds(temp_medium["temp"])
        
        init_env.apply(model)
        solution = CAFBA(model,   objective="Growth")    
        
        second_gr = copy.deepcopy(solution.fobj)
        
        if reverse:
            if second_gr!=0:
                ratio_of_pair = first_gr/second_gr
                results=[medium,f.split("/")[-1],s.split("/")[-1],ratio_of_pair]
                return(results) 
            else:       
                return([medium,f.split("/")[-1],s.split("/")[-1],"NA"]) # retornar vacío si tal.
        else:
            ratio_of_pair = second_gr/first_gr
            results=[medium,f.split("/")[-1],s.split("/")[-1],ratio_of_pair]
            return(results)

    else:
        print("The first species ("+s+") does not grow or produce any extracellular compounds. Moving to the next pair...")
        return([medium,f.split("/")[-1],s.split("/")[-1],"NA"]) # retornar vacío si tal.

def chunks(L, n): return [L[x: x+n] for x in range(0, len(L), n)]

def process(pairlist, medium, media_db, 
            outputdir=".",outputname="ratiocalcresults",
            include_original=False, to_remove=None, reverse=False):
    results=[]
    for pair in pairlist:
        results.append(ratio_calculator(pair, medium, media_db, 
            outputname=outputname,outputdir=wd,include_original=True, 
            to_remove=to_remove, reverse=reverse))    
    
    return results    
    

# =============================================================================
# run
# =============================================================================
os.chdir(wd)

models1 = load_xml(node1,wd)
models2 = load_xml(node2,wd)

media_db = load_media_db(mediadb, compound_col="compound")

init_env = Environment.from_compounds(media_db[medium])

pairs = ratio_calculator_pairs(models1,models2, medium, media_db, outputdir=wd,outputname=outputname,include_original=True, to_remove=to_remove, reverse=reverse)

for chunk_of_pairs in chunks(pairs,8400):
    all_results = Parallel(n_jobs=14)(delayed(process)(pairlist, medium, media_db, 
            outputname=outputname,outputdir=wd,include_original=True, 
            to_remove=to_remove, reverse=reverse) for pairlist in chunks(chunk_of_pairs,600))
    
    final_results = [item for sublist in all_results for item in sublist]
    all_data=pd.DataFrame(final_results)
    
    all_data.to_csv(wd+"/"+outputname,index=False,mode="a",header=False)
