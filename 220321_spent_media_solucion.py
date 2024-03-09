import  os, cplex
os.chdir("/home/silvia/AAA/cit_check_TRUE/Node35562")

from reframed import load_cbmodel, FBA, CAFBA, Environment, plot_flux_envelope
from reframed import *
from carveme.reconstruction.utils import load_media_db

medium="M9[cit]"
media_db = load_media_db("../../my_media.tsv", compound_col="compound")
init_env = Environment.from_compounds(media_db[medium])


# =============================================================================
# Creo dos modelos y dos soluciones
# =============================================================================
# Entero
f="691423.xml"

model = load_cbmodel(f,flavor="fbc2")
init_env.apply(model)
solution = CAFBA(model)  

# Pseudo
g="../Node28866/4419276.xml"

model2 = load_cbmodel(g,flavor="fbc2")
init_env.apply(model2)
solution2 = CAFBA(model2)  


# A partir de aquí exploro...
for model in [model,model2]:
    model.summary()
    # model.get_compartment_metabolites()
    # model.get_metabolite_consumers()
    
    model.get_external_metabolites()
    
    model.metabolites
    
    plot_flux_envelope(model, model.biomass_reaction, 'R_EX_cit_e')


# =============================================================================
# Clasificar reacciones (el metodo no me funciona, me lo hago yo)
# =============================================================================
my_r=[ r for r in model.reactions.keys() if  model.reactions[r].reaction_type == ReactionType.EXCHANGE]
# model.get_reactions_by_type('')
my_r2= model.get_exchange_reactions()

# =============================================================================
# Obtener metabolitos que se liberan (?) al medio.
# =============================================================================
# model.compartments.keys()
# Out[4]: odict_keys(['C_c', 'C_p', 'C_e'])
my_metab = [m for m in model.metabolites.keys() if model.metabolites[m].compartment == "C_e"]


# Vemos que no es lo mismo
# len(my_metab)
# Out[27]: 205

# len(my_r)
# Out[28]: 199
# len(model.get_exchange_reactions())
# Out[36]: 199


# =============================================================================
# There are three different types of pre-defined boundary reactions: exchange, 
# demand, and sink reactions. All of them are unbalanced pseudo reactions, that 
# means they fulfill a function for modeling by adding to or removing metabolites
# from the model system but are not based on real biology. An exchange reaction 
# is a reversible reaction that adds to or removes an extracellular metabolite 
# from the extracellular compartment. A demand reaction is an irreversible 
# reaction that consumes an intracellular metabolite. A sink is similar to an 
# exchange but specifically for intracellular metabolites, i.e., a reversible 
# reaction that adds or removes an intracellular metabolite.
# =============================================================================

# =============================================================================
# --> las sink no me interesan, solo me interesan metabs que están fuera
# --> los metabolitos externos solo me interesan si se estan produciendo/liberando
# --> tengo que coger aquellos metabs externos 
#       1) que tengan reaccion de intercambio
#       2) y cuya reaccion de intercambio(no sink) tenga un flujo positivo.
# =============================================================================

# =============================================================================
# explorar SOLUCION CAFBA
# =============================================================================
df=solution.to_dataframe()

positive_fluxes= df[df["value"]>0]

solution.get_metabolites_turnover(model) #Metabolic turnover is the sum of all 
# biochemical reactions – the metabolic network – within a living cell that are 
# involved in the conversion of substrates, particularly carbon, nitrogen and 
# energy sources into Gibbs free energy, secreted metabolites and cell mass. 
# Analysis of metabolic turnover involves identification of which reactions are 
# active

# el máximo es 10 por defecto y vemos que es 10 en citrato. Esto es el grado de 
# CONSUMO, NO me interesa


