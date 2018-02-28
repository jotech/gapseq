from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from cobra.flux_analysis.gapfilling import GapFiller
from pandas import read_csv
from cobra import Model, Reaction, Metabolite
import re
import os
import sys

if len(sys.argv) < 2 or len(sys.argv) > 3:
  print("usage:", sys.argv[0], "file.xml", "[agora_organism.xml]")
  print("agora_organism can be used to get biomass from agora organism and get reactions to produce biomass with rich medium")
  sys.exit(0)

dir = os.path.dirname(__file__)
mod = read_sbml_model(sys.argv[1])

dic = read_csv(dir+"/../dat/SEED2VMH_translation.csv", header=None, names=["seed","vmh"])
dic["seed"] = dic["seed"].str.replace("\\(e\\)","") # needed to match exchange reactions

# find biomass
if len(sys.argv) == 3:
    refmod = read_sbml_model(sys.argv[2])
    for m in refmod.metabolites: # compartment tag inconsistency
        m.id = re.sub("_(.)$","[\\1]", m.id)
    for ex in refmod.exchanges:
        ex.id = ex.id.replace("_LPAREN_","(").replace("_RPAREN_",")").replace("_DASH_","_") # compartment tag for exchange reactions
    print "Read reference model:", refmod.name
    biomassDB = read_csv(dir+"/../dat/agora_biomass.csv", header=None, names=["org","biomass"])
    try:
        biomass = biomassDB["biomass"][biomassDB["org"]==refmod.id].values[0]
    except IndexError:
        biomass = "biomass0"
else:
    biomass = "biomass0"

print "Take biomass reaction:", biomass


dicSeed = dic["seed"].tolist()
reactions = [re.sub("_.0$","", r.id) for r in mod.reactions]
matched = [dic["vmh"][dic["seed"]==r].tolist()[0] for r in reactions if r in dicSeed]

vmh_rea = read_csv(dir+"/../dat/vmh_reactions.csv")
vmh_met = read_csv(dir+"/../dat/vmh_metabolites.csv")

newmod = Model('NEW')
newmod.compartments = {'c':'cytosol','e':'extracellular'}

matched += [biomass, "EX_biomass(e)"]
#matched += [biomass, "EX_biomass(e)", "26DAPt2r", "HEMEti", "GCALDt", "FE2abc"]
for r in matched:
    if r in newmod.reactions:
        continue
    line = vmh_rea[vmh_rea["abbreviation"]==r]
    formula = line["formula"].values[0]
    name  = line["description"].values[0]
    rea = Reaction(id=r,name=name)
    newmod.add_reactions(rea)
    rea.build_reaction_from_string(formula, verbose=0)

for m in newmod.metabolites:
    mid = re.sub("\\[.\\]","",m.id)
    line = vmh_met[vmh_met["abbreviation"]==mid]
    m.name = line["fullName"].values[0]
    m.charge = line["charge"].values[0]
    m.formula = line["chargedFormula"].values[0]

# set objective
newmod.objective=biomass
newmod.repair()
#newmod.optimize()

sol = newmod.optimize()
if 'refmod' in globals() and (sol.status != "optimal" or round(sol.objective_value,6) == 0):
    print "no biomass production possible try to fix"
    tmpmod = newmod.copy()
    for ex in tmpmod.exchanges:
        ex.lower_bound=-1000
    try:
        gapsol = GapFiller(tmpmod, refmod, demand_reactions=False, integer_threshold=1e-16).fill()
    except RuntimeError:
        pass
    except:
        pass
    if len(gapsol[0]) > 0:
        print "Add reactions from reference to be able to produce biomass:", ",".join([r.id for r in gapsol[0]])
        newmod.add_reactions([r for r in gapsol[0]])
fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
write_sbml_model(newmod, filename=fileID+"_vmh.xml")
