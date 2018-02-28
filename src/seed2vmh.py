from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from pandas import read_csv
from cobra import Model, Reaction, Metabolite
import re
import os
import sys

if len(sys.argv) < 2 or len(sys.argv) > 3:
  print("usage:", sys.argv[0], "file.xml", "[agora_organism]")
  print("agora_organism can be used to get biomass from agora organism")
  sys.exit(0)

dir = os.path.dirname(__file__)
mod = read_sbml_model(sys.argv[1])

dic = read_csv(dir+"/../dat/SEED2VMH_translation.csv", header=None, names=["seed","vmh"])
dic["seed"] = dic["seed"].str.replace("\\(e\\)","") # needed to match exchange reactions

# find biomass
if len(sys.argv) == 3:
    biomassDB = read_csv(dir+"/../dat/agora_biomass.csv", header=None, names=["org","biomass"])
    try:
        biomass = biomassDB["biomass"][biomassDB["org"]==sys.argv[2]].values[0]
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
#newmod.optimize()

newmod.repair()
fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
write_sbml_model(newmod, filename=fileID+"_vmh.xml")
