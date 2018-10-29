from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from pandas import read_csv
from cobra import Model, Reaction, Metabolite
import re
import os
import sys
import collections
from helper import addFromSeedDB

if len(sys.argv) < 2:
  print("usage:", sys.argv[0], "file.xml")
  sys.exit(0)

dir = os.path.dirname(__file__)
print("working on:", sys.argv[1])
mod = read_sbml_model(sys.argv[1])

seed_rea = read_csv(dir+"/../dat/seed_reactions_corrected.tsv", sep="\t")
#seed_met = read_csv(dir+"/../dat/seed_metabolites.tsv")

fix_dir = 0
fix_stoich = 0
fix_rm = 0
na_rm = 0
new_met = 0
obsolete = 0 # obsolete
comp = {'0':'_c0', '1':'_e0', '2':'_p0'}
obj = [r for r in mod.reactions if r.objective_coefficient!=0] # dont remove the objective 
for r in set(mod.reactions).difference(mod.exchanges + obj):
    rshortID = re.sub("_.0$","", r.id)
    line = seed_rea[seed_rea["id"]==rshortID]
    if len(line.index) == 0: # reactions which are marked as duplicated but are not considered yet (i.e. status of not.assessed)
        #print r
        obsolete += 1
        r.remove_from_model()
        continue
    
    # check curration status of reaction
    status = line["gapseq.status"].values[0]
    if status == "removed":
        fix_rm += 1
        r.remove_from_model()
        continue
    
    if status == "not.assessed": # reactions which are not assessed now
        na_rm += 1
        r.remove_from_model()
        continue

    # check if stoichiometry has to be updated
    formula = line["stoichiometry"].values[0]
    org = {m.id:float(r.metabolites[m]) for m in r.metabolites}
    new = {m.split(":")[1]+comp[m.split(":")[2]]:float(m.split(":")[0]) for m in formula.split(";")} 
    if org != new:
        fix_stoich += 1
        try:
            for m in new:
                if m not in mod.metabolites:
                    mod.add_metabolites(Metabolite(m))
                    new_met += 1
            replace = {mod.metabolites.get_by_id(m):new[m] for m in new}
        except:
            print( "ERROR", sys.exc_info()[0])
            print( r)
            print( org)
            print( new)
            sys.exit(1)
        r.subtract_metabolites(r.metabolites)
        r.add_metabolites(replace)
        #if round(mod.slim_optimize(),6) == 0:
        #    print "Cannot grow anymore after fixing stoichiometry:", r
        #    sys.exit(0)
    
    # check reversibility
    reversibility = line["reversibility"].values[0]
    if reversibility=='=':
        lb =-1000
        ub = 1000
    if reversibility=='>':
        lb = 0
        ub = 1000
    if reversibility=='<':
        lb =-1000
        ub =  0
    if r.lower_bound != lb or r.upper_bound != ub:
        fix_dir += 1
        r.lower_bound = lb
        r.upper_bound = ub
        #if round(mod.slim_optimize(),6) == 0:
        #    print "Cannot grow anymore after fixing reversibility:", r
        #    sys.exit(0)

print( "\tfixed reversibility:", fix_dir)
print( "\tfixed stoichiometry:", fix_stoich)
print( "\tremoved reactions  :", fix_rm)
print( "\tNot in DB:", obsolete)
print( "\tNot assessed:", na_rm)
print( "\tNew metabobolites:", new_met)

#print mod.summary()

fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
write_sbml_model(mod, filename=fileID+"_corrected.xml")
