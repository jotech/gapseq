from cobra.flux_analysis import gapfill
from cobra.flux_analysis.gapfilling import GapFiller
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
#from cobra import Model, Reaction, Metabolite
import pandas
import sys
import re
import os
import gapfill_reference
import seed_fix

if len(sys.argv) != 3:
  print("usage:", sys.argv[0], "file.xml reaction.lst")
  sys.exit(0)

print " ".join(sys.argv)

#seedrDB = pandas.read_csv("/home/jo/uni/gapseq/dat/seed_reactions.tsv", sep="\t")
#seedmDB = pandas.read_csv("/home/jo/uni/gapseq/dat/seed_metabolites.tsv", sep="\t")
dir = os.path.dirname(__file__)
substances=pandas.read_csv(dir+"/../dat/sub2pwy.csv", sep=",")

mod = gapfill_reference.repair_mass_balance(read_sbml_model(sys.argv[1]), verbose=False)
with open(sys.argv[2]) as f:
    lines = f.read()
newR = lines.rstrip("\n").split(" ") # remove linebreak at file end

#
# Create reference model with reactions pool
#

refmod = gapfill_reference.get_reference(mod, newR, delete_unbalanced=True, verbose=False)
refmod = seed_fix.seed_fix(refmod, rmLoops=False)


# Add some reactions, which have evidence by sequence and which are not covered by gapfilling
NRea = len(mod.reactions)
NeededReactions = ["rxn10122_c0", # respiratory complex I
                   "rxn13689_c0", # respiratory complex III
                   "rxn10043_c0"]  # respiratory complex IV
for r in NeededReactions:
    if r in refmod.reactions and r not in mod.reactions:
        mod.add_reaction(refmod.reactions.get_by_id(r))
print "Added needed reactions from reference:", len(mod.reactions) - NRea
mod = seed_fix.seed_fix(mod)


ignore = ["H", "Li", "Na", "Mg", "K", "Ca", "P", "Fe", "Cl", "Zn", "Mn", "Mo", "Se", "Co", "Cu", "Ni", "W", "H2O"]
# prepare model with minimal medium
minmed = ["EX_cpd00007_e0", "EX_cpd00001_e0", # o2, h2o
          "EX_cpd00013_e0", # nh3
          "EX_cpd00048_e0", # sulfate
          # minerals
          "EX_cpd00009_e0", "EX_cpd00012_e0", "EX_cpd00030_e0", "EX_cpd00034_e0", "EX_cpd00058_e0", "EX_cpd00063_e0", "EX_cpd00067_e0", "EX_cpd00099_e0", "EX_cpd00149_e0", "EX_cpd00205_e0", "EX_cpd00254_e0", "EX_cpd10515_e0", "EX_cpd00971_e0", "EX_cpd01012_e0", "EX_cpd10516_e0", "EX_cpd11574_e0"]
# try to find a carbon source
cs = ["EX_cpd00027_e0", # glucose
      "EX_cpd00020_e0", # pyruvate
      "EX_cpd00106_e0", # fumarat
      "EX_cpd00179_e0"] # maltose
found = [c for c in cs if c in [ex.id for ex in mod.exchanges]]
if len(found) > 0:
    carbon = found[0]
    print "Try doing biomass gapfilling with carbon source", carbon, mod.reactions.get_by_id(carbon).name
    minmed += [carbon]
else:
    found2= substances["exid_seed"].loc[substances["exid_seed"].isin([ex.id for ex in mod.exchanges])].tolist()
    if len(found2) > 0:
        carbon  = found2[0]
        print "Try doing biomass gapfilling with carbon source", carbon, mod.reactions.get_by_id(carbon).name
        minmed += [carbon]
    else:
        print "Error: no carbon source found for organism"
        sys.exit(0)


def set_medium(model, medium):
    mod = model.copy()
    for ex in mod.exchanges:
        ex.lower_bound=0

    for m in medium: # set medium
        if m not in mod.reactions:
            #ex = "EX_"+m+"_LPAREN_e_RPAREN_" # vmh
            ex = "EX_"+m+"_e0"
        else:
            ex = m
        #ex = m
        if ex in mod.reactions:
            mod.reactions.get_by_id(ex).lower_bound=-1000
        else:
            print "compound not found in model:", m
    return(mod)

def checkProduction(model, metabolite):
    modtmp = set_medium(model, minmed)
    #refmod.reactions.rxn05115.upper_bound=1000
    #modnew.add_reactions([refmod.reactions.rxn00533])

    newM = []
    for m in refmod.metabolites:
        if m not in modtmp.metabolites:
            newM.append(m.copy())
    modtmp.add_metabolites(newM)

    m = modtmp.metabolites.get_by_id(metabolite)
    modtmp.objective = modtmp.add_boundary(m, type='demand')
    sol = modtmp.slim_optimize()
    if round(sol,12) > 0:
        print m.id, m.name, "can be produced:", sol
    else:
        print m.id, m.name, "CANNOT be produced:", sol
        print "\t --> try to find reactions for gapfilling"
        #print gapfill(modtmp, refmod)   
        filler = GapFiller(modtmp, refmod, integer_threshold=1e-32).fill()
        print filler
        for r in filler[0]:
            print "\t", r.id, r.reaction
        return filler

def checkBiomass(model, fill_gaps=True):
    modmin = set_medium(model, minmed)
    if fill_gaps:
        modnew = model.copy() # will be returned
    newM = []
    for m in refmod.metabolites:
        if m not in modmin.metabolites:
            newM.append(m.copy())
    # add all metabolites from reference (needed for gapfilling => otherwise reactions sometimes not found!)
    modmin.add_metabolites(newM) 
    
    Call=0; Cwork=0; Cfix=0
    fixLater = []
    fail = []
    #relMet = ["cpd00288_c0"]
    allMet = model.reactions.bio1.reactants #+ [modmin.metabolites.get_by_id(m) for m in relMet]
    #allMet = model.reactions.Biomass.reactants # vmh
    #bioaa = [m for m in model.metabolites if m.name.startswith("L_") and m in model.reactions.bio1.metabolites]
    #for m in bioaa:
    for m in allMet:
        if m.formula in ignore: # ignore minerals, etc
            continue
        Call += 1
        modtmp = modmin.copy()
        modtmp.objective = modtmp.add_boundary(m, type='demand')
        obj = modtmp.slim_optimize()
        if obj > 0:
            #print "\t => could be produced:"
            Cwork += 1
        elif fill_gaps: # gapfilling
            print obj, m.id, m.name
            try:
                #gapsol = gapfill(modtmp, refmod) # integer_threshold=1e-06
                gapsol = GapFiller(modtmp, refmod, integer_threshold=1e-16).fill()
            except RuntimeError:
                print "\t => Runtime error: lowering the integer_threshold?"
                fixLater.append(m)
                continue
            except:
                print "\t => failed:", sys.exc_info()[0]
                fail.append(m)
                continue
            if len(gapsol[0]) > 0:
                Cfix += 1
                print "\t => could be fixed:", ",".join([r.id for r in gapsol[0]])
                modnew.add_reactions([r for r in gapsol[0] if r not in modnew.reactions])
    print "\nTotal compounds:", Call, "\t can be produced:", Cwork, "\t could be fixed", Cfix, "\t altogether:", round(100*float(Cwork+Cfix)/Call,1), "%", " (before:",round(100*float(Cwork)/Call,1),"% )"
    if fill_gaps:
        return(modnew)


modnew = checkBiomass(mod, fill_gaps=True)
fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
write_sbml_model(modnew, filename=fileID+"_gapfilled.xml")
