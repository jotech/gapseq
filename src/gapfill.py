from cobra.flux_analysis import gapfill
from cobra.flux_analysis.gapfilling import GapFiller
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
#from cobra import Model, Reaction, Metabolite
#import pandas
import sys
import re
import os
import gapfill_reference

if len(sys.argv) != 3:
  print("usage:", sys.argv[0], "file.xml reaction.lst")
  sys.exit(0)

#seedrDB = pandas.read_csv("/home/jo/uni/gapseq/dat/seed_reactions.tsv", sep="\t")
#seedmDB = pandas.read_csv("/home/jo/uni/gapseq/dat/seed_metabolites.tsv", sep="\t")

mod = gapfill_reference.repair_mass_balance(read_sbml_model(sys.argv[1]), verbose=False)
with open(sys.argv[2]) as f:
    lines = f.read()
newR = lines.rstrip("\n").split(" ") # remove linebreak at file end

#
# Create reference model with reactions pool
#

refmod = gapfill_reference.get_reference(mod, newR, delete_unbalanced=True, verbose=False)


# CORRECTION IN SEED DATABASE
if "rxn05115_c0" in refmod.reactions:
    refmod.reactions.rxn05115_c0.upper_bound=1000 # wrong direction
if "rxn03031_c0" in mod.reactions:
    mod.reactions.rxn03031_c0.lower_bound=-1000 # bidirectional http://www.rhea-db.org/reaction?id=17325


ignore = ["H", "Li", "Na", "Mg", "K", "Ca", "P", "Fe", "Cl", "Zn", "Mn", "Mo", "Se", "Co", "Cu", "Ni", "W", "H2O"]
# prepare model with minimal medium
minmed = ["cpd00027", "cpd00007", "cpd00001", # glc, o2, h2o
          #"cpd00020", # pyruvate
          "cpd00106", # fumarat # myb10
          "cpd00013", # nh3
          "cpd00048", # sulfate
          # minerals
          "cpd00009", "cpd00012", "cpd00030", "cpd00034", "cpd00058", "cpd00063", "cpd00067", "cpd00099", "cpd00149", "cpd00205", "cpd00254", "cpd10515", "cpd00971", "cpd01012", "cpd10516", "cpd11574"]

def set_medium(model, medium):
    mod = model.copy()
    for ex in mod.exchanges:
        ex.lower_bound=0

    for m in medium: # set medium
        #ex = "EX_"+m+"_LPAREN_e_RPAREN_" # vmh
        ex = "EX_"+m+"_e0"
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
