from cobra.flux_analysis import gapfill
from cobra.flux_analysis.gapfilling import GapFiller
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from cobra import Model, Reaction, Metabolite
from pandas import read_csv
import sys
import re
import os
import gapfill_reference
from helper import getReferenceMod
#import seed_fix

if len(sys.argv) < 3 or len(sys.argv) > 4:
  print("usage:", sys.argv[0], "file.xml reaction.lst [transporter.lst]")
  sys.exit(0)

print " ".join(sys.argv)

#seedrDB = pandas.read_csv("/home/jo/uni/gapseq/dat/seed_reactions.tsv", sep="\t")
#seedmDB = pandas.read_csv("/home/jo/uni/gapseq/dat/seed_metabolites.tsv", sep="\t")
dir = os.path.dirname(__file__)
seed_rea = read_csv(dir+"/../dat/seed_reactions_corrected.tsv", sep="\t")
substances=read_csv(dir+"/../dat/sub2pwy.csv", sep=",")

#mod = gapfill_reference.repair_mass_balance(read_sbml_model(sys.argv[1]), verbose=False)
mod = read_sbml_model(sys.argv[1])
with open(sys.argv[2]) as f:
    lines = f.read()
newR = lines.rstrip("\n").split(" ") # remove linebreak at file end

if len(sys.argv) == 4:
    with open(sys.argv[3]) as f:
        lines = f.read()
    newR += lines.rstrip("\n").split(" ") # remove linebreak at file end



#newR.extend(["rxn08020", "rxn12239", "rxn08131", "rxn00957"]) # TODO: workaround!
#print(newR)
#sys.exit()

#
# Create reference model with reactions pool
#


refmod = read_sbml_model(dir+"/../dat/seedMetaModel.xml") 
#refmod = getReferenceMod(newR)



# Add some reactions, which have evidence by sequence and which are not covered by gapfilling
NRea = len(mod.reactions)
NeededReactions = ["rxn10122_c0", # respiratory complex I
                   "rxn13689_c0", # respiratory complex III
                   "rxn10043_c0"]  # respiratory complex IV
for r in NeededReactions:
    if r in refmod.reactions and r not in mod.reactions:
        mod.add_reaction(refmod.reactions.get_by_id(r))
print "Added needed reactions from reference:", len(mod.reactions) - NRea
#mod = seed_fix.seed_fix(mod)


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


def set_medium(model, medium, verbose=True):
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
            if verbose:
                print "compound not found in model:", m
    return(mod)

def checkProduction(model, metabolite, fill=True, addMets=False):
    modtmp = set_medium(model, minmed, verbose=False)
    
    if addMets:
        newM = []
        for m in refmod.metabolites:
            if m not in modtmp.metabolites:
                newM.append(m.copy())
        modtmp.add_metabolites(newM)
    
    rea1 = Reaction("Refresh_ATP", lower_bound=0, upper_bound=1000)
    rea1.add_metabolites({modtmp.metabolites.cpd00002_c0:-1, modtmp.metabolites.cpd00001_c0:-1, modtmp.metabolites.cpd00008_c0:1,modtmp.metabolites.cpd00009_c0:1, modtmp.metabolites.cpd00067_c0:1})
    modtmp.add_reaction(rea1)

    rea2 = Reaction("Refresh_NADH", lower_bound=0, upper_bound=1000)
    rea2.add_metabolites({modtmp.metabolites.cpd00004_c0:-1, modtmp.metabolites.cpd00067_c0:1, modtmp.metabolites.cpd00003_c0:1})
    modtmp.add_reaction(rea2)
    
    rea3 = Reaction("Refresh_NADPH", lower_bound=0, upper_bound=1000)
    rea3.add_metabolites({modtmp.metabolites.cpd00005_c0:-1, modtmp.metabolites.cpd00067_c0:1, modtmp.metabolites.cpd00006_c0:1})
    modtmp.add_reaction(rea3)

    modtmp.add_boundary(modtmp.metabolites.cpd00012_c0, type='demand') # ppi
    modtmp.add_boundary(modtmp.metabolites.cpd00009_c0, type='demand') # pi
    modtmp.add_boundary(modtmp.metabolites.cpd00067_c0, type='demand') # h+

    #modtmp.add_boundary(modtmp.metabolites.cpd00010_c0, type='exchange', lb=-1000)

    #rea4 = Reaction("Refresh_Pi", lower_bound=0, upper_bound=1000)
    #rea4.add_metabolites({modtmp.metabolites.cpd00009_c0:-1})
    #modtmp.add_reaction(rea4)


    m = modtmp.metabolites.get_by_id(metabolite)
    modtmp.objective = modtmp.add_boundary(m, type='demand')
    #sol = modtmp.slim_optimize()
    sol = modtmp.optimize()
    if round(sol.objective_value,12) > 0:
        print m.id, m.name, "can be produced:", sol.objective_value
        return
    else:
        print m.id, m.name, "CANNOT be produced:", sol.objective_value
        print "\t --> try to find reactions for gapfilling"
        #print gapfill(modtmp, refmod)   
        try:
            gapsol = GapFiller(modtmp, refmod, integer_threshold=1e-32).fill()
        except RuntimeError:
            print "\t => Runtime error: lowering the integer_threshold?"
            return
        except:
            print "\t => failed:", sys.exc_info()[0]
            return
        if len(gapsol[0]) > 0:
            print "\t => could be fixed:", ",".join([r.id for r in gapsol[0]])
            if fill:
                model.add_reactions([r for r in gapsol[0] if r not in model.reactions])


def checkBiomass(model, fill_gaps=True, addMets=False):
    #modmin = set_medium(model, minmed, verbose=False)
    if fill_gaps:
        modnew = model.copy() # will be returned
    
        
    if addMets: # add all metabolites from reference (needed for gapfilling => otherwise reactions sometimes not found!)
        newM = []
        for m in refmod.metabolites:
            if m not in modnew.metabolites:
                newM.append(m.copy())
        #modmin.add_metabolites(newM) 
        modnew.add_metabolites(newM) 
    
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
            print m.id, m.name, m.formula
        Call += 1
        
        #modtmp = modmin.copy()
        modtmp = set_medium(modnew, minmed, verbose=False)
    
#        rea1 = Reaction("Refresh_ATP", lower_bound=0, upper_bound=1000)
#        rea1.add_metabolites({modtmp.metabolites.cpd00002_c0:-1, modtmp.metabolites.cpd00001_c0:-1, modtmp.metabolites.cpd00008_c0:1,modtmp.metabolites.cpd00009_c0:1, modtmp.metabolites.cpd00067_c0:1})
#        modtmp.add_reaction(rea1)
#
#        rea2 = Reaction("Refresh_NADH", lower_bound=0, upper_bound=1000)
#        rea2.add_metabolites({modtmp.metabolites.cpd00004_c0:-1, modtmp.metabolites.cpd00067_c0:1, modtmp.metabolites.cpd00003_c0:1})
#        modtmp.add_reaction(rea2)
#        
#        rea3 = Reaction("Refresh_NADPH", lower_bound=0, upper_bound=1000)
#        rea3.add_metabolites({modtmp.metabolites.cpd00005_c0:-1, modtmp.metabolites.cpd00067_c0:1, modtmp.metabolites.cpd00006_c0:1})
#        modtmp.add_reaction(rea3)
#
#        modtmp.add_boundary(modtmp.metabolites.cpd00012_c0, type='demand') # ppi
#        modtmp.add_boundary(modtmp.metabolites.cpd00009_c0, type='demand') # pi
#        modtmp.add_boundary(modtmp.metabolites.cpd00067_c0, type='demand') # h+

        modtmp.objective = modtmp.add_boundary(m, type='demand')
        obj = modtmp.slim_optimize()
        if round(obj,6) > 0:
            print obj, m.id, m.name
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
    print "Could not synthesis:", ",".join([m.id for m in fixLater+fail])
    print "Could not synthesis:", ",".join([m.name for m in fixLater+fail])
    if fill_gaps:
        return(modnew)

def checkReaction(model, reaction, direction=1):
    modtmp = set_medium(model, minmed)
    rea = modtmp.reactions.get_by_id(reaction)

    modtmp.objective = {rea:direction}
    sol = modtmp.optimize()
    if round(sol.objective_value ,12) > 0:
        print rea.id, rea.name, "can be used:", sol.objective_value
        print "\t Carbon source usage", carbon, sol.fluxes[carbon]
        return sol
    else:
        print rea.id, rea.name, "CANNOT be used:", sol.objective_value
        print "\t --> try to find reactions for gapfilling"
        try:
            gapsol = GapFiller(modtmp, refmod, demand_reactions=False).fill()
        except RuntimeError:
            print "\t => Runtime error: lowering the integer_threshold?"
            return
        except:
            print "\t => failed:", sys.exc_info()[0]
            return
        if len(gapsol[0]) > 0:
            print "\t => could be fixed:", ",".join([r.id for r in gapsol[0]])
            model.add_reactions([r for r in gapsol[0] if r not in model.reactions])
            print "\t Carbon source usage", carbon, sol.fluxes[carbon]


checkProduction(mod, "cpd00002_c0") # atp
#checkProduction(mod, "cpd00002_c0", addMets=True) # atp

sys.exit(0)
mod.add_boundary(mod.metabolites.cpd00229_c0, type="demand") # glycoaldehyde sink (~thf biosynthesis)
mod.add_reaction(refmod.reactions.rxn05039_c0) # riboflavin is not found by gapfiller
checkProduction(mod, "cpd02201_c0") # phosphopantothenate needed otherwise CoA biosynthesis does not work
checkProduction(mod, "cpd00264_c0", addMets=True) # spermidine


modnew = checkBiomass(mod, fill_gaps=True)
fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
write_sbml_model(modnew, filename=fileID+"_gapfilled.xml")


#for b in mod.reactions.bio1.reactants:
#    tmp = mod.copy()
#    rea = tmp.reactions.bio1
#    biom= tmp.metabolites.get_by_id(b.id)
#    rea.subtract_metabolites({biom:rea.metabolites[biom]})
#    if round(tmp.slim_optimize(),3) > 0:
#        print biom
#cpd = mod.metabolites.get_by_id("cpd00063_c0")
#mod.reactions.bio1.subtract_metabolites({cpd:mod.reactions.bio1.metabolites[cpd]})
