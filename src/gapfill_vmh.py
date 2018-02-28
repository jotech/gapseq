from cobra.flux_analysis import gapfill
from cobra.flux_analysis.gapfilling import GapFiller
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
import pandas
import sys
import re
import os
import gapfill_reference


if len(sys.argv) != 3:
  print("usage:", sys.argv[0], "file.xml referenceReactions.lst/referenceModel.xml")
  sys.exit(0)

print " ".join(sys.argv)

#gapfill_reference.get_reference_vmh(mod, newR)

dir = os.path.dirname(__file__)
substances=pandas.read_csv(dir+"/../dat/sub2pwy.csv", sep=",")

mod = read_sbml_model(sys.argv[1])

if sys.argv[2].endswith(".lst"):
    with open(sys.argv[2]) as f:
        lines = f.read()
    newR = lines.rstrip("\n").split(" ") # remove linebreak at file end
    # Create reference model with reactions pool
    refmod = gapfill_reference.get_reference_vmh(mod, newR, verbose=False)
else:
    refmod = read_sbml_model(sys.argv[2])
    for m in refmod.metabolites: # compartment tag inconsistency
        m.id = re.sub("_(.)$","[\\1]", m.id)
    for ex in refmod.exchanges:
        ex.id = ex.id.replace("_LPAREN_","(").replace("_RPAREN_",")").replace("_DASH_","_") # compartment tag for exchange reactions


ignore = ["H", "Li", "Na", "Mg", "K", "Ca", "P", "Fe", "Cl", "Zn", "Mn", "Mo", "Se", "Co", "Cu", "Ni", "W", "H2O"]
# prepare model with minimal medium
minmed = ["EX_o2(e)", "EX_h2o(e)", # o2, h2o
          "EX_nh4(e)", # nh4
          # minerals
          "EX_ca2(e)","EX_pi(e)","EX_mg2(e)","EX_na1(e)","EX_k(e)","EX_cl(e)","EX_fe2(e)","EX_fe3(e)","EX_mn2(e)","EX_zn2(e)","EX_cu2(e)"]

# try to find a carbon source
cs = ["EX_glc(e)", # glucose
      "EX_pyr(e)", # pyruvate
      "EX_pro_L(e)",
      "EX_fum(e)", # fumarat
      "EX_mal(e)"] # maltose
found = [c for c in cs if c in [ex.id for ex in mod.exchanges]]
if len(found) > 0:
    carbon = found[0]
    print "Try doing biomass gapfilling with carbon source", carbon, mod.reactions.get_by_id(carbon).name
    minmed += [carbon]
else:
    found2= substances["exid"].loc[substances["exid"].isin([ex.id for ex in mod.exchanges])].tolist()
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
            ex = "EX_"+m+"_LPAREN_e_RPAREN_" # vmh
        else:
            ex = m
        if ex in mod.reactions:
            mod.reactions.get_by_id(ex).lower_bound=-1000
        else:
            print "compound not found in model:", m
    return(mod)

def checkProduction(model, metabolite):
    #modtmp = set_medium(model, minmed)
    modtmp = model.copy()
    #print modtmp.medium
    #print refmod.metabolites.get_by_id("atp[c]").reactions
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
    #modmin = model.copy()
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
    biomass = [r.id for r in model.reactions if r.objective_coefficient!=0][0]
    allMet = model.reactions.get_by_id(biomass).reactants 
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


mod.add_reactions([ex for ex in refmod.exchanges if ex not in mod.reactions])
checkProduction(mod, "atp[c]")
sol = mod.optimize()
if sol.status != "optimal" or round(sol.objective_value,6) == 0:
    print "no biomass production possible try to fix"
    tmpmod = mod.copy()
    for ex in tmpmod.exchanges:
        ex.lower_bound=-1000
    gapsol = GapFiller(tmpmod, refmod, demand_reactions=False, integer_threshold=1e-16).fill()
    print gapsol

modnew = checkBiomass(mod, fill_gaps=True)
fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
write_sbml_model(modnew, filename=fileID+"_gapfilled.xml")
