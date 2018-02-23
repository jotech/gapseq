from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from cobra import Model, Reaction, Metabolite
import pandas
import sys
import re
from cobra.flux_analysis import gapfill
from cobra.flux_analysis.gapfilling import GapFiller
import os
import seed_fix
import gapfill_reference


if len(sys.argv) != 4:
  print("usage:", sys.argv[0], "file.xml transporter.lst degradation.lst")
  sys.exit(0)


def add_Exchanges(mod, exids):
    model = mod.copy()
    for ex in [ex for ex in exids if ex not in model.reactions]:
        mid =re.findall("cpd[0-9]{5}", ex)[0]
        hit = seedmDB.loc[seedmDB["id"]==mid]
        m = Metabolite(id=mid+"_e0", name=hit["name"].values[0], formula=hit["formula"].values[0], charge=hit["charge"].values[0])
        r = Reaction(id=ex, name="Exchange for "+m.name, lower_bound=-1000, upper_bound=1000)
        r.add_metabolites({m:-1})
        model.add_reaction(r)
    return(model)

minmed = ["cpd00007", "cpd00001", # o2, h2o
          #"cpd00027", # glc
          #"cpd00020", # pyruvate
          #"cpd00106", # fumarat # myb10
          "cpd00013", # nh3
          #"cpd00048", # sulfate
          #"cpd11606", # menaquinone7" myb71
          #"cpd00305", # thiamine
          #"cpd01188",# lanosterol
          #"cpd00220", # riboflavin
          #"cpd00028", # heme myb11
          #"cpd00264", # spermidin myb10 # attention: spermidine is itself a carbon source!!
          # minerals
          "cpd00009", "cpd00030", "cpd00034", "cpd00058", "cpd00063", "cpd00067", "cpd00099", "cpd00149","cpd00205", "cpd00254", "cpd00971", "cpd01012", "cpd10515", "cpd10516", "cpd11574"]

def set_medium(model, medium, verbose=True):
    mod = model.copy()
    for ex in mod.exchanges:
        ex.lower_bound=0

    for m in medium: # set medium
        #ex = "EX_"+m+"_LPAREN_e_RPAREN_" # vmh
        if not m.startswith("EX_"):
            ex = "EX_"+m+"_e0"
        else:
            ex = m
        if ex in mod.reactions:
            mod.reactions.get_by_id(ex).lower_bound=-1000
        else:
            if verbose:
                print "compound not found in model:", m
    return(mod)

def set2respiration(model):
    mod = model.copy()
    mql = mod.metabolites.get_by_id("cpd15499_c0")
    mqn = mod.metabolites.get_by_id("cpd15500_c0")
    uql = mod.metabolites.get_by_id("cpd15561_c0")
    uqn = mod.metabolites.get_by_id("cpd15560_c0")
    h   = mod.metabolites.get_by_id("cpd00067_c0")
    esp1 = Reaction(id="ESP1", lower_bound=0, upper_bound=1000)
    esp2 = Reaction(id="ESP2", lower_bound=0, upper_bound=1000)
    esp1.add_metabolites({mql:-1,h:2,mqn:1})
    esp2.add_metabolites({uql:-1,h:2,uqn:1})
    mod.add_reactions([esp1,esp2])
    mod.objective = {esp1:1,esp2:1}
    return(mod)


mod = read_sbml_model(sys.argv[1])
fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
dir = os.path.dirname(__file__)
seedmDB = pandas.read_csv(dir+ "/../dat/seed_metabolites.tsv", sep="\t")
seedrDB = pandas.read_csv(dir+ "/../dat/seed_reactions.tsv", sep="\t")
substances=pandas.read_csv(dir+"/../dat/sub2pwy.csv", sep=",")

with open(sys.argv[2]) as f:
    lines = f.read()
ReaTransport = lines.rstrip("\n").split(" ") # remove linebreak at file end
with open(sys.argv[3]) as f:
    lines = f.read()
ReaDegradation = lines.rstrip("\n").split(" ") # remove linebreak at file end

newmod = seed_fix.seed_fix(mod)

# 1) Add diffusion reactions
ReaDiffusion = ["EX_cpd00011", "rxn05467", # co2
             "EX_cpd00007","rxn09031",  # o2
             "EX_cpd11640","rxn10542",  # h2
             "EX_cpd00001","rxn08686",  # h2o
             "EX_cpd03481","rxn13196"]  # etoh
refDiffusion = gapfill_reference.get_reference(mod, ReaDiffusion, delete_unbalanced=False, verbose=0)
Nrea = len(newmod.reactions)
newmod.add_reactions([r for r in refDiffusion.reactions if r not in mod.reactions])
print "Added diffusion reactions:", len(newmod.reactions)-Nrea

# 2) Add exchanges
ReaExchange = [ex for ex in ReaTransport if ex.startswith("EX_")]
Nrea = len(newmod.reactions)
newmod = add_Exchanges(newmod, ReaExchange)
print "Added exchange reactions:", len(newmod.reactions)-Nrea, "\n"

# 3) Get reference reactions for gapfilling
ReaRef = ReaDegradation + list(set(ReaTransport).difference(set(ReaExchange)))
refmod = gapfill_reference.get_reference(mod, ReaRef, delete_unbalanced=True, verbose=1)
refmod = seed_fix.seed_fix(refmod)

# 4) set model objecrive to respiration (easier minimal medium)
tmpmod = set2respiration(newmod)

# ASSERTION 1: check for growth without carbon source
sol = set_medium(tmpmod, minmed, verbose=False).optimize()
print "Negative growth control", sol.status, sol.objective_value
if sol.status == "optimal" and round(sol.objective_value,3) > 0:
    print "ATTENTION: Found growth without carbon source!"
    print set_medium(tmpmod, minmed, verbose=False).summary()
    for r, flux in sol.fluxes.iteritems():
        if abs(flux) >= 1000:
            print flux, r, tmpmod.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True)
    sys.exit()

# ASSERTION 2: Adding some dummy carbon source (glc,fum,pyr) => there should be growth
sol = set_medium(tmpmod, minmed+["cpd00027", "cpd00020", "cpd00106"], verbose=False).optimize()
print "Positive growth control", sol.status, sol.objective_value
if sol.status == "optimal" and round(sol.objective_value,3) < 0:
    print "No growth possible, check minimal medium?"
    sys.exit()


# 5) Try Gapfilling for each exchange reaction and carbon source
print "\nStarting gapfilling"
tmpmod = set_medium(tmpmod, minmed, verbose=False)
debugmod = tmpmod.copy()
debugmod.add_reactions([r for r in refmod.reactions if r not in debugmod.reactions])
write_sbml_model(debugmod, filename=fileID+"_debug.xml")
Nrea = len(tmpmod.reactions)
s1 = [ex for ex in substances["exid_seed"].dropna().values]
s2 = [ex.id for ex in tmpmod.exchanges]
csources = set(s1+s2)
tmpmod = add_Exchanges(tmpmod, set(s1).difference(set(s2))) # add exchange reactions
newmod = add_Exchanges(newmod, set(s1).difference(set(s2))) # add exchange reactions
Nfix = 0
for cs in csources:
    csname = tmpmod.reactions.get_by_id(cs).metabolites.keys()[0].name
    med = minmed+[cs]
    tmpmod = set_medium(tmpmod, med, verbose=False)
    sol = tmpmod.slim_optimize()
    if round(sol,12) > 0:
        pass
        print csname, "can be used:", sol
    else:
        print csname, "CANNOT be used:", sol
        print "\t --> try to find reactions for gapfilling"    
        try:
            #gapsol = GapFiller(tmpmod, refmod, demand_reactions=False, integer_threshold=1e-16).fill()
            gapsol = GapFiller(tmpmod, refmod, demand_reactions=False).fill()
        except RuntimeError:
            print "\t => Runtime error: lowering the integer_threshold?"
            newmod.remove_reactions(cs) # remove exchange reaction for compounds that cannot be used
            continue
        except:
            print "\t => failed:", sys.exc_info()[0]
            newmod.remove_reactions(cs) # remove exchange reaction for compounds that cannot be used
            continue
        if len(gapsol[0]) > 0:
            Nfix += 1            
            print "\t => could be fixed:", ",".join([r.id for r in gapsol[0]])
            newmod.add_reactions([r for r in gapsol[0] if r not in newmod.reactions])
print "Fixed growth for compounds:", Nfix, "by adding reactions:", len(newmod.reactions) - Nrea

write_sbml_model(newmod, filename=fileID+"_csfilled.xml")
