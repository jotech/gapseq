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


if len(sys.argv) != 3:
  print("usage:", sys.argv[0], "file.xml transporter.lst")
  sys.exit(0)

dir = os.path.dirname(__file__)
seedmDB = pandas.read_csv(dir+"/../dat/seed_metabolites.tsv", sep="\t")
seedrDB = pandas.read_csv(dir+"/../dat/seed_reactions.tsv", sep="\t")
substances = pandas.read_csv(dir+"/../dat/sub2pwy.csv", sep=",")

mod = read_sbml_model(sys.argv[1])
with open(sys.argv[2]) as f:
    lines = f.read()
newR = lines.rstrip("\n").split(" ") # remove linebreak at file end

# diffusion reactions should be added always!
diffusion = ["EX_cpd00011", "rxn05467", # co2
             "EX_cpd00007","rxn09031",  # o2
             "EX_cpd11640","rxn10542",  # h2
             "EX_cpd00001","rxn08686",  # h2o
             "EX_cpd03481","rxn13196"]  # etoh
newR = newR + diffusion

def get_reference(mod, newR, delete_unbalanced, verbose):
    #refdb  = seedrDB.loc[seedrDB['abbreviation'].isin(newR)] # vmh
    refdb  = seedrDB.loc[seedrDB['id'].isin(newR)]
    refmod = Model("reaction_database")
    Rmod   = [re.sub("_.*$", "", r.id) for r in mod.reactions]
    print "Consider", len(refdb.index), "reactions"
    Cold = 0
    Calready = 0
    for i,row in refdb.iterrows():
        #rid = row["abbreviation"] # vmh
        rid = row["id"]
        if row["is_obsolete"] == 1: # do not add old reactions
            Cold += 1
            continue
        #if rid in Rmod: # vmh
        elif rid in Rmod:
            Calready += 1
            continue
        r = Reaction(rid+"_c0") # default seed model have compartment tag
        refmod.add_reaction(r)
        #rstr = row["formula"] # vmh
        rstr = row["equation"]
        #rstr = row["code"]
        rstr = rstr.replace("[0]", "_c0").replace("[1]", "_e0").replace("[3]", "_p0")
        if "[2]" in rstr: # TODO: consider all compartments!
            continue
        r.build_reaction_from_string(rstr, verbose=False)
        #r.reaction = rstr
    for m in refmod.metabolites:
        mid =re.sub("_.*$","", m.id)
        hit = seedmDB.loc[seedmDB["id"]==mid]
        m.name = hit["name"].values[0]
        m.formula = hit["formula"].values[0]
        m.charge = hit["charge"].values[0]

    print Calready, "reactions already in the model"
    print Cold, "removed deprecated reactions"
    #refmod = repair_mass_balance(refmod, delete_unbalanced, verbose)
    print len(refmod.reactions), "remaining reaction in reference database:", 
    return(refmod)

refmod = get_reference(mod, newR, delete_unbalanced=True, verbose=False)

def add_Exchanges(mod, exids):
    model = mod.copy()
    for ex in exids:
        mid =re.findall("cpd[0-9]{5}", ex)[0]
        hit = seedmDB.loc[seedmDB["id"]==mid]
        m = Metabolite(id=mid+"_e0", name=hit["name"].values[0], formula=hit["formula"].values[0], charge=hit["charge"].values[0])
        r = Reaction(id=ex, name="Exchange for "+m.name, lower_bound=-1000, upper_bound=1000)
        r.add_metabolites({m:-1})
        model.add_reaction(r)
    return(model)

refmod = add_Exchanges(refmod, [ex for ex in newR if ex.startswith("EX_")])

newmod = mod.copy()
reallynew = [r for r in refmod.reactions if r not in mod.reactions]
newmod.add_reactions(reallynew)
newmod = seed_fix.seed_fix(newmod)
print "added: ", len(newmod.reactions) - len(mod.reactions), "reactions"
print "already in model:", ",".join(set(newR).difference(set([r.id.replace("_c0","") for r in reallynew])))


fileID = os.path.splitext(os.path.basename(sys.argv[1]))[0]
write_sbml_model(newmod, filename=fileID+"_transporter.xml")

