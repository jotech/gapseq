def addFromSeedDB(mod, lineFromSeedDB):
    import pandas
    from cobra import Model, Reaction, Metabolite
    comp = {'0':'_c0', '1':'_e0', '2':'_p0'}
    model = mod.copy()
    formula = lineFromSeedDB["stoichiometry"].values[0]
    rea_met = {m.split(":")[1]+comp[m.split(":")[2]]:float(m.split(":")[0]) for m in formula.split(";")}
    rea_id  = lineFromSeedDB["id"].values[0] + "_c0" # always the case??
    rea_name= lineFromSeedDB["name"].values[0]
    rea = Reaction(id=rea_id, name=rea_name)
    for m in rea_met:
        if m not in model.metabolites:
            model.add_metabolites(Metabolite(m))
    model.add_reactions(rea)
    rea_stoich = {model.metabolites.get_by_id(m):rea_met[m] for m in rea_met}
    rea.add_metabolites(rea_stoich)
    return(model)

def getReferenceMod(ReaList):
    from pandas import read_csv
    from cobra import Model, Reaction, Metabolite
    import os
    dir = os.path.dirname(__file__)
    seed_rea = read_csv(dir+"/../dat/seed_reactions_corrected.tsv", sep="\t")
    refmod = Model("Reference model")
    comp = {'0':'_c0', '1':'_e0', '2':'_p0'}
    for r in ReaList: 
        line = seed_rea[seed_rea["id"]==r]
        if len(line.index) == 0:
            continue
        
        # check curration status of reaction
        status = line["gapseq.status"].values[0]
        if status == "not.assessed" or status == "removed":
            continue 

        # stoichiometry
        formula = line["stoichiometry"].values[0]
        rea_met = {m.split(":")[1]+comp[m.split(":")[2]]:float(m.split(":")[0]) for m in formula.split(";")}
        rea_id  = line["id"].values[0] + "_c0" # always the case??
        rea_name= line["name"].values[0]
        rea = Reaction(id=rea_id, name=rea_name)
        for m in rea_met:
            if m not in refmod.metabolites:
                refmod.add_metabolites(Metabolite(m))
        if rea not in refmod.reactions:
            refmod.add_reactions(rea)
        rea_stoich = {refmod.metabolites.get_by_id(m):rea_met[m] for m in rea_met}
        rea.add_metabolites(rea_stoich)
        
        # reversibility
        reversibility = line["reversibility"].values[0]
        if reversibility=='=':
            lb =-1000
            ub = 1000
        elif reversibility=='>':
            lb = 0
            ub = 1000
        elif reversibility=='<':
            lb =-1000
            ub =  0
        rea.lower_bound = lb
        rea.upper_bound = ub

    print "Reference model for gapfilling with", len(refmod.reactions), "reactions and", len(refmod.metabolites), "metabolites"
    return(refmod)
