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
