def correct_mass(dic): # how to handle wrong charge????
    if dic == {}:
        return True
    else:
        return False

def only_missing_H(dic):
    if len(dic)==1 and "H" in dic:
        return True
    else:
        return False

def repair_mass_balance(model, delete_unbalanced=True, verbose=1):
    from cobra import Model, Reaction, Metabolite
    import re
    import pandas
    import os
    dir = os.path.dirname(__file__)
    seedmDB = pandas.read_csv(dir+"/../dat/seed_metabolites.tsv", sep="\t")

    mod = model.copy()
    Cwrong = 0
    Ccorrected = 0
    if "cpd00067_c0" not in mod.metabolites:
        h = Metabolite(id="cpd00067_c0", name="H+", formula="H", charge=1)
        mod.add_metabolites(h)
    
    rea_rm = []
    
    mets = [m for m in mod.metabolites if m.formula==""]
    for m in mets:
        mid =re.sub("_.*$","", m.id)
        hit = seedmDB.loc[seedmDB["id"]==mid]
        m.formula = hit["formula"].values[0]
        m.charge = hit["charge"].values[0]
    print len(mets), "compounds without formula, fixed from total", len(mod.metabolites) 
    
    for r in mod.reactions:
        if r in mod.exchanges or r.id=="bio1":
            continue
        
        mbal = r.check_mass_balance()
        if "charge" in mbal:
            mbal.pop("charge") # ToDo: how to handle charge? (is not always in agreement with H!)
        if not correct_mass(mbal): 
            if verbose == 2:
                print r.id, r.reaction
            Cwrong += 1
            if not only_missing_H(mbal):
                for m in r.metabolites: # check if one substance causes miss balance => remove
                    tmp = r.copy()
                    coef = tmp.get_coefficient(m.id)
                    #print m.id, coef
                    tmp.subtract_metabolites({m:coef})
                    mbal2 = tmp.check_mass_balance()
                    if "charge" in mbal2:
                        mbal2.pop("charge")
                    if only_missing_H(mbal2):
                        if verbose == 2:
                            print "removed", m.id,m.name, "from", r.id
                            print mbal, mbal2
                        mbal = mbal2
                        r.subtract_metabolites({m:-1})
                        break
            if only_missing_H(mbal): # check again if corrections could work
                miss = mbal["H"]
            else:
                if verbose == 2:
                    print "could not be fixed"
                    print "\t",r.id, r.name, r.reaction
                    print "\t",mbal
                if delete_unbalanced:
                    rea_rm.append(r)
                    #mod.remove_reactions(r)
                continue
            if miss != round(miss): # check if integer
                print "error in", r.id, "mass imbalance:", mbal
                if delete_unbalanced:
                    rea_rm.append(r)
                    #mod.remove_reactions(r)
                continue
            h = mod.metabolites.get_by_id("cpd00067_c0")
            if miss < 0:
                r.add_metabolites({h: -miss})
            if miss > 0:
                r.subtract_metabolites({h: miss})
            #print mbal, r.check_mass_balance()
            mbal3 = r.check_mass_balance()
            if "charge" in mbal3:
                mbal3.pop("charge")
            if not correct_mass(mbal3):
                if verbose == 2:
                    print r.id, r.name, r.reaction
                    print "old:",mbal, "\tnew:",mbal3
                if delete_unbalanced:
                    rea_rm.append(r)
                    #mod.remove_reactions(r)
            else:
                if verbose == 2:
                    print r.id, "fixed"
                Ccorrected += 1
    mod.remove_reactions(rea_rm)
    #print ",".join([r.name for r in rea_rm])
    #print [r.check_mass_balance() for r in rea_rm]
    mod.repair()
    if verbose > 0:
        print "Unbalanced reactions:", Cwrong, "\tCould be corrected:", Ccorrected, "\tRemoved reactions:", len(model.reactions)-len(mod.reactions)
    return(mod)

def get_reference(mod, newR, delete_unbalanced, verbose=1):
    from cobra import Model, Reaction, Metabolite
    import pandas
    import re
    import os
    dir = os.path.dirname(__file__)
    seedmDB = pandas.read_csv(dir+"/../dat/seed_metabolites.tsv", sep="\t")
    seedrDB = pandas.read_csv(dir+"/../dat/seed_reactions.tsv", sep="\t")

    #refdb  = seedrDB.loc[seedrDB['abbreviation'].isin(newR)] # vmh
    refdb  = seedrDB.loc[(seedrDB['id'].isin(newR)) & (seedrDB['status']=="OK")] # only choose reaction which have status okay
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
        rstr = rstr.replace("[0]", "_c0").replace("[1]", "_e0").replace("[3]", "_p0").replace("[2]", "_m0")
        r.build_reaction_from_string(rstr, verbose=0)
        #r.reaction = rstr
    for m in refmod.metabolites:
        if verbose == 2:
            print m.id
        mid =re.sub("_.*$","", m.id)
        hit = seedmDB.loc[seedmDB["id"]==mid]
        m.name = hit["name"].values[0]
        m.formula = hit["formula"].values[0]
        m.charge = hit["charge"].values[0]
    
    if verbose > 0:
        print Calready, "reactions already in the model"
        print Cold, "removed deprecated reactions"
    refmod = repair_mass_balance(refmod, delete_unbalanced, verbose)
    if verbose > 0:
        print len(refmod.reactions), "remaining reaction in reference database:", 
    return(refmod)
