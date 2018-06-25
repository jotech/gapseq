# needs: pathway-tools -python-local-only-non-strict -lisp

import pythoncyc
import sys

meta = pythoncyc.select_organism('meta')
pythoncyc.sendQueryToPTools("(select-organism :org-id 'META)")

def getReaInfos(pwy):
    if meta[pwy] == None:
        print(pwy, "Pathway does not exist")
        return([[],[]])
    rea_list = meta[pwy]["reaction_list"]

    superpathway = filter(lambda x:'|PWY' in x,rea_list)
    if len(superpathway) > 0:
        rea_list_new = []
        for sp in superpathway:
            rea_list_new.extend(meta[sp]["reaction_list"])
        rea_list = rea_list_new

    ec_list  = []
    reaid_list = []
    for r in rea_list:
        if meta[r]["spontaneous_p"]:
            continue
        ec  = meta[r]["ec_number"]
        reaid = r
        if ec == None:
            # there are reaction without EC number in metacyc which have a EC number assigned in the pathway-overview
            rea_related = filter(lambda x:'RXN' in x,meta[r]["in_pathway"])
            if len(rea_related) == 0:
                continue
            if len(rea_related) > 1:
                # if there are more than 1 related reaction, it is unclear which one to choose
                print("unclear ecs substitution:", r, pwy)
                continue
            r2 = rea_related[0]
            ec = meta[r2]["ec_number"]
            if ec == None:
                continue
        ec_nr=len(ec)
        ec = str(",".join(ec)).replace("|","").replace("EC-","")
        ec_list.append(ec)
        reaid = reaid.replace("|","")
        reaid_list.append(",".join([reaid for i in range(0,ec_nr)]))
    ec_col = (",".join(ec_list))
    rea_col= (",".join(reaid_list)) 
    return([ec_col,rea_col])


ofile = open("./meta_pwy.tbl", "w")
ofile.write("id" + "\t" + "name" + "\t" + "altname" + "\t" + "hierarchy" + "\t" + "taxrange" + "\t" + "reaId" + "\t" + "reaEc" + "\t" + "keyRea" + "\n")
for p in meta.all_pathways():
    pwy = meta[p]
    qry = "(get-instance-all-types '"+p+")"
    hierarchy = ",".join(pythoncyc.sendQueryToPTools(qry))
    name   = pwy.common_name.replace("|","")
    altname= ",".join(pwy.names)
    reaInf = getReaInfos(p)
    reaEc = reaInf[0]
    reaId = reaInf[1]
    # TODO: causes problems somehow?
    #taxrange = ",".join([meta[t].common_name.replace("|","") for t in pwy.taxonomic_range])
    if pwy.taxonomic_range != None:
        taxrange = ",".join(pwy.taxonomic_range)
    else:
        taxrange = ""
    if pwy.key_reactions != None:
        keyRea = ",".join(pwy.key_reactions).replace("|","")
    else:
        keyRea = ""
    ofile.write(p + "\t" + name + "\t" + altname + "\t" + hierarchy + "\t" + taxrange + "\t" + reaId + "\t" + reaEc + "\t" + keyRea  +"\n")
ofile.close()
# sort file!!
# sort -k 1 dat/meta_pwy.tbl
