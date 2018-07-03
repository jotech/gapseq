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
    
    keyRea_list = []
    if meta[pwy].key_reactions != None: # get key enzymes
        keyRea_list = meta[pwy].key_reactions
    
    isSuperPwy = False
    superpathway = filter(lambda x:'PWY' in x,rea_list)
    if len(superpathway) > 0:
        isSuperPwy = True
        rea_list_new = []
        for sp in superpathway:
            rea_list_new.extend(meta[sp]["reaction_list"])
            if meta[sp].key_reactions != None: # get key enzymes
                keyRea_list.extend(meta[sp].key_reactions)
        rea_list = rea_list_new
    keyRea = ",".join(keyRea_list).replace("|","")
    #print(keyRea)

    ec_list  = []
    reaId_list = []
    reaName_list = []
    rea_counter = 0 # count non-spontanous reactions
    rea_spont_counter = 0
    
    for r in rea_list:
        if meta[r]["spontaneous_p"]:
            rea_spont_counter +=1
            continue
        rea_counter +=1 
        ec  = meta[r]["ec_number"]
        reaId = r
        if meta[r].common_name != None:
            reaName = meta[r].common_name
        else:
            set1 = set(meta.enzymes_of_pathway(pwy))
            set2 = set(meta.enzymes_of_reaction(r))
            enzyme = {meta.enzyme_activity_name(e, pwy) for e in set1.intersection(set2)} 
            if len(enzyme) > 1:
                print("more than one enzyme-reaction match:",pwy,r)
                #sys.exit(0)
                reaName = enzyme.pop()
            if len(enzyme) == 1:
                reaName = "".join(enzyme)
            if len(enzyme) == 0:
                reaName = reaId
            #if meta[r].enzymatic_reaction != None:
            #    if meta[meta[r].enzymatic_reaction[0]].common_name != None:
            #        reaName = meta[meta[r].enzymatic_reaction[0]].common_name
            #    else:
            #        reaName = reaId
            #else:
            #    reaName = reaId
        if ec == None:
            # there are reaction without EC number in metacyc which have a EC number assigned in the pathway-overview
            rea_related = filter(lambda x:'RXN' in x,meta[r]["in_pathway"])
            if len(rea_related) == 0:
                ec = [""]
            elif len(rea_related) > 1:
                # if there are more than 1 related reaction, it is unclear which one to choose
                print("unclear ecs substitution:", r, pwy)
                ec = [""]
            else:
                r2 = rea_related[0]
                ec = meta[r2]["ec_number"]
                if ec == None:
                    ec = [""]
        ec_nr=len(ec)
        ec = str(",".join(ec)).replace("|","").replace("EC-","") # more than one EC possible
        ec_list.append(ec)
        reaId = reaId.replace("|","")
        reaId_list.append(",".join([reaId for i in range(0,ec_nr)]))
        #reaName = reaName.replace("|","")
        reaName_list.append(";".join([reaName for i in range(0,ec_nr)]))
    ec_col = (",".join(ec_list))
    rea_col= (",".join(reaId_list)) 
    reaName_col= (";".join(reaName_list)) 
    status = rea_counter == len(ec_list) # check if all reactions are covered by EC numbers
#    print("ec_list:", ec_list)
#    print("reaName_list:", reaName_list)
#    print(rea_counter, len(ec_list))
    ec_counter = len([ec for ec in ec_list if ec!=""])
    return([ec_col, rea_col, status, rea_counter, isSuperPwy, keyRea, reaName_col, ec_counter])

#getReaInfos("|PWY-6381|")
#getReaInfos("|p381-PWY|")
#getReaInfos("|BRANCHED-CHAIN-AA-SYN-PWY|")
#sys.exit(0)


ofile = open("./meta_pwy.tbl", "w")
ofile.write("id" + "\t" + "name" + "\t" + "altname" + "\t" + "hierarchy" + "\t" + "taxrange" + "\t" + "reaId" + "\t" + "reaEc" + "\t" + "keyRea" + "\t"+ "reaName" + "\t" + "reaNr" + "\t" + "ecNr" + "\t" + "superpathway" + "\t" "status" + "\n")
for p in meta.all_pathways():
    pwy = meta[p]
    qry = "(get-instance-all-types '"+p+")"
    hierarchy = ",".join(pythoncyc.sendQueryToPTools(qry))
    name    = pwy.common_name.replace("|","")
    altname = ",".join(pwy.names)
    reaInf  = getReaInfos(p)
    reaEc   = reaInf[0]
    reaId   = reaInf[1]
    status  = reaInf[2]
    reaN    = reaInf[3]
    superPwy=reaInf[4]
    keyRea  = reaInf[5]
    reaName = reaInf[6]
    ecN     = reaInf[7]
    # TODO: causes problems somehow?
    #taxrange = ",".join([meta[t].common_name.replace("|","") for t in pwy.taxonomic_range])
    if pwy.taxonomic_range != None:
        taxrange = ",".join(pwy.taxonomic_range)
    else:
        taxrange = ""
    ofile.write(p + "\t" + name + "\t" + altname + "\t" + hierarchy + "\t" + taxrange + "\t" + reaId + "\t" + reaEc + "\t" + keyRea + "\t" + reaName + "\t" + str(reaN) + "\t" + str(ecN) + "\t" + str(superPwy) +"\t" + str(status)  + "\n")
ofile.close()
# sort file!!
# sort -k 1 dat/meta_pwy.tbl
