# needs: pathway-tools -python-local-only-non-strict -lisp

import pythoncyc
import sys
import re

meta = pythoncyc.select_organism('meta')
pythoncyc.sendQueryToPTools("(select-organism :org-id 'META)")



ofile = open("./meta_rea.tbl", "w")
ofile.write("id" + "\t" + "name" + "\t" + "altname" + "\t" + "ec" + "\t" + "kegg" + "\t" + "pwy" + "\t" + "enzymatic" + "\t" + "gibbs" + "\t" + "direction" + "\n")
for r in meta.all_reactions():
    rea = meta[r]
    #qry = "(get-instance-all-types '"+p+")"
    #hierarchy = ",".join(pythoncyc.sendQueryToPTools(qry))
    if rea.common_name != None:
        name   = rea.common_name
    else:
        name = ""
    if rea.names != None:
        altname= ",".join(rea.names)
    else:
        altname = ""
    if rea.ec_number != None:
        ec = ",".join([ec.replace("|","") for ec in rea.ec_number])
    else:
        ec = ""
    kegghit = re.findall("R[0-9]{5}",str(rea.dblinks))
    if kegghit != []:
        kegg = ",".join(kegghit)
    else:
        kegg = ""
    if rea.in_pathway != None:
        pwy = ",".join([p.replace("|","") for p in rea.in_pathway])
    else:
        pwy = ""
    if rea.enzymatic_reaction != None and len(rea.enzymatic_reaction) > 0:
        enzymatic = True
    else:
        enzymatic = False
    if rea.gibbs_0 != None:
        gibbs = rea.gibbs_0
    else:
        gibbs = ""
    if rea.reaction_direction != None:
        dir_dic = {"|PHYSIOL-LEFT-TO-RIGHT|":">", "|PHYSIOL-RIGHT-TO-LEFT|":"<", "|LEFT-TO-RIGHT|":">", "|RIGHT-TO-LEFT|":"<", "|REVERSIBLE|":"=", "|IRREVERSIBLE-RIGHT-TO-LEFT|":"<", "|IRREVERSIBLE-LEFT-TO-RIGHT|":">"}
        direction = dir_dic[rea.reaction_direction]
    else:
        direction = "?"
    ofile.write(r + "\t" + name + "\t" + altname + "\t" + ec + "\t" + kegg + "\t" + pwy + "\t" + str(enzymatic) + "\t" + str(gibbs) + "\t" + direction + "\n")
ofile.close()

