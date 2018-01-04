# needs: pathway-tools -python-local-only-non-strict -lisp

import pythoncyc
import sys
import re

meta = pythoncyc.select_organism('meta')
pythoncyc.sendQueryToPTools("(select-organism :org-id 'META)")



ofile = open("./meta_rea.csv", "w")
ofile.write("id" + "\t" + "name" + "\t" + "altname" + "\t" + "ec" + "\t" + "kegg" + "\t" + "pwy" + "\n")
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
    kegghit = re.findall("R[0-9]{5}",str(rea.dblinks))
    if kegghit != []:
        kegg = ",".join(kegghit)
    else:
        kegg = ""
    if rea.in_pathway != None:
        pwy = ",".join([p.replace("|","") for p in rea.in_pathway])
    ofile.write(r + "\t" + name + "\t" + altname + "\t" + ec + "\t" + kegg + "\t" + pwy + "\n")
ofile.close()

