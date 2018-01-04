# needs: pathway-tools -python-local-only-non-strict -lisp

import pythoncyc

meta = pythoncyc.select_organism('meta')
pythoncyc.sendQueryToPTools("(select-organism :org-id 'META)")

#query = "|Amino-Acid-Biosynthesis|"
#query = "|Nucleotide-Biosynthesis|"
#query = "|Cofactor-Biosynthesis|"
#query = "|Carbohydrates-Degradation|"
query = "|Polyamine-Biosynthesis|"

instances = meta.get_class_all_instances(query)


ofile = open("/home/jo/pwy.csv", "w")
ofile.write("id" + "\t" + "name" + "\t" + "ec" + "\t" + "rea" + "\n")
for idx, pwy in enumerate(instances):
    if meta[pwy] == None:
        print(pwy, "Pathway does not exist")
        continue
    rea_list = meta[pwy]["reaction_list"]
    ec_list  = []
    reaid_list = []
    for r in rea_list:
        ec  = meta[r]["ec_number"]
        reaid = r
        if ec == None:
            continue
        ec_nr=len(ec)
        ec = str(",".join(ec)).replace("|","").replace("EC-","")
        ec_list.append(ec)
        reaid = reaid.replace("|","")
        reaid_list.append(",".join([reaid for i in range(0,ec_nr)]))
    ec_col = (",".join(ec_list))
    rea_col= (",".join(reaid_list)) 
    name   = meta[pwy].common_name.replace("|","")
    ofile.write(pwy + "\t" + name + "\t" + ec_col + "\t" + rea_col + "\n")
ofile.close()
