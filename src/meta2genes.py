import pythoncyc
import pandas

genes   = []
rxn     = []
ncbi    = []
pwy     = []
uniprot = []
location= []
go      = []
org     = []

for r in meta.all_reactions():
    #r="RXN-12309"
    for g in meta.genes_of_reaction(r):
        gene = meta[g]
        genes.append(gene.common_name)
        rxn.append(r)
        pwy.append(";".join(meta.pathways_of_gene(g)))
        
        try:
            ncbi.append(gene.dblinks["|ENTREZ|"][0])
        except:
            ncbi.append("")
        try:
            uniprot.append(";".join([meta[p].dblinks['|UNIPROT|'][0] for p in gene.product]))
        except:
            uniprot.append("")
        try:
            location.append(";".join([";".join(meta[p].locations) for p in gene.product]))
        except:
            location.append("")
        try:
            go.append(";".join([";".join(meta[p].go_terms) for p in gene.product]))
        except:
            go.append("")
        try:
            org.append(";".join([";".join(meta[p].species) for p in gene.product]))
        except:
            org.append("")


df = pandas.DataFrame()
df["rxn"] = rxn
df["genes"] = genes
df["pwy"] = pwy
df["ncbi"] = ncbi
df["uniprot"] = uniprot
df["location"] = location
df["go"] = go
df["org"] = org

df.to_csv("meta_genes.csv", index=False)
