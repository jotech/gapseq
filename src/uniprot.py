from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from bioservices import UniProt
import pandas
from io import StringIO
import sys
import os
import re

# TODO: do not search for incomplete EC numbers
# TODO: should a limit be set? u.search(..., limit=10) ~ blast speed

dir_path = os.path.dirname(os.path.realpath(__file__))

# usage
if len(sys.argv) != 2:
  print("usage:", sys.argv[0], Pathway/subsystem/EC)
  print("Pathway/subsystem/EC: ID from metacyc hierarchy, pathway ID or EC number")
  sys.exit(0)


def download_EC(ec):
    results = u.search(ec+"+and+reviewed:yes", columns="id,entry name, protein names, sequence", limit=20) # uniprot swissprot db (reviewed)
    if len(results) == 0:
        results = u.search(ec, columns="id,entry name, protein names, sequence", limit=20) # uniprot trembl db (unreviewed) limit=10
    if len(results) == 0:
        print "\t\tNo entry found for:", ec
        return
    print "\t\t downloading ... ", ec, "..."
    df = pandas.read_csv(StringIO(results), sep="\t")
    records = []
    for index, row2 in df.iterrows():
        record = SeqRecord(Seq(row2['Sequence'], IUPAC.protein), id=row2['Entry'], description=row2['Protein names'])
        records.append(record)
    SeqIO.write(records, dir_path+"/../dat/seq/"+ec+".fasta", "fasta")



if sys.argv[1].startswith("|"):
    key = sys.argv[1].replace("|","")
else:
    key = sys.argv[1]

print "\t --> Trying to download sequencing data for: ", key

hit = re.match("\d{1,2}(\.\d{1,2}){3}", key)
if hit:
    ec = hit.group()
    u = UniProt()
    download_EC(ec)
    sys.exit(0)

pwydf_src = pandas.read_csv(dir_path+"/../dat/meta_pwy.tbl", sep="\t")
pwydf = pwydf_src[pwydf_src["hierarchy"].str.contains(key)]
if len(pwydf.index) == 0:
    pwydf = pwydf_src[pwydf_src["id"] == "|"+key+"|"]
print "\tfound:", len(pwydf.index), "\n"
if len(pwydf.index) == 0:
    sys.exit(0)

u = UniProt()
for index,row in pwydf.iterrows():
    print "\t", row["id"], row["name"]
    if row["reaEc"]== "" or row["reaEc"] == None or pandas.isnull(row["reaEc"]):
        continue
    ec_list = row["reaEc"].split(",")
    for ec in ec_list:
        download_EC(ec)
