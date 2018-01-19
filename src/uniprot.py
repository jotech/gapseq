from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from bioservices import UniProt
import pandas
from io import StringIO
import sys
import os

# TODO: do not search for incomplete EC numbers
# TODO: should a limit be set? u.search(..., limit=10) ~ blast speed

dir_path = os.path.dirname(os.path.realpath(__file__))

if len(sys.argv) != 2:
  print("usage:", sys.argv[0], Pathway/subsystem)
  print("Pathway/subsystem: ID from metacyc hierarchy")
  sys.exit(0)

print "Looking for pathways: ", sys.argv[1]
pwydf_src = pandas.read_csv(dir_path+"/../dat/meta_pwy.tbl", sep="\t")
pwydf = pwydf_src[pwydf_src["hierarchy"].str.contains(sys.argv[1])]
if len(pwydf.index) == 0:
    pwydf = pwydf_src[pwydf_src["id"].str.contains(sys.argv[1])]
print "\tfound:", len(pwydf.index), "\n"
if len(pwydf.index) == 0:
    sys.exit(0)

u = UniProt()
for index,row in pwydf.iterrows():
    print row["id"], row["name"]
    if row["reaEc"]== "" or row["reaEc"] == None or pandas.isnull(row["reaEc"]):
        continue
    ec_list = row["reaEc"].split(",")
    for ec in ec_list:
        results = u.search(ec+"+and+reviewed:yes", columns="id,entry name, protein names, sequence", limit=20) # uniprot swissprot db (reviewed)
        if len(results) == 0:
            results = u.search(ec, columns="id,entry name, protein names, sequence", limit=20) # uniprot trembl db (unreviewed) limit=10
        if len(results) == 0:
            print("No entry found for:", row["name"],row["id"], ec)
            continue
        df = pandas.read_csv(StringIO(results), sep="\t")
        records = []
        for index, row2 in df.iterrows():
            record = SeqRecord(Seq(row2['Sequence'], IUPAC.protein), id=row2['Entry'], description=row2['Protein names'])
            records.append(record)
        SeqIO.write(records, dir_path+"/../dat/seq/"+ec+".fasta", "fasta")
