from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import pairwise2
from bioservices import UniProt
import pandas
from io import StringIO
import sys
import os
import re


# TODO: cluster similar sequences (seq identity implemented but not working. global/local alignment, cutoff??)
# TODO: should a limit be set? u.search(..., limit=10) ~ blast speed

dir_path = os.path.dirname(os.path.realpath(__file__))

# usage
if len(sys.argv) < 3:
  print("usage:", sys.argv[0], "Pathway/subsystem/EC, taxonomy [max_download]")
  print("Pathway/subsystem/EC: ID from metacyc hierarchy, pathway ID or EC number")
  print("Taxonomy: taxonomic range (default: Bacteria)")
  print("max_download: How many sequences should be download (default: 20)")
  sys.exit(0)

if len(sys.argv) == 4:
    max_download = sys.argv[3]
else: 
    max_download = 20

def identity(str1,str2):
    if not len(str1) == len(str2):
        print "WARNING: string length mismatch"
        return 0
    idC = 0
    l = len(str1)
    for i in range(l):
        if str1[i] == str2[i]:
            idC += 1
    return(idC/float(l))

def mean_identity(records, newseq):
    identities = []
    if len(records) == 0:
        return 0
    else:
        for r in records:
            ali = pairwise2.align.globalxx(r.seq, newseq.seq)[0]
            identities.append(identity(ali[0],ali[1]))
    Msum = 0
    for n in identities:
        Msum += n
    return(Msum/len(identities))


def download_EC(ec):
    results = u.search("ec:"+ec+"+and+reviewed:yes"+"+and+taxonomy:"+taxonomy, columns="id,entry name, protein names, sequence", limit=max_download) # uniprot swissprot db (reviewed)
    Nhits = len(results.split("\n")) - 2 # minus header & empty last line
    if Nhits < 15: # if only some sequences found get more
        results = u.search("ec:"+ec+"+and+taxonomy:"+taxonomy, columns="id,entry name, protein names, sequence", limit=max_download) # uniprot trembl db (unreviewed) limit=10
    if len(results) == 0:
        print "\t\tNo entry found for:", ec
        return
    print "\t\t downloading ... ", ec, "..."
    df = pandas.read_csv(StringIO(results), sep="\t")
    records = []
    for index, row2 in df.iterrows():
        record = SeqRecord(Seq(row2['Sequence'], IUPAC.protein), id=row2['Entry'], description=row2['Protein names'])
        #print mean_identity(records, record)
        records.append(record)
    SeqIO.write(records, dir_path+"/../dat/seq/"+ec+".fasta", "fasta")



if sys.argv[1].startswith("|"):
    key = sys.argv[1].replace("|","")
else:
    key = sys.argv[1]
taxonomy = sys.argv[2]

print "\t --> Trying to download sequencing data for: ", key

hit = re.match("[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+", key)
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
