# Growth medium definition

#### Growth medium file format

A user gapfill medium can be defined in a 3-column comma-(csv) or tab-(tsv)-separated table.

- **1st column** "compounds": Should list the compounds' identifier
- **2nd column** "name": List of compound names
- **3rd column** "maxFlux": Defines the maximum uptake rate of the compound in `mmol/gDW/hr`  

Example: M9-like medium (aerobic) with limited Glucose+Oxygen uptake and virtually unlimited salts/ions.

```
compounds,name,maxFlux
cpd00001,H2O,100
cpd00007,O2,10
cpd00009,Phosphate,100
cpd00027,D-Glucose,5
cpd00030,Mn2+,100
cpd00034,Zn2+,100
cpd00048,Sulfate,100
cpd00058,Cu2+,100
cpd00063,Ca2+,100
cpd00067,H+,100
cpd00099,Cl-,100
cpd00149,Co2+,100
cpd00205,K+,100
cpd00254,Mg,100
cpd00531,Hg2+,100
cpd00971,Na+,100
cpd01012,Cd2+,100
cpd01048,Arsenate,100
cpd10515,Fe2+,100
cpd10516,fe3,100
cpd11595,chromate,100
cpd00013,NH3,100
```

#### List of potential nutrients

A table with all nutrients (IDs and names), that can be part of a growth medium definition is stored in your gapseq installation (`dat/nutrients.tsv`) or accessed also [directly on github](https://github.com/jotech/gapseq/blob/master/dat/nutrients.tsv).

A list of database-specific compounds, which might diverge from metabolite entries in other databases, is given in the [database section](../database/biochemistry).
