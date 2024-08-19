# Biochemistry

The gapseq database for chemical compounds and reactions originated from the [SEED database](https://modelseed.org/biochem/) which was curated and extended.
At the moment, around 1500 reactions were changed in their reversibility or stoichiometry and more than 100 reactions were newly added.
The overall goal is to provide a harmonized and completely futile cycle free database that facilitates the gapfilling process.

The help package [gapseq biochemistry](https://github.com/Waschina/gapseq.biochem) can be used to explore the gapseq database!

## List of database-specific compound definitions

###### Miscellaneous

- Hydrogen sulfite (H2S, neutral charge) cpd000239 is actually HS- (-1 charge)
- Sulfite (SO3, -2 charge) cpd00081 is actually: HSO3 (-1 charge)
- Ammonia: cpd00013 (labeled NH3 in seed with formula NH4...) do not use cpd19013

###### Electron transferring compounds

- Ferredoxin transfers each 1e- (oxidized: cpd11621, reduced: cpd11620)
- Rubredoxin transfers each 1e- (oxidized: cpd11681, reduced: cpd11651)
- Electron-Transferring-Flavoprotein (ETF) transfers 2e- (oxidized: cpd27005, reduced: cpd27006)

###### Polymers

- Starch: cpd90003 (27 x D-glucose)
- Pectin: cpd27519 (3 x Methylesterified-homogalacturonate)
- Cellulose: cpd90020 (500 x glucose)
- D-Xylan (beta-1,4): cpd90021 (500 x D-Xylose)
- 2'-Fucosyllactose: cpd90001
- Inulin: cpd28763 (11 x D-fructose + 1 x D-glucose)
