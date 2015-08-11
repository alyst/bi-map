For protein modules inference from AP-MS data BI-MAP requires the following information:
**description of experiments** measurements data
**properties of proteins.**

## Plain text tables data format ##

One could provide the necessary information via plain-text table files (tab-separated). The corresponding `BIMAP-sampler` options are given below

**--exp\_design\_file <experimental design filename>**
Experimental design table should consist of the following columns:
**Bait Accession Code** Biological Sample Identifier
**MS run Identifier** MS run normalization constant

_Example_
```
Bait_AC  Sample_ID  MSrun_ID     mult
Prot1    Prot1_WT   Prot1_WT-1   1.0
Prot1    Prot1_WT   Prot1_WT-2   0.9
Prot1    Prot1_Mut  Prot1_Mut-1  0.5
Prot2    Prot2_WT   Prot2_WT-1   1.0
```

**--measurements\_file <AP-MS measurements file>**
**MS run identifier** Prey protein Accession Code
**spectral counts**

_Example_
```
MSrun_ID    Prey_AC  SC
Prot1_WT-1  Prot1    500
Prot1_WT-1  Prot3    10
Prot2_WT-1  Prot1    5
Prot2_WT-1  Prot2    600
```
/Note/: for correct mapping between the bait and prey proteins the same accession codes should be used in measurements and experimental design tables

**--proteins\_file <proteins file>**
**Accession Code (same as for `Prey_AC`)** AA Sequence Length

_Example_
```
AC     SeqLength
Prot1  1000
Prot2  1200
```

## OPA XML file format ##

**--input\_file <OPA file>**

The data could also be provided in OPA XML file format. This file could be easily prepared in R with the help of `OPASaveData()` function from `RBIMAP` library. ...

See `samples/tip49/analysis.R` script for further details of preparing AP-MS data and running BI-MAP within R session.

## Mapping between baits and preys ##
Depending on AP-MS experiments baits could be protein as well as the other molecules (DNA). In the latter case the mapping between proteins and baits is not possible. To disable the mapping an option
```
--map_baits_to_preys=0
```
has to be specified. This also instructs BI-MAP not to use the experimental design component of the likelihood function.