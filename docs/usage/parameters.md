# Parameters for network reconstruction



**gapseq** provides default values for parameters relevant for the model reconstruction process. Those defaults usually produce good reconstruction results as shown in the publication of the gapseq software. The user has the full opportunity to specify those parameters if wanted. The following notes serve as an illustration and detailed explanation of the effect of the different parameters on gapseq reconstructions. 



##### Absence/presence of reactions

If a reaction is added to a metabolic network model depends mainly on sequence homology criteria (i.e. bitscore) on the basis of sequence comparisons to reference databases. To avoid a single arbitrary cutoff  value, on which the decision of reaction absence or presence is made, gapseq uses a different approach that aims to add reactions to the final gapfilled model on the basis of the combination of sequence homology evidence and network topology. This is achieved and user-adjustable by supplying four parameters, that are used in the gapseq work-flow:



![](https://raw.githubusercontent.com/jotech/gapseq/master/docs/gfx/parameters.svg)



#### Parameter `-b` in `gapseq find` and `gapseq find-transporter`  

If a reference sequence (gene/proteine) results in a local alignment with the specified genome sequence and a bitscore higher or equal than **`b `**(and a query coverage >= 75%), associated reactions or transporters are considered to be present as likely metabolic capabilities of the organism. The default of **`b`** is 200. All bitscore statistics of this step, including *hits* with a bitscore lower than **`b`** are further user for the draft network model reconstruction and gapfilling (see below).

#### Parameters `-l` and `-u` in `gapseq draft`

`gapseq draft` constructs the first draft network model. Therefore, it adds all reactions that have been predicted as present in the previous step, either through direct sequence homology or indirectly by pathway prediction (see manuscript for details), to the network. Additionally, the command `gapseq draft`  constructs an ordered '*waiting list*' for all remaining reactions by prioritizing reactions based on their sequence evidence score. This list is the basis for the subsequent gapfilling step of gapseq.

Therefore, the bitscore is transformed to a *reaction weight* as outlined in the graph above. Reactions with an associated bitscore equal or lower than **`l`** are considered as reactions without any sequence evidence and are assigned the maximum reaction weight `w_max` . Reactions with an associated bitscore equal or above **`u`** are considered reactions with the maximum sequence evidence and, thus, are assigned a minimum weight `w_min`. In case the associated bitscore is between **`l`** and **`u`** the weights are calculated using a linear function with weight(`l`) = `w_max` and weight(`u`) = `w_min`.

We recommend using a value between 50 (default) and 100 for **`l`** and 200 (default) for **`u`**. If you want to use the default, no indication of **`l`** and **`u`** are required in the command line. If you wish to adjust the parameters, we further recommend using for **`u`** in `gapseq draft` the same value as for the parameter **`b`** in `gapseq find` (defaults: 200).

#### Parameter `-b` in `gapseq fill`

The parameter **`b`** in `gapseq fill` denotes the bitscore threshold above which associated reactions are considered as core candidate reactions for gapfilling. Those reactions are considered for the gapfilling steps 1 (enabling growth), 2/2b (Biomass component biosynthesis gapfilling), 3 (alternative energy sources gapfilling), and 4 (potential metabolic (by)-products). Reactions with bitscores below **`b`** are only considered for the initial gapfilling step 1 (enable growth).

For details concerning the gapfilling algorithm, please refer to the gapseq manuscript and/or source code. 