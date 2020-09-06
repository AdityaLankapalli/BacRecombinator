# Ttrip-BAC : Tree Topology based Recombination Identification Program for BACteria
Ttrip-BAC tries to detect recombination patterns in bacterial enomes. It is an automated pipeline that estimates incongruous phylogenetic patterns across a bacterial genome. Here the pipeline runs on a specific bacterial genome (target genome) and tries to quantify and visualise the mosaic genomic regions.
It utilises RAxML to build a maximum likelihood tree (backbone tree) for a set of representative genomes. Tree topology models are generated based on the backbone tree for each target genome. TREE-PUZZLE generates maximum likelihood trees for each alignment given a set of tree topology models. These likelihood values are compared and evaluated. 

Ttrip-BAC is a simple, fast and first hand pipeline to quanitfy any variation in phylogenetic patterns. It has been applied to multiple datasets from diverse contexts. The pipeline is specifically applicable to ancient pathogen genomes, where recombination hinders a clonal inference. 

## unpublished

Usage of the pipeline
```
./main.py -i <fasta file> -g <genomes list> -pre <prefix> -sw <num> -ovl <num> -nt <num> -s <seed num>

Example:
./main.py -i fullAlignment.fasta -g genomes.txt -pre PeCan4RUN -sw 1000 -ovl 1000 -nt 4 -s 18494

```
## help and options 

```
./main.py -h

-h, --help            show this help message and exit
-i INPUT_FASTA, --input INPUT_FASTA
                      Input multifasta
-t BACKBONE_TREE, --trees BACKBONE_TREE
                      Backbone topology/topologies file (Not applicable if raxml is active)
-g REP_GENOMES, --rep_genomes REP_GENOMES
                      List of represntative genomes for backbone topology and target genome
-b BOOTSTRAP, --bootstrap BOOTSTRAP
                      Bootstrap replicates for Backbone Topology
-p PARAMS_FILE, --parameterfile PARAMS_FILE
                      Parameters file (Not applicable if raxml active)
-pre PREFIX, --prefix PREFIX
                      prefix for output
-sw SLIDING_WINDOW, --windowlength SLIDING_WINDOW
                      Size of sliding window for each alginment block
-ovl SLIDING_OVERLAP, --overlap SLIDING_OVERLAP
                      Overlap size of sliding windows
-s RANDOMSEED, --seed RANDOMSEED
                      Seed value

-op OPTIONS, --options OPTIONS
                      Additional options from the user
-pro OPTIONS, --program OPTIONS
                      Program for Backbone Topology and topology testing
-nt THREADS, --threads THREADS
                      Number of threads
-r, --redo            Re run the analaysis
-incr, --rootinclude  include root tree
```


## Work under progress

