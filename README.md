# social-infection

This repository supports [Facultative symbiont virulence determines horizontal transmission rate without host strain specificity in _Dictyostelium discoideum_ social amoebas](https://doi.org/10.1093/evlett/qrae001)

It contains two directories:

**experiment**, which contains
*   **bonniea_transmission_2023.R** - R code used to run statistical analyses using the data files and generate figures
*   **host_fitness.clean.20220711.txt** - Data from host fitness experiment with columns in the following order
    *   host identity
    *   symbiont identity
    *   MOI (multiplicity of infection) of host-symbiont pairing
    *   date of experiment
    *   estimated total number of host spores produced by sample
    *   percent of host spored infected by symbiont in sample
    *   percent of spores produced by sample relative to mean of uninfected controls
    *   type of host-symbiont pairing
*   **symbiont_transmission.clean.20220711.txt** - Data from symbiont transmission experiment with columns in the following order
    *   host identity
    *   symbiont identity
    *   MOI of host-symbiont pairing
    *   date of experiment
    *   percent of host spores infected by symbiont in sample
    *   percent of infected spores that were previously uninfected
    *   type of host-symbiont pairing

**genomics**, which contains
*   **GFF** files for the three *P. bonniea* genomes - Predicted intact gene for each symbiont genome
*   **assembly_to_set_operations.txt** - Description and command-line code for how the GFF files were generated including long read assembly, short read polishing, gene annotation, pseudogene prediction, and set operations
*   **virulence_candidates.faa** - Protein sequences for the candidate virulence factors identified in the study
