# food-web-and-ecosystem-service-robustness
code for manuscript resubmitted to Nature Communications
NCOMMS-20-38839

These R code files were used to generate the main results for our manuscript and can be edited to reproduce the results in our supplementary analyses. 
Co-author, Allison Barner, intends to publish her robustness functions as an R package in the future. Please contact the authors for access to source function files and with any questions.

All of this code was written (and tested on) for R Studio version 3.6.2. Beyond the packages included in this code, there are no external non-standard hardwares.


We have included 5 R files with code, see NCOMMS-20-38839 for data.
R files:
1. Robustness-fw-es.R includes the code for all sequences of species losses and robustness calculations for both the food webs and ecosystem services (aggregate, binary links). The code starts by reading in the data (including food webs and sequences), and creating network objects. We simulated 12 extinction sequences, some of which are specified in other data files (not included here). You can still test this code on the Most-to-least connected and random sequences, as these do not rely on exterior files. For each extinction sequence, we simulate species losses first for the food web, and second for the ecosystem services. The output of this file is 12 data files (.csv) with tracking information for when species are removed and lost. You also can get the robustness values for all food webs and services across all the sequences. This code runs once for each salt marsh - and as such, the file names should be updated to match the system of interest.

2. INDIVIDUAL_ES_ROBUSTNESS.R includes the code for all sequences of species losses and robustness calculations for the individual ecosystem service robustness. This code uses the output files that are saved from the Robustness-fw-es.R file. This code provides a robustness value for each of the ecosystem servies in each of the systems across all 12 sequences of extinctions.

3. WEIGHTED_AGGREGATE.R and WEIGHTED_INDIVIDUAL.R include the code fo all sequences of species losses and robustness calculations for the aggregate and individual weighted ecosystem service robustness, respectively. These calculations factor in species' biomasses for their interaction weights. This code provides a robustness value for all ecosystem services in aggregate and each of the ecosystem services individually (respectively) in each of the systems across all 12 sequences of extinctions.

4. PARASITE-ROBUSTNESS.R includes the code for our supplementary analysis on whether or not the inclusion of parasites and their interactions changes our results. This file follows the same format as the previous robustness code, but uses food webs with parasites (not included here).This code provides a robustness value for the food webs and ecosystem services in each of the systems across all 12 sequences of extinctions.

5. PAGERANK.R includes the code used to identify supporting species. This code runs the personalized PageRank for each ecosystem service (where the food web directed edges are flipped). The results of this code inform the indirect, supporting species sequence in the robustness analyses. Therefore, if you run this code first, you can also simulate the supporting species sequences in our robustness files. This code provides you with a list of species and their personalized PageRank score.
