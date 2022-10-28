#TITLE   Pocket invRinacc value for residues
#OUTPUTFILENAME    ./results/GHECOM/B2/B2.res
#DATE    Oct 28,2022 14:53:34
#COMMAND ghecom -M M -ipdb ./data/B2.pdb -opocpdb ./results/GHECOM/B2/B2.pocket.pdb -opdb ./results/GHECOM/B2/B2.pocketness.pdb -ores ./results/GHECOM/B2/B2.res -atmhet B -gw 0.8 -rlx 10.0
#COMMENT INPUT_RECEPTOR_PDB_FILE:./data/B2.pdb
#COMMENT NATOM_OF_RECEPTOR:226
#COMMENT grid_width  0.800 A VdW volume: 4580 grids, 2344.96 AAA
#COMMENT MIN_RLARGE:2.000000
#COMMENT MAX_RLARGE:10.000000
#COMMENT BIN_RLARGE:0.500000
#COMMENT RSMALL: 1.870
#MULSC_PROBE: NRADIUS 17
#MULSC_PROBE: RADIUS   1 th  2.000 Ang  0.500 1/Ang 
#MULSC_PROBE: RADIUS   2 th  2.500 Ang  0.400 1/Ang 
#MULSC_PROBE: RADIUS   3 th  3.000 Ang  0.333 1/Ang 
#MULSC_PROBE: RADIUS   4 th  3.500 Ang  0.286 1/Ang 
#MULSC_PROBE: RADIUS   5 th  4.000 Ang  0.250 1/Ang 
#MULSC_PROBE: RADIUS   6 th  4.500 Ang  0.222 1/Ang 
#MULSC_PROBE: RADIUS   7 th  5.000 Ang  0.200 1/Ang 
#MULSC_PROBE: RADIUS   8 th  5.500 Ang  0.182 1/Ang 
#MULSC_PROBE: RADIUS   9 th  6.000 Ang  0.167 1/Ang 
#MULSC_PROBE: RADIUS  10 th  6.500 Ang  0.154 1/Ang 
#MULSC_PROBE: RADIUS  11 th  7.000 Ang  0.143 1/Ang 
#MULSC_PROBE: RADIUS  12 th  7.500 Ang  0.133 1/Ang 
#MULSC_PROBE: RADIUS  13 th  8.000 Ang  0.125 1/Ang 
#MULSC_PROBE: RADIUS  14 th  8.500 Ang  0.118 1/Ang 
#MULSC_PROBE: RADIUS  15 th  9.000 Ang  0.111 1/Ang 
#MULSC_PROBE: RADIUS  16 th  9.500 Ang  0.105 1/Ang 
#MULSC_PROBE: RADIUS  17 th 10.000 Ang  0.100 1/Ang 
#OUTSIDE_OF_MAX_RADIUS      10.500 Ang  0.095 1/Ang 
#COLUMN  1|RNUM              |Residue Number
#COLUMN  2|CHAIN             |Chain Identifier
#COLUMN  3|RES               |Three-letter residue name
#COLUMN  4|shellAcc          |shell accessibility (%)
#COLUMN  5|Rinacc            |averaged Rinaccess (A)
#COLUMN  6|Natom             |Number of atoms
#COLUMN  7|Natom_contact     |Number of contacting atoms with ligand
#COLUMN  8|pocketness        |sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%)
#COLUMN  9|pocketness_clus 1 |sum of 1/[Rpocket_for_pocketcluster 1] /(1/[Rmin]*[vol of shell]) (%)
#COLUMN 10|pocketness_clus 2 |sum of 1/[Rpocket_for_pocketcluster 2] /(1/[Rmin]*[vol of shell]) (%)
#COLUMN 11|pocketness_clus 3 |sum of 1/[Rpocket_for_pocketcluster 3] /(1/[Rmin]*[vol of shell]) (%)
#COLUMN 12|pocketness_clus 4 |sum of 1/[Rpocket_for_pocketcluster 4] /(1/[Rmin]*[vol of shell]) (%)
#COLUMN 13|pocketness_clus 5 |sum of 1/[Rpocket_for_pocketcluster 5] /(1/[Rmin]*[vol of shell]) (%)
   1  - UNK  60.79  4.804 14  0   1.61   1.61   0.00   0.00   0.00   0.00
   2  - UNK  42.34  2.854  1  0   0.63   0.63   0.00   0.00   0.00   0.00
   1  - UNK  56.03  4.668  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   2  - UNK  54.18  4.413  3  0   1.13   1.13   0.00   0.00   0.00   0.00
   3  - UNK  33.67  1.760  1  0   1.59   1.59   0.00   0.00   0.00   0.00
   1  - UNK  62.24  4.908 14  0   1.55   1.55   0.00   0.00   0.00   0.00
   2  - UNK  59.83  4.896  2  0   0.98   0.98   0.00   0.00   0.00   0.00
   1  - UNK  76.30  6.651 32  0   1.64   1.64   0.00   0.00   0.00   0.00
   2  - UNK  56.51  4.895  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   1  - UNK  48.24  3.726  1  0   0.58   0.58   0.00   0.00   0.00   0.00
   4  - UNK  57.47  4.753  1  0   0.19   0.19   0.00   0.00   0.00   0.00
   1  - UNK  53.45  4.074  1  0   1.77   1.77   0.00   0.00   0.00   0.00
   2  - UNK  67.61  5.455  1  0   0.24   0.24   0.00   0.00   0.00   0.00
   3  - UNK  62.96  4.728  1  0   1.72   1.72   0.00   0.00   0.00   0.00
   1  - UNK  64.41  4.931 23  0   2.52   2.52   0.00   0.00   0.00   0.00
   5  - UNK  54.31  3.160  1  0   6.23   6.23   0.00   0.00   0.00   0.00
   6  - UNK  54.31  3.309  1  0   3.99   3.99   0.00   0.00   0.00   0.00
   5  - UNK  57.09  3.398  1  0   5.68   5.68   0.00   0.00   0.00   0.00
   1  - UNK  71.35  5.716  7  0   3.42   3.42   0.00   0.00   0.00   0.00
   7  - UNK  50.04  2.477  1  0   5.49   5.49   0.00   0.00   0.00   0.00
   1  - UNK  54.05  3.216  4  0   5.51   5.51   0.00   0.00   0.00   0.00
   2  - UNK  74.68  6.591  4  0   0.72   0.72   0.00   0.00   0.00   0.00
   1  - UNK  85.35  7.826 68  0   1.52   1.52   0.00   0.00   0.00   0.00
   5  - UNK  57.52  3.785  8  0   4.39   4.39   0.00   0.00   0.00   0.00
   1  - UNK  74.90  6.868  5  0   0.20   0.20   0.00   0.00   0.00   0.00
   2  - UNK  70.57  6.192  3  0   0.53   0.53   0.00   0.00   0.00   0.00
   6  - UNK  61.58  4.837  1  0   3.88   3.88   0.00   0.00   0.00   0.00
   1  - UNK  59.49  4.970  2  0   0.48   0.48   0.00   0.00   0.00   0.00
   2  - UNK  65.18  5.060  1  0   1.25   1.25   0.00   0.00   0.00   0.00
   3  - UNK  64.58  5.170  1  0   0.62   0.62   0.00   0.00   0.00   0.00
   1  - UNK  60.66  4.237  3  0   3.69   3.69   0.00   0.00   0.00   0.00
   5  - UNK  54.51  3.294  3  0   5.86   5.86   0.00   0.00   0.00   0.00
   6  - UNK  66.09  5.677  2  0   2.84   2.84   0.00   0.00   0.00   0.00
   1  - UNK  54.12  3.255  3  0   6.00   6.00   0.00   0.00   0.00   0.00
   2  - UNK  66.19  5.533  1  0   0.57   0.57   0.00   0.00   0.00   0.00
   1  - UNK  81.92  7.546  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   8  - UNK  74.72  7.336  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   1  - UNK  66.14  6.032  2  0   0.00   0.00   0.00   0.00   0.00   0.00
   5  - UNK  55.25  3.338  4  0   5.68   5.68   0.00   0.00   0.00   0.00
   1  - UNK  43.51  2.726  1  0   2.14   2.14   0.00   0.00   0.00   0.00
