#TITLE   Pocket invRinacc value for residues
#OUTPUTFILENAME    ./results/GHECOM/B8/B8.res
#DATE    Jan 28,2023 11:51:46
#COMMAND ghecom -M M -ipdb ./hosts/B8.pdb -opocpdb ./results/GHECOM/B8/B8.pocket.pdb -opdb ./results/GHECOM/B8/B8.pocketness.pdb -ores ./results/GHECOM/B8/B8.res -atmhet B -gw 0.8 -rlx 10
#COMMENT INPUT_RECEPTOR_PDB_FILE:./hosts/B8.pdb
#COMMENT NATOM_OF_RECEPTOR:240
#COMMENT grid_width  0.800 A VdW volume: 5274 grids, 2700.29 AAA
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
   1  - UNK  88.14  7.492 96  0   3.74   3.57   0.17   0.00   0.00   0.00
   2  - UNK  68.47  4.864  2  0   6.41   6.41   0.00   0.00   0.00   0.00
   1  - UNK  73.43  5.556  6  0   4.81   4.81   0.00   0.00   0.00   0.00
   2  - UNK  67.84  4.539  2  0   6.39   6.39   0.00   0.00   0.00   0.00
   1  - UNK  77.48  6.366  7  0   3.11   3.11   0.00   0.00   0.00   0.00
   2  - UNK  70.90  5.100  3  0   5.58   5.58   0.00   0.00   0.00   0.00
   1  - UNK  80.98  6.824  4  0   2.44   2.44   0.00   0.00   0.00   0.00
   2  - UNK  71.72  5.292  3  0   5.19   5.19   0.00   0.00   0.00   0.00
   1  - UNK  75.62  5.702  2  0   5.07   5.07   0.00   0.00   0.00   0.00
   2  - UNK  74.79  5.754  2  0   5.05   5.05   0.00   0.00   0.00   0.00
   1  - UNK  78.14  6.318  3  0   4.16   4.16   0.00   0.00   0.00   0.00
   2  - UNK  76.21  6.061  3  0   5.08   5.08   0.00   0.00   0.00   0.00
   1  - UNK  90.79  8.672  2  0   0.00   0.00   0.00   0.00   0.00   0.00
   2  - UNK  89.73  7.753 105  0   3.38   3.38   0.00   0.00   0.00   0.00
