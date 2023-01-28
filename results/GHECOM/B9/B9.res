#TITLE   Pocket invRinacc value for residues
#OUTPUTFILENAME    ./results/GHECOM/B9/B9.res
#DATE    Jan 28,2023 11:51:46
#COMMAND ghecom -M M -ipdb ./hosts/B9.pdb -opocpdb ./results/GHECOM/B9/B9.pocket.pdb -opdb ./results/GHECOM/B9/B9.pocketness.pdb -ores ./results/GHECOM/B9/B9.res -atmhet B -gw 0.8 -rlx 10
#COMMENT INPUT_RECEPTOR_PDB_FILE:./hosts/B9.pdb
#COMMENT NATOM_OF_RECEPTOR:156
#COMMENT grid_width  0.800 A VdW volume: 3428 grids, 1755.14 AAA
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
   1  - UNK  70.33  5.566  5  0   4.94   4.94   0.00   0.00   0.00   0.00
   2  - UNK  71.04  5.663  5  0   4.99   4.99   0.00   0.00   0.00   0.00
   3  - UNK  69.60  5.534  5  0   5.02   5.02   0.00   0.00   0.00   0.00
   4  - UNK  70.07  5.616  5  0   5.24   5.24   0.00   0.00   0.00   0.00
   1  - UNK  68.68  5.449  8  0   4.91   4.91   0.00   0.00   0.00   0.00
   2  - UNK  68.05  5.412  8  0   5.03   5.03   0.00   0.00   0.00   0.00
   3  - UNK  68.74  5.514  8  0   4.85   4.85   0.00   0.00   0.00   0.00
   4  - UNK  68.60  5.451  8  0   5.01   5.01   0.00   0.00   0.00   0.00
   1  - UNK  69.67  5.351 10  0   4.86   4.86   0.00   0.00   0.00   0.00
   2  - UNK  69.52  5.367 10  0   4.96   4.96   0.00   0.00   0.00   0.00
   3  - UNK  70.05  5.488 10  0   4.70   4.70   0.00   0.00   0.00   0.00
   4  - UNK  67.67  5.134  9  0   5.01   5.01   0.00   0.00   0.00   0.00
   3  - UNK  60.61  4.522  3  0   2.55   2.55   0.00   0.00   0.00   0.00
   4  - UNK  66.88  5.416  3  0   0.83   0.83   0.00   0.00   0.00   0.00
   2  - UNK  66.58  5.515  3  0   0.90   0.90   0.00   0.00   0.00   0.00
   3  - UNK  83.99  8.368  2  0   0.00   0.00   0.00   0.00   0.00   0.00
   2  - UNK  73.60  6.447  2  0   0.15   0.15   0.00   0.00   0.00   0.00
   4  - UNK  81.91  8.048  2  0   0.00   0.00   0.00   0.00   0.00   0.00
   1  - UNK  83.09  8.218  2  0   0.00   0.00   0.00   0.00   0.00   0.00
   2  - UNK  56.04  3.797  2  0   2.89   2.89   0.00   0.00   0.00   0.00
   3  - UNK  70.58  5.879  3  0   3.19   3.19   0.00   0.00   0.00   0.00
   4  - UNK  61.60  4.473  2  0   3.99   3.99   0.00   0.00   0.00   0.00
   3  - UNK  80.71  7.927  2  0   0.00   0.00   0.00   0.00   0.00   0.00
   2  - UNK  65.75  5.271  2  0   5.69   5.69   0.00   0.00   0.00   0.00
   4  - UNK  72.62  6.253  4  0   3.31   3.31   0.00   0.00   0.00   0.00
   1  - UNK  67.27  5.338  5  0   3.43   3.43   0.00   0.00   0.00   0.00
   3  - UNK  73.13  6.700  3  0   0.72   0.72   0.00   0.00   0.00   0.00
   2  - UNK  71.34  6.446  3  0   0.65   0.65   0.00   0.00   0.00   0.00
   1  - UNK  59.72  4.976  1  0   1.02   1.02   0.00   0.00   0.00   0.00
   4  - UNK  55.73  3.778  1  0   1.91   1.91   0.00   0.00   0.00   0.00
   2  - UNK  73.78  6.988  2  0   0.07   0.07   0.00   0.00   0.00   0.00
   4  - UNK  78.26  7.244  3  0   0.10   0.10   0.00   0.00   0.00   0.00
   2  - UNK  72.78  7.001  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   3  - UNK  62.55  4.822  1  0   0.96   0.96   0.00   0.00   0.00   0.00
   1  - UNK  74.00  6.482  4  0   0.26   0.26   0.00   0.00   0.00   0.00
   4  - UNK  71.84  6.858  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   3  - UNK  70.44  5.942  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   1  - UNK  71.95  6.829  2  0   0.00   0.00   0.00   0.00   0.00   0.00
   2  - UNK  84.81  8.477  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   1  - UNK  69.93  6.406  2  0   0.21   0.21   0.00   0.00   0.00   0.00
   3  - UNK  79.23  7.276  1  0   0.00   0.00   0.00   0.00   0.00   0.00
   4  - UNK  62.08  5.366  1  0   3.15   3.15   0.00   0.00   0.00   0.00
