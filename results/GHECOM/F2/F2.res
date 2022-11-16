#TITLE   Pocket invRinacc value for residues
#OUTPUTFILENAME    ./results/GHECOM/F2/F2.res
#DATE    Nov 16,2022 11:3:45
#COMMAND ghecom -M M -ipdb ./hosts/F2.pdb -opocpdb ./results/GHECOM/F2/F2.pocket.pdb -opdb ./results/GHECOM/F2/F2.pocketness.pdb -ores ./results/GHECOM/F2/F2.res -atmhet B -gw 0.8 -rlx 20.0
#COMMENT INPUT_RECEPTOR_PDB_FILE:./hosts/F2.pdb
#COMMENT NATOM_OF_RECEPTOR:1290
#COMMENT grid_width  0.800 A VdW volume: 31703 grids, 16231.94 AAA
#COMMENT MIN_RLARGE:2.000000
#COMMENT MAX_RLARGE:20.000000
#COMMENT BIN_RLARGE:0.500000
#COMMENT RSMALL: 1.870
#MULSC_PROBE: NRADIUS 37
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
#MULSC_PROBE: RADIUS  18 th 10.500 Ang  0.095 1/Ang 
#MULSC_PROBE: RADIUS  19 th 11.000 Ang  0.091 1/Ang 
#MULSC_PROBE: RADIUS  20 th 11.500 Ang  0.087 1/Ang 
#MULSC_PROBE: RADIUS  21 th 12.000 Ang  0.083 1/Ang 
#MULSC_PROBE: RADIUS  22 th 12.500 Ang  0.080 1/Ang 
#MULSC_PROBE: RADIUS  23 th 13.000 Ang  0.077 1/Ang 
#MULSC_PROBE: RADIUS  24 th 13.500 Ang  0.074 1/Ang 
#MULSC_PROBE: RADIUS  25 th 14.000 Ang  0.071 1/Ang 
#MULSC_PROBE: RADIUS  26 th 14.500 Ang  0.069 1/Ang 
#MULSC_PROBE: RADIUS  27 th 15.000 Ang  0.067 1/Ang 
#MULSC_PROBE: RADIUS  28 th 15.500 Ang  0.065 1/Ang 
#MULSC_PROBE: RADIUS  29 th 16.000 Ang  0.062 1/Ang 
#MULSC_PROBE: RADIUS  30 th 16.500 Ang  0.061 1/Ang 
#MULSC_PROBE: RADIUS  31 th 17.000 Ang  0.059 1/Ang 
#MULSC_PROBE: RADIUS  32 th 17.500 Ang  0.057 1/Ang 
#MULSC_PROBE: RADIUS  33 th 18.000 Ang  0.056 1/Ang 
#MULSC_PROBE: RADIUS  34 th 18.500 Ang  0.054 1/Ang 
#MULSC_PROBE: RADIUS  35 th 19.000 Ang  0.053 1/Ang 
#MULSC_PROBE: RADIUS  36 th 19.500 Ang  0.051 1/Ang 
#MULSC_PROBE: RADIUS  37 th 20.000 Ang  0.050 1/Ang 
#OUTSIDE_OF_MAX_RADIUS      20.500 Ang  0.049 1/Ang 
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
   1  - UNK 100.00 14.847 1290  0   8.72   8.72   0.00   0.00   0.00   0.00
