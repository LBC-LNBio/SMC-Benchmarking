#set dir "/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/CAVER/O2/caver_output/1/inputs"

mol load pdb ../data/O2.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 
