#set dir "/home/jvsguerra/remote-repos/SMC-Benchmarking/results/CAVER/N1/caver_output/1/inputs"

mol load pdb ../data/N1.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 
