set dir "/home/jvsguerra/remote-repos/SMC-Benchmarking/results/CAVER/B3/caver_output/1/inputs"

mol load pdb ${dir}/B3.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

