import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/W1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/fpocket/W1/pocket7_vert.pqr", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/fpocket/W1/pocket7_atm.pdb", quiet=False)

cmd.hide("sticks", "pocket7_vert")
cmd.show("surface", "pocket7_vert")
cmd.color("white", "pocket7_vert")
cmd.rebuild()

cmd.orient()

cmd.save("W1.pse")
