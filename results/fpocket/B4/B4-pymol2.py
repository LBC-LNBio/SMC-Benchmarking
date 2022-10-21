import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/B4.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/fpocket/B4/pocket1_vert.pqr", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/fpocket/B4/pocket1_atm.pdb", quiet=False)

cmd.hide("sticks", "pocket1_vert")
cmd.show("surface", "pocket1_vert")
cmd.color("white", "pocket1_vert")
cmd.rebuild()

cmd.orient()
