import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B14.pdb", quiet=False)
cmd.load("pocket3_vert.pqr", quiet=False)
cmd.load("pocket3_atm.pdb", quiet=False)

cmd.hide("sticks", "pocket3_vert")
cmd.show("surface", "pocket3_vert")
cmd.color("white", "pocket3_vert")
cmd.rebuild()

cmd.orient()

cmd.save("B14.pse")
