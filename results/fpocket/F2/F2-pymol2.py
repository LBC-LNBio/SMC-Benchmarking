import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/F2.pdb", quiet=False)
cmd.load("pocket37_vert.pqr", quiet=False)
cmd.load("pocket37_atm.pdb", quiet=False)

cmd.hide("sticks", "pocket37_vert")
cmd.show("surface", "pocket37_vert")
cmd.color("white", "pocket37_vert")
cmd.rebuild()

cmd.orient()

cmd.save("F2.pse")
