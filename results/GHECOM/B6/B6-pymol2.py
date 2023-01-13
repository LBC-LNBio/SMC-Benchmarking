import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("B6.pocketness.pdb", quiet=False)

cmd.load("B6.pocket.pdb", quiet=False)
cmd.hide("sticks", "B6.pocket")
cmd.show("spheres", "B6.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B6.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B6.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "B6.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj B6.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B6.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "B6.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B6.pse")
