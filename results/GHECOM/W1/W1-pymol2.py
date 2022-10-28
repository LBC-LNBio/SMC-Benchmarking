import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("W1.pocketness.pdb", quiet=False)

cmd.load("W1.pocket.pdb", quiet=False)
cmd.hide("sticks", "W1.pocket")
cmd.show("spheres", "W1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj W1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "W1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "W1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj W1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "W1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "W1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("W1.pocketness.pse")
