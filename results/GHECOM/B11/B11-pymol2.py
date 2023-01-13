import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("B11.pocketness.pdb", quiet=False)

cmd.load("B11.pocket.pdb", quiet=False)
cmd.hide("sticks", "B11.pocket")
cmd.show("spheres", "B11.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B11.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B11.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "B11.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj B11.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B11.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "B11.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B11.pse")
