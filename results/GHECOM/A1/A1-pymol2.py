import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("A1.pocketness.pdb", quiet=False)

cmd.load("A1.pocket.pdb", quiet=False)
cmd.hide("sticks", "A1.pocket")
cmd.show("spheres", "A1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj A1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "A1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "A1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj A1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "A1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "A1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("A1.pse")
