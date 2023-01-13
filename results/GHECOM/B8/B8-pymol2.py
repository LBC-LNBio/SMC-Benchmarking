import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("B8.pocketness.pdb", quiet=False)

cmd.load("B8.pocket.pdb", quiet=False)
cmd.hide("sticks", "B8.pocket")
cmd.show("spheres", "B8.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B8.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B8.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "B8.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj B8.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B8.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "B8.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B8.pse")
