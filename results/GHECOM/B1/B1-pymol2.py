import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("B1.pocketness.pdb", quiet=False)

cmd.load("B1.pocket.pdb", quiet=False)
cmd.hide("sticks", "B1.pocket")
cmd.show("spheres", "B1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "B1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj B1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "B1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()
