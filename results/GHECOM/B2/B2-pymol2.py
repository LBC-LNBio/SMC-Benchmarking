import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("B2.pocketness.pdb", quiet=False)

cmd.load("B2.pocket.pdb", quiet=False)
cmd.hide("sticks", "B2.pocket")
cmd.show("spheres", "B2.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B2.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B2.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "B2.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj B2.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B2.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "B2.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()
