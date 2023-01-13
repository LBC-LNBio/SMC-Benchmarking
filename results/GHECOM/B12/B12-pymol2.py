import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("B12.pocketness.pdb", quiet=False)

cmd.load("B12.pocket.pdb", quiet=False)
cmd.hide("sticks", "B12.pocket")
cmd.show("spheres", "B12.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B12.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B12.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "B12.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj B12.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B12.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "B12.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B12.pse")
