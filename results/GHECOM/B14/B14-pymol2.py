import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("B14.pocketness.pdb", quiet=False)

cmd.load("B14.pocket.pdb", quiet=False)
cmd.hide("sticks", "B14.pocket")
cmd.show("spheres", "B14.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B14.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B14.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "B14.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj B14.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B14.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "B14.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B14.pse")
