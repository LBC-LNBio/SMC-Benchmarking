import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("C1.pocketness.pdb", quiet=False)

cmd.load("C1.pocket.pdb", quiet=False)
cmd.hide("sticks", "C1.pocket")
cmd.show("spheres", "C1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj C1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "C1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "C1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj C1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "C1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "C1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("C1.pse")
