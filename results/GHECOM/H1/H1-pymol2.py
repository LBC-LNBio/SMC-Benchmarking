import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("H1.pocketness.pdb", quiet=False)

cmd.load("H1.pocket.pdb", quiet=False)
cmd.hide("sticks", "H1.pocket")
cmd.show("spheres", "H1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj H1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "H1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "H1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj H1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "H1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "H1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("H1.pocketness.pse")
