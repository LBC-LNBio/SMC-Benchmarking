import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("F1.pocketness.pdb", quiet=False)

cmd.load("F1.pocket.pdb", quiet=False)
cmd.hide("sticks", "F1.pocket")
cmd.show("spheres", "F1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj F1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "F1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "F1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj F1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "F1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "F1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("F1.pocketness.pse")
