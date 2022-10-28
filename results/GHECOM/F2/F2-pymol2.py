import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("F2.pocketness.pdb", quiet=False)

cmd.load("F2.pocket.pdb", quiet=False)
cmd.hide("sticks", "F2.pocket")
cmd.show("spheres", "F2.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj F2.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "F2.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "F2.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj F2.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "F2.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "F2.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("F2.pocketness.pse")
