import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("N1.pocketness.pdb", quiet=False)

cmd.load("N1.pocket.pdb", quiet=False)
cmd.hide("sticks", "N1.pocket")
cmd.show("spheres", "N1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj N1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "N1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "N1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj N1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "N1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "N1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()
