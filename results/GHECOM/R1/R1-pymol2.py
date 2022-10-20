import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("R1.pocketness.pdb", quiet=False)

cmd.load("R1.pocket.pdb", quiet=False)
cmd.hide("sticks", "R1.pocket")
cmd.show("spheres", "R1.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj R1.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "R1.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "R1.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj R1.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "R1.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "R1.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()
