import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("O2.pocketness.pdb", quiet=False)

cmd.load("O2.pocket.pdb", quiet=False)
cmd.hide("sticks", "O2.pocket")
cmd.show("spheres", "O2.pocket")
cmd.alter("resn GRD", "vdw=0.4")
cmd.alter("resn CEN", "vdw=0.1")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj O2.pocket)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "O2.pocket", [min(stored.b), max(stored.b)])
cmd.ramp_new("Rinaccess", "O2.pocket", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

stored.b = []
cmd.iterate("(obj O2.pocketness)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "O2.pocketness", [min(stored.b), max(stored.b)])
cmd.ramp_new("Pocketness", "O2.pocketness", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("O2.pse")
