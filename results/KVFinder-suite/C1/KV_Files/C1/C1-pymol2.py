import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../../.././hosts/C1.pdb", quiet=False)
cmd.load("C1.KVFinder.output.pdb", quiet=False)

cmd.hide("sticks", "C1.KVFinder.output")
cmd.show("spheres", "C1.KVFinder.output")
cmd.alter("obj C1.KVFinder.output", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj C1.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "C1.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "C1.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("C1.pse")
