import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../../.././hosts/B2.pdb", quiet=False)
cmd.load("B2.KVFinder.output.pdb", quiet=False)

cmd.hide("sticks", "B2.KVFinder.output")
cmd.show("spheres", "B2.KVFinder.output")
cmd.alter("obj B2.KVFinder.output", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B2.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B2.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "B2.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B2.pse")
