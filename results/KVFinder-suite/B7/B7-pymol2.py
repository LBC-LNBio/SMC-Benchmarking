import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../.././hosts/B7.pdb", quiet=False)
cmd.load("B7.KVFinder.output.pdb", quiet=False)
cmd.load("B7.surface.pdb", quiet=False)

cmd.hide("sticks", "B7.KVFinder.output")
cmd.show("spheres", "B7.KVFinder.output")
cmd.alter("obj B7.KVFinder.output", "vdw=0.125")
cmd.alter("obj B7.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B7.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B7.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "B7.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B7.pse")
