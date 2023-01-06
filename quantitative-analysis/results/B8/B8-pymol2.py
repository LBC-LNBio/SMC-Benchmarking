import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../.././hosts/B8.pdb", quiet=False)
cmd.load("B8.KVFinder.output.pdb", quiet=False)
cmd.load("B8.surface.pdb", quiet=False)

cmd.hide("sticks", "B8.KVFinder.output")
cmd.show("spheres", "B8.KVFinder.output")
cmd.alter("obj B8.KVFinder.output", "vdw=0.125")
cmd.alter("obj B8.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B8.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B8.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "B8.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B8.pse")
