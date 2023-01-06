import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../.././hosts/B10.pdb", quiet=False)
cmd.load("B10.KVFinder.output.pdb", quiet=False)
cmd.load("B10.surface.pdb", quiet=False)

cmd.hide("sticks", "B10.KVFinder.output")
cmd.show("spheres", "B10.KVFinder.output")
cmd.alter("obj B10.KVFinder.output", "vdw=0.125")
cmd.alter("obj B10.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B10.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B10.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "B10.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B10.pse")
