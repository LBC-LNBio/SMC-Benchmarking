import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../.././hosts/B12.pdb", quiet=False)
cmd.load("B12.KVFinder.output.pdb", quiet=False)
cmd.load("B12.surface.pdb", quiet=False)

cmd.hide("sticks", "B12.KVFinder.output")
cmd.show("spheres", "B12.KVFinder.output")
cmd.alter("obj B12.KVFinder.output", "vdw=0.125")
cmd.alter("obj B12.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B12.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B12.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "B12.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B12.pse")
