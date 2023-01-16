import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../.././hosts/B11.pdb", quiet=False)
cmd.load("B11.KVFinder.output.pdb", quiet=False)
cmd.load("B11.surface.pdb", quiet=False)

cmd.hide("sticks", "B11.KVFinder.output")
cmd.show("spheres", "B11.KVFinder.output")
cmd.alter("obj B11.KVFinder.output", "vdw=0.125")
cmd.alter("obj B11.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B11.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B11.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "B11.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B11.pse")