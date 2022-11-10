import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../.././hosts/F1.pdb", quiet=False)
cmd.load("F1.KVFinder.output.pdb", quiet=False)
cmd.load("F1.surface.pdb", quiet=False)

cmd.hide("sticks", "F1.KVFinder.output")
cmd.show("spheres", "F1.KVFinder.output")
cmd.alter("obj F1.KVFinder.output", "vdw=0.125")
cmd.alter("obj F1.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj F1.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "F1.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "F1.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("F1.pse")
