import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/A1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/KVFinder-suite/A1/KV_Files/A1/A1.KVFinder.output.pdb", quiet=False)

cmd.hide("sticks", "A1.KVFinder.output")
cmd.show("spheres", "A1.KVFinder.output")
cmd.alter("obj A1.KVFinder.output", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj A1.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "A1.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "A1.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("A1.pse")
