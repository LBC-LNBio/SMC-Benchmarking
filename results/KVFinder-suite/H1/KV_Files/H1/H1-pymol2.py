import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/H1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/KVFinder-suite/H1/KV_Files/H1/H1.KVFinder.output.pdb", quiet=False)

cmd.hide("sticks", "H1.KVFinder.output")
cmd.show("spheres", "H1.KVFinder.output")
cmd.alter("obj H1.KVFinder.output", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj H1.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "H1.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "H1.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("H1.pse")
