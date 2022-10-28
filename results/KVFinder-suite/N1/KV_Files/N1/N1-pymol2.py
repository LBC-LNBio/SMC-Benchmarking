import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/N1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/KVFinder-suite/N1/KV_Files/N1/N1.KVFinder.output.pdb", quiet=False)

cmd.hide("sticks", "N1.KVFinder.output")
cmd.show("spheres", "N1.KVFinder.output")
cmd.alter("obj N1.KVFinder.output", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj N1.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "N1.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "N1.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("N1.pse")
