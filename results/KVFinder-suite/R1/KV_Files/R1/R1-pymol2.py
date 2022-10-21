import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/R1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/KVFinder-suite/R1/KV_Files/R1.KVFinder.output.pdb", quiet=False)

cmd.hide("sticks", "R1.KVFinder.output")
cmd.show("spheres", "R1.KVFinder.output")
cmd.alter("obj R1.KVFinder.output", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj R1.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "R1.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "R1.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()