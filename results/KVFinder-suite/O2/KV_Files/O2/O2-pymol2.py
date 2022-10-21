import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/O2.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/KVFinder-suite/O2/KV_Files/O2/O2.KVFinder.output.pdb", quiet=False)

cmd.hide("sticks", "O2.KVFinder.output")
cmd.show("spheres", "O2.KVFinder.output")
cmd.alter("obj O2.KVFinder.output", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj O2.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "O2.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "O2.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()
