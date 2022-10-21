import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/F2.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/KVFinder-suite/F2/F2.KVFinder.output.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/KVFinder-suite/F2/F2.surface.pdb", quiet=False)

cmd.hide("sticks", "F2.KVFinder.output")
cmd.show("spheres", "F2.KVFinder.output")
cmd.alter("obj F2.KVFinder.output", "vdw=0.3")
cmd.alter("obj F2.surface", "vdw=0.3")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj F2.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "F2.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "F2.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()
