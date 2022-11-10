import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/hosts/B4.pdb", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/KVFinder-suite/B4/B4.KVFinder.output.pdb", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/KVFinder-suite/B4/B4.surface.pdb", quiet=False)

cmd.hide("sticks", "B4.KVFinder.output")
cmd.show("spheres", "B4.KVFinder.output")
cmd.alter("obj B4.KVFinder.output", "vdw=0.125")
cmd.alter("obj B4.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj B4.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "B4.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "B4.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("B4.pse")
