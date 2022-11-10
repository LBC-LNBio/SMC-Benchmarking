import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/hosts/O1.pdb", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/KVFinder-suite/O1/O1.KVFinder.output.pdb", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/KVFinder-suite/O1/O1.surface.pdb", quiet=False)

cmd.hide("sticks", "O1.KVFinder.output")
cmd.show("spheres", "O1.KVFinder.output")
cmd.alter("obj O1.KVFinder.output", "vdw=0.125")
cmd.alter("obj O1.surface", "vdw=0.125")
cmd.rebuild()

stored.b = []
cmd.iterate("(obj O1.KVFinder.output)", "stored.b.append(b)")
cmd.spectrum("b", "blue_white_red", "O1.KVFinder.output", [min(stored.b), max(stored.b)])
cmd.ramp_new("Depth", "O1.KVFinder.output", [min(stored.b), max(stored.b)], ["blue", "white", "red"])

cmd.orient()

cmd.save("O1.pse")
