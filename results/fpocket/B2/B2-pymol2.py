import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/hosts/B2.pdb", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/fpocket/B2/pocket1_vert.pqr", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/fpocket/B2/pocket1_atm.pdb", quiet=False)

cmd.hide("sticks", "pocket1_vert")
cmd.show("surface", "pocket1_vert")
cmd.color("white", "pocket1_vert")
cmd.rebuild()

cmd.orient()

cmd.save("B2.pse")
