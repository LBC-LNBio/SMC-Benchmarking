import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/hosts/F2.pdb", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/fpocket/F2/pocket37_vert.pqr", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/fpocket/F2/pocket37_atm.pdb", quiet=False)

cmd.hide("sticks", "pocket37_vert")
cmd.show("surface", "pocket37_vert")
cmd.color("white", "pocket37_vert")
cmd.rebuild()

cmd.orient()

cmd.save("F2.pse")
