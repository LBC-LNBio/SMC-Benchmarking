import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/hosts/W1.pdb", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/fpocket/W1/pocket7_vert.pqr", quiet=False)
cmd.load("/home/ABTLUS/joao.guerra/CNPEM/gyorgy.szaloki/benchmarking/results/fpocket/W1/pocket7_atm.pdb", quiet=False)

cmd.hide("sticks", "pocket7_vert")
cmd.show("surface", "pocket7_vert")
cmd.color("white", "pocket7_vert")
cmd.rebuild()

cmd.orient()

cmd.save("W1.pse")
