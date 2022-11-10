import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/O2.pdb", quiet=False)
cmd.load("O2.pywindow.pdb", quiet=False)

stored.b = []
cmd.iterate("(resn COM+W1+W2)", "stored.b.append(b)")
cmd.alter("resn COM", "vdw=stored.b[0]")
cmd.alter("resn W1", "vdw=stored.b[1]")
cmd.alter("resn W2", "vdw=stored.b[2]")
cmd.rebuild()

cmd.spectrum("b", "blue_white_red", "O2.pywindow", [0, max(stored.b)])
cmd.ramp_new("radius", "O2.pywindow", [0, max(stored.b)], ["blue", "white", "red"])
cmd.orient()

cmd.save("O2.pse")
