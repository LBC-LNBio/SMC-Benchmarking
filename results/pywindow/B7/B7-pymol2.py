import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B7.pdb", quiet=False)
cmd.load("B7.pywindow.pdb", quiet=False)

stored.b = []
cmd.iterate("(resn COM+W1)", "stored.b.append(b)")
cmd.alter("resn COM", "vdw=stored.b[0]")
cmd.alter("resn W1", "vdw=stored.b[1]")
cmd.rebuild()

cmd.spectrum("b", "blue_white_red", "B7.pywindow", [0, max(stored.b)])
cmd.ramp_new("radius", "B7.pywindow", [0, max(stored.b)], ["blue", "white", "red"])
cmd.orient()

cmd.save("B7.pse")
