import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B11.pdb", quiet=False)
cmd.load("B11.pywindow.pdb", quiet=False)

stored.b = []
cmd.iterate("(resn COM)", "stored.b.append(b)")
cmd.alter("resn COM", "vdw=stored.b[0]")
cmd.rebuild()

cmd.spectrum("b", "blue_white_red", "B11.pywindow", 0, max(stored.b))
cmd.ramp_new("radius", "B11.pywindow", [0, max(stored.b)], ["blue", "white", "red"])
cmd.orient()

cmd.save("B11.pse")
