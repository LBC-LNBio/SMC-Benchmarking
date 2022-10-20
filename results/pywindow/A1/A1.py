import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../data/A1.pdb")
cmd.load("A1.pywindow.pdb")

stored.b = []
cmd.iterate("(resn COM+W1+W1+W2+W3+W4+W5+W6+W7+W8)", "stored.b.append(b)")
cmd.alter("resn COM", "vdw=stored.b[0]")
cmd.alter("resn W1", "vdw=stored.b[1]")
cmd.alter("resn W2", "vdw=stored.b[2]")
cmd.alter("resn W3", "vdw=stored.b[3]")
cmd.alter("resn W4", "vdw=stored.b[4]")
cmd.alter("resn W5", "vdw=stored.b[5]")
cmd.alter("resn W6", "vdw=stored.b[6]")
cmd.alter("resn W7", "vdw=stored.b[7]")
cmd.alter("resn W8", "vdw=stored.b[8]")

cmd.rebuild()

cmd.spectrum("b", "blue_white_red", "A1.pywindow", [0, max(stored.b)])
cmd.ramp_new("radius", "A1.pywindow", [0, max(stored.b)], ["blue", "white", "red"])
