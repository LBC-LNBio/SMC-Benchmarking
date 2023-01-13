import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/N1.pdb", quiet=False)
cmd.load("N1.pywindow.pdb", quiet=False)

stored.b = []
cmd.iterate("(resn COM+W1+W2+W3+W4+W5)", "stored.b.append(b)")
cmd.alter("resn COM", "vdw=stored.b[0]")
cmd.alter("resn W1", "vdw=stored.b[1]")
cmd.alter("resn W2", "vdw=stored.b[2]")
cmd.alter("resn W3", "vdw=stored.b[3]")
cmd.alter("resn W4", "vdw=stored.b[4]")
cmd.alter("resn W5", "vdw=stored.b[5]")
cmd.rebuild()

cmd.spectrum("b", "blue_white_red", "N1.pywindow", [0, max(stored.b)])
cmd.ramp_new("radius", "N1.pywindow", [0, max(stored.b)], ["blue", "white", "red"])
cmd.orient()

cmd.save("N1.pse")
