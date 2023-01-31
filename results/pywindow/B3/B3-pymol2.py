import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B3.pdb", quiet=False)
cmd.load("B3.pywindow.pdb", quiet=False)

stored.b = []
cmd.iterate("(resn COM+W1+W2+W3+W4)", "stored.b.append(b)")
cmd.alter("resn COM", "vdw=stored.b[0]")
cmd.alter("resn W1", "vdw=stored.b[1]")
cmd.alter("resn W2", "vdw=stored.b[2]")
cmd.alter("resn W3", "vdw=stored.b[3]")
cmd.alter("resn W4", "vdw=stored.b[4]")
cmd.rebuild()

cmd.spectrum("b", "blue_white_red", "B3.pywindow", 0, max(stored.b))
cmd.ramp_new("radius", "B3.pywindow", [0, max(stored.b)], ["blue", "white", "red"])
cmd.orient()

cmd.save("B3.pse")
