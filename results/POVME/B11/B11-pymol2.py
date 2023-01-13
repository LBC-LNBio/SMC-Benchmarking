import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B11.pdb", quiet=False)
cmd.load("B11_frameInfo/B11_inclusion.pdb", quiet=False)
cmd.load("B11_frameInfo/B11_volumetric_density.dx", quiet=False)
cmd.load("B11_frameInfo/B11_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B11_volumetric_density")

cmd.alter("obj B11_inclusion", "vdw=0.5")
cmd.alter("obj B11_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B11.pse")
