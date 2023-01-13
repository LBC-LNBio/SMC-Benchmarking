import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B4.pdb", quiet=False)
cmd.load("B4_frameInfo/B4_inclusion.pdb", quiet=False)
cmd.load("B4_frameInfo/B4_volumetric_density.dx", quiet=False)
cmd.load("B4_frameInfo/B4_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B4_volumetric_density")

cmd.alter("obj B4_inclusion", "vdw=0.5")
cmd.alter("obj B4_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B4.pse")
