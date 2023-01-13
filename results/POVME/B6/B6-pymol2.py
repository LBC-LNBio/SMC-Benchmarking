import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B6.pdb", quiet=False)
cmd.load("B6_frameInfo/B6_inclusion.pdb", quiet=False)
cmd.load("B6_frameInfo/B6_volumetric_density.dx", quiet=False)
cmd.load("B6_frameInfo/B6_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B6_volumetric_density")

cmd.alter("obj B6_inclusion", "vdw=0.25")
cmd.alter("obj B6_frame_1_surface", "vdw=0.25")
cmd.rebuild()

cmd.orient()

cmd.save("B6.pse")
