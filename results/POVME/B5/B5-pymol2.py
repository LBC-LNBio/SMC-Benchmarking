import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B5.pdb", quiet=False)
cmd.load("B5_frameInfo/B5_inclusion.pdb", quiet=False)
cmd.load("B5_frameInfo/B5_volumetric_density.dx", quiet=False)
cmd.load("B5_frameInfo/B5_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B5_volumetric_density")

cmd.alter("obj B5_inclusion", "vdw=0.25")
cmd.alter("obj B5_frame_1_surface", "vdw=0.25")
cmd.rebuild()

cmd.orient()

cmd.save("B5.pse")
