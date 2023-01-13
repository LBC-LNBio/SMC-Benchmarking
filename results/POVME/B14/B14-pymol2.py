import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B14.pdb", quiet=False)
cmd.load("B14_frameInfo/B14_inclusion.pdb", quiet=False)
cmd.load("B14_frameInfo/B14_volumetric_density.dx", quiet=False)
cmd.load("B14_frameInfo/B14_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B14_volumetric_density")

cmd.alter("obj B14_inclusion", "vdw=0.5")
cmd.alter("obj B14_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B14.pse")
