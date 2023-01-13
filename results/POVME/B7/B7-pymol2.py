import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B7.pdb", quiet=False)
cmd.load("B7_frameInfo/B7_inclusion.pdb", quiet=False)
cmd.load("B7_frameInfo/B7_volumetric_density.dx", quiet=False)
cmd.load("B7_frameInfo/B7_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B7_volumetric_density")

cmd.alter("obj B7_inclusion", "vdw=0.5")
cmd.alter("obj B7_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B7.pse")
