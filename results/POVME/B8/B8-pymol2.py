import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B8.pdb", quiet=False)
cmd.load("B8_frameInfo/B8_inclusion.pdb", quiet=False)
cmd.load("B8_frameInfo/B8_volumetric_density.dx", quiet=False)
cmd.load("B8_frameInfo/B8_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B8_volumetric_density")

cmd.alter("obj B8_inclusion", "vdw=0.5")
cmd.alter("obj B8_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B8.pse")