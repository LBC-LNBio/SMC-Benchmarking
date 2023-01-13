import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B9.pdb", quiet=False)
cmd.load("B9_frameInfo/B9_inclusion.pdb", quiet=False)
cmd.load("B9_frameInfo/B9_volumetric_density.dx", quiet=False)
cmd.load("B9_frameInfo/B9_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B9_volumetric_density")

cmd.alter("obj B9_inclusion", "vdw=0.5")
cmd.alter("obj B9_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B9.pse")
