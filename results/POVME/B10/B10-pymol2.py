import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B10.pdb", quiet=False)
cmd.load("B10_frameInfo/B10_inclusion.pdb", quiet=False)
cmd.load("B10_frameInfo/B10_volumetric_density.dx", quiet=False)
cmd.load("B10_frameInfo/B10_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B10_volumetric_density")

cmd.alter("obj B10_inclusion", "vdw=0.5")
cmd.alter("obj B10_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B10.pse")
