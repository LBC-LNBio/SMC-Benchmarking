import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B13.pdb", quiet=False)
cmd.load("B13_frameInfo/B13_inclusion.pdb", quiet=False)
cmd.load("B13_frameInfo/B13_volumetric_density.dx", quiet=False)
cmd.load("B13_frameInfo/B13_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B13_volumetric_density")

cmd.alter("obj B13_inclusion", "vdw=0.5")
cmd.alter("obj B13_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B13.pse")
