import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B3.pdb", quiet=False)
cmd.load("B3_frameInfo/B3_inclusion.pdb", quiet=False)
cmd.load("B3_frameInfo/B3_volumetric_density.dx", quiet=False)
cmd.load("B3_frameInfo/B3_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B3_volumetric_density")

cmd.alter("obj B3_inclusion", "vdw=0.5")
cmd.alter("obj B3_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B3.pse")
