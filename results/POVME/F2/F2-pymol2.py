import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/F2.pdb", quiet=False)
cmd.load("F2_frameInfo/F2_inclusion.pdb", quiet=False)
cmd.load("F2_frameInfo/F2_volumetric_density.dx", quiet=False)
cmd.load("F2_frameInfo/F2_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "F2_volumetric_density")

cmd.alter("obj F2_inclusion", "vdw=1.0")
cmd.alter("obj F2_frame_1_surface", "vdw=1.0")
cmd.rebuild()

cmd.orient()

cmd.save("F2.pse")
