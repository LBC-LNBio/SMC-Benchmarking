import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/C1.pdb", quiet=False)
cmd.load("C1_frameInfo/C1_inclusion.pdb", quiet=False)
cmd.load("C1_frameInfo/C1_volumetric_density.dx", quiet=False)
cmd.load("C1_frameInfo/C1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "C1_volumetric_density")

cmd.alter("obj C1_inclusion", "vdw=0.5")
cmd.alter("obj C1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("C1.pse")
