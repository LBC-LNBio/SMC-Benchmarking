import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/O2.pdb", quiet=False)
cmd.load("O2_frameInfo/O2_inclusion.pdb", quiet=False)
cmd.load("O2_frameInfo/O2_volumetric_density.dx", quiet=False)
cmd.load("O2_frameInfo/O2_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "O2_volumetric_density")

cmd.alter("obj O2_inclusion", "vdw=0.25")
cmd.alter("obj O2_frame_1_surface", "vdw=0.25")
cmd.rebuild()

cmd.orient()

cmd.save("O2.pse")
