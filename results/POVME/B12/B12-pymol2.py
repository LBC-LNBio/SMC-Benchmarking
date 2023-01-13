import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B12.pdb", quiet=False)
cmd.load("B12_frameInfo/B12_inclusion.pdb", quiet=False)
cmd.load("B12_frameInfo/B12_volumetric_density.dx", quiet=False)
cmd.load("B12_frameInfo/B12_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B12_volumetric_density")

cmd.alter("obj B12_inclusion", "vdw=0.5")
cmd.alter("obj B12_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B12.pse")
