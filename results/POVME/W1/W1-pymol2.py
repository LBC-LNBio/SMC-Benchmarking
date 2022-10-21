import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/W1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/W1/W1_frameInfo/W1_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/W1/W1_frameInfo/W1_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/W1/W1_frameInfo/W1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "W1_volumetric_density")

cmd.alter("obj W1_inclusion", "vdw=0.5")
cmd.alter("obj W1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()
