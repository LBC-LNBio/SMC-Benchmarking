import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/O1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/O1/O1_frameInfo/O1_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/O1/O1_frameInfo/O1_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/O1/O1_frameInfo/O1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "O1_volumetric_density")

cmd.alter("obj O1_inclusion", "vdw=0.5")
cmd.alter("obj O1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()