import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/B1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/B1/B1_frameInfo/B1_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/B1/B1_frameInfo/B1_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/B1/B1_frameInfo/B1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B1_volumetric_density")

cmd.alter("obj B1_inclusion", "vdw=0.5")
cmd.alter("obj B1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()
