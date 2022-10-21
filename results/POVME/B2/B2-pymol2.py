import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/data/B2.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/B2/B2_frameInfo/B2_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/B2/B2_frameInfo/B2_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/moc-benchmarking/results/POVME/B2/B2_frameInfo/B2_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B2_volumetric_density")

cmd.alter("obj B2_inclusion", "vdw=0.5")
cmd.alter("obj B2_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()
