import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/B3.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/B3/B3_frameInfo/B3_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/B3/B3_frameInfo/B3_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/B3/B3_frameInfo/B3_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "B3_volumetric_density")

cmd.alter("obj B3_inclusion", "vdw=0.5")
cmd.alter("obj B3_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("B3.pse")
