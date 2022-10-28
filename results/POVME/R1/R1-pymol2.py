import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/R1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/R1/R1_frameInfo/R1_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/R1/R1_frameInfo/R1_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/R1/R1_frameInfo/R1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "R1_volumetric_density")

cmd.alter("obj R1_inclusion", "vdw=0.5")
cmd.alter("obj R1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("R1.pse")
