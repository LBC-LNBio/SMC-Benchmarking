import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/N1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/N1/N1_frameInfo/N1_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/N1/N1_frameInfo/N1_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/N1/N1_frameInfo/N1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "N1_volumetric_density")

cmd.alter("obj N1_inclusion", "vdw=0.5")
cmd.alter("obj N1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("N1.pse")
