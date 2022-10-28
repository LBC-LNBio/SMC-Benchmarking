import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/F1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/F1/F1_frameInfo/F1_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/F1/F1_frameInfo/F1_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/F1/F1_frameInfo/F1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "F1_volumetric_density")

cmd.alter("obj F1_inclusion", "vdw=0.5")
cmd.alter("obj F1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("F1.pse")
