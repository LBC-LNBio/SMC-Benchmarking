import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/data/A1.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/A1/A1_frameInfo/A1_inclusion.pdb", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/A1/A1_frameInfo/A1_volumetric_density.dx", quiet=False)
cmd.load("/home/jvsguerra/remote-repos/SMC-Benchmarking/results/POVME/A1/A1_frameInfo/A1_frame_1_surface.pdb", quiet=False)

cmd.show("dots", "A1_volumetric_density")

cmd.alter("obj A1_inclusion", "vdw=0.5")
cmd.alter("obj A1_frame_1_surface", "vdw=0.5")
cmd.rebuild()

cmd.orient()

cmd.save("A1.pse")
