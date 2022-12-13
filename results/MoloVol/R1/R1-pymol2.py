import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/R1.pdb", quiet=False)
cmd.load("R1_surface-map_grid-0.6_rad1-1.4_rad2-7.0_full-structure.dx", quiet=False)

cmd.isomesh("cavities", "R1_surface-map_grid-0.6_rad1-1.4_rad2-7.0_full-structure", level=3.0) # Isolated cavitycmd.orient()

cmd.save("R1.pse")