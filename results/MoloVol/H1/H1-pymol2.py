import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/H1.pdb", quiet=False)
cmd.load("H1_surface-map_grid-0.6_rad1-1.4_rad2-5.0_full-structure.dx", quiet=False)

cmd.isomesh("cavities", "H1_surface-map_grid-0.6_rad1-1.4_rad2-5.0_full-structure", level=3.0) # Isolated cavitycmd.orient()

cmd.save("H1.pse")
