import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/A1.pdb", quiet=False)
cmd.load("A1_surface-map_grid-0.6_rad1-1.4_rad2-6.0_full-structure.dx", quiet=False)

cmd.isomesh("cavities", "A1_surface-map_grid-0.6_rad1-1.4_rad2-6.0_full-structure", level=3.0) # Isolated cavity
cmd.orient()

cmd.save("A1.pse")
