import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/O2.pdb", quiet=False)
cmd.load("O2_surface-map_grid-0.6_rad1-1.4_rad2-15.0_full-structure.dx", quiet=False)

cmd.isomesh("cavities", "O2_surface-map_grid-0.6_rad1-1.4_rad2-15.0_full-structure", level=3.0) # Isolated cavitycmd.orient()

cmd.save("O2.pse")
