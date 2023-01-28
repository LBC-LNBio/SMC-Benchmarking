import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/B12.pdb", quiet=False)
cmd.load("B12_surface-map_grid-0.6_rad1-1.4_rad2-5.0_full-structure.dx", quiet=False)

cmd.isomesh("cavities", "B12_surface-map_grid-0.6_rad1-1.4_rad2-5.0_full-structure", level=3.0) # Isolated cavity
cmd.orient()

cmd.save("B12.pse")
