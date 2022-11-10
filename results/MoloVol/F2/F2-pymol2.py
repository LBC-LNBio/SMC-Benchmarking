import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/F2.pdb", quiet=False)
cmd.load("F2_surface-map_grid-0.6_rad1-1.4_rad2-20.0_full-structure.dx", quiet=False)

cmd.isomesh("cavities", "F2_surface-map_grid-0.6_rad1-1.4_rad2-20.0_full-structure", level=1.5) # Molecular surfacecmd.orient()

cmd.save("F2.pse")
