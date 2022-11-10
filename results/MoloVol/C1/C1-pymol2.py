import pymol
from pymol import cmd, stored

pymol.finish_launching(["pymol", "-q"])

cmd.load("../../../hosts/C1.pdb", quiet=False)
cmd.load("C1_surface-map_grid-0.6_rad1-1.4_rad2-5.0_full-structure.dx", quiet=False)

cmd.isomesh("cavities", "C1_surface-map_grid-0.6_rad1-1.4_rad2-5.0_full-structure", level=1.5) # Molecular surfacecmd.orient()

cmd.save("C1.pse")
