"""
Computes vorticity from given U and V fluid velocity fields and outputs on the same grid as the
given temperature field. Intended for output from the GYRE configuration.

Requires a mesh_mask.nc file, which can be output from NEMO by setting the `ln_meshmask` namelist
parameter in `&namdom` to `.true.`.
"""

from cdfcurl import cdfcurl
from iris import load_cube, save
from argparse import ArgumentParser
from cf_units import Unit

parser = ArgumentParser(description="Computes vorticity from given U and V fluid velocity fields" +
                        " and outputs on the same grid as the given temperature field. Intended" +
                        " for output from the GYRE configuration.")
parser.add_argument("u_file", type=str, help="The name of the file containing the U field")
parser.add_argument("v_file", type=str, help="The name of the file containing the V field")
parser.add_argument("t_file", type=str,
                    help="The name of the file containing the temperature field (copied as a" +
                    " template for the vorticity file)")
parser.add_argument("mesh_mask_file", type=str,
                    help="The name of the mesh mask file")
parser.add_argument("output_filename", type=str, help="The output file")
parser.add_argument("--u_varname", dest="u_varname", type=str, default="vozocrtx",
                    help="The varname of the U field (default: vozocrtx)")
parser.add_argument("--v_varname", dest="v_varname", type=str, default="vomecrty",
                    help="The varname of the V field (default: vomecrty)")
parser.add_argument("--u_stanname", dest="u_stanname", type=str, default="sea_water_x_velocity",
                    help="The standard name of the U field (default: sea_water_x_velocity)")
parser.add_argument("--v_stanname", dest="v_stanname", type=str, default="sea_water_y_velocity",
                    help="The standard name of the V field (default: sea_water_y_velocity)")
parser.add_argument("--t_stanname", dest="t_stanname", type=str, default="sea_surface_temperature",
                    help="The standard name of the V field (default: sea_surface_temperature)")
args = parser.parse_args()

# Load temperature cube and copy it
t_cube = load_cube(args.t_file, args.t_stanname)
vor_cube = t_cube[:, ...].copy()
vor_cube.rename("vorticity")
vor_cube.units = Unit("s**-1")

# Compute vorticity for each time step in the output cube
for t in range(t_cube.shape[0]):
    print(f"Processing timestep {t+1} of {t_cube.shape[0]}")
    _, _, vor_cube.data[t, :, :], _ = cdfcurl(args.u_file, args.u_varname,
                                              args.v_file, args.v_varname,
                                              args.mesh_mask_file,
                                              kt=t,
                                              loverf=True)

save(vor_cube, args.output_filename)
