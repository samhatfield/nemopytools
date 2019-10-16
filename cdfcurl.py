from numpy import zeros, pi, sin
from netCDF4 import Dataset


def cdfcurl(u_file, u_var, v_file, v_var, mesh_mask_file, **kwargs):
    """
    Note: stolen from Julian Mak at https://github.com/julianmak/NEMO-related.

    Compute the curl on something F-points from U and V variables
    EDIT: unlike cdfcurl.f90 I did not put an option to work out
        things like wind-stress curl for A-grid data
        (e.g. CORE-II forcing, where both components are on the same point)
    Inputs:
    u_file         = string for file with u in
    u_var          = string for u variable name
    v_file         = string for file with v in
    v_var          = string for v variable name
    mesh_mask_file = string for mesh mask file
    Optional arguments (as of 10 Apr 2018):
    kt       = number for using a specified time entry
    kz       = vertical level to use
    lprint   = True   for printing out variable names in netcdf file
    lperio   = True   for imposing E-W periodicity
    loverf   = True   for having curl / f
    Returns:
    lonF (glamf), latF (gphif), curl for plotting, opt_dic for record
    """

    # Some defaults for optional arguments
    opt_dic = {"kt": 0,
               "kz": 0,
               "lprint": False,
               "lperio": False,
               "loverf": False}

    # Overwrite the options by cycling through the input dictionary
    for key in kwargs:
        opt_dic[key] = kwargs[key]

    # Open some files and pull variables out
    cn_mask = Dataset(mesh_mask_file)
    npiglo = len(cn_mask.dimensions["x"])
    npjglo = len(cn_mask.dimensions["y"])
    glamf = cn_mask.variables["glamf"][0, :, :]
    gphif = cn_mask.variables["gphif"][0, :, :]
    e1u = cn_mask.variables["e1u"][0, :, :]
    e1f = cn_mask.variables["e1f"][0, :, :]
    e2v = cn_mask.variables["e2v"][0, :, :]
    e2f = cn_mask.variables["e2f"][0, :, :]
    fmask = cn_mask.variables["fmask"][0, opt_dic["kz"], :, :]
    cn_mask.close()

    cf_ufil = Dataset(u_file)
    if opt_dic["lprint"]:
        print(cf_ufil)

    # Check the dimension of the u field and load accordingly
    depth_var = False
    dimen = cf_ufil.variables[u_var].dimensions
    for key in dimen:
        if "depth" in key:
            print("a depth variable detected in u_file, loading with kz option")
            depth_var = True
            break
    if depth_var:
        print("a depth variable detected in u_file, loading with kz option")
        zu = cf_ufil.variables[u_var][opt_dic["kt"], opt_dic["kz"], :, :]
    else:
        print("no depth variable detected in u_file, loading only with kt option")
        zu = cf_ufil.variables[u_var][opt_dic["kt"], :, :]
    cf_ufil.close()

    cf_vfil = Dataset(v_file)
    if opt_dic["lprint"]:
        print(cf_vfil)
    # Check the dimension of the v field and load accordingly
    depth_var = False
    dimen = cf_vfil.variables[v_var].dimensions
    for key in dimen:
        if "depth" in key:
            print("a depth variable detected in v_file, loading with kz option")
            depth_var = True
            break
    if depth_var:
        print("a depth variable detected in v_file, loading with kz option")
        zv = cf_vfil.variables[v_var][opt_dic["kt"], opt_dic["kz"], :, :]
    else:
        print("no depth variable detected in v_file, loading only with kt option")
        zv = cf_vfil.variables[v_var][opt_dic["kt"], :, :]
    cf_vfil.close()

    # Check periodicity
    if (glamf[0, 0] - glamf[0, npiglo-2] == 0):
        print("E-W periodicity detected!")
    if not opt_dic["lperio"]:
        print("--> forcing lperio = True")
        opt_dic["lperio"] = True

    curl = zeros((npjglo, npiglo))

    for ji in range(npiglo-1):
        for jj in range(npjglo-1):
            curl[jj, ji] = (  e2v[jj  , ji+1]*zv[jj  , ji+1]
                            - e2v[jj  , ji  ]*zv[jj  , ji  ]
                            - e1u[jj+1, ji  ]*zu[jj+1, ji  ]
                            + e1u[jj  , ji  ]*zu[jj  , ji  ])\
                            * fmask[jj, ji] / (e1f[jj, ji] * e2f[jj, ji])

    if opt_dic["lperio"]:
        curl[:, npiglo-1] = curl[:, 1]

    if opt_dic["loverf"]:
        dl_Ω = (2 * pi / 86400.0)
        dl_ff = 2 * dl_Ω * sin(gphif * pi / 180.0)
        # Avoid divide by zero
        curl[(dl_ff == 0)] = 0.0
        dl_ff[(dl_ff == 0)] = 1.0
        curl = curl / dl_ff

    return glamf, gphif, curl, opt_dic
