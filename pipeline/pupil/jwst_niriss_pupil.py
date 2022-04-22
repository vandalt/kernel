# %%
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from xara import kpi
from xara.core import create_discrete_model, symetrizes_model

import webbpsf

# %%
# Parameters.
step = 0.2  # m, grid step size
# tmin = 1e-2  # minimum transmission for a sub-aperture to be considered
tmin = 0.8  # minimum transmission for a sub-aperture to be considered
# tmin = 5e-1  # minimum transmission for a sub-aperture to be considered
binary = True  # binary or grey mask
textpath = "niriss_clear_pupil.txt"
fitspath = "niriss_clear_pupil.fits"
bmax = None

# %%
pupil_path = Path(os.getenv("WEBBPSF_PATH")) / "NIRISS/optics/MASKCLEAR.fits.gz"
hdul = fits.open(pupil_path)
apert_clear = hdul["PRIMARY"].data
ppscale = hdul["PRIMARY"].header["REALSCAL"]  # pupil pixel size in meters

PSZ = apert_clear.shape[0]  # pupil image size

# %%
clearp = webbpsf.optics.NIRISS_CLEARP()
# Rescale clearp to match pix scale
clearp_pup = np.round(clearp.sample(npix=1024, grid_size=1024 * ppscale))

# %%
aper = apert_clear * np.flipud(clearp_pup)

# %%
# Create discrete pupil model using XARA.
model = create_discrete_model(aper, ppscale, step, binary=binary, tmin=tmin)
model = symetrizes_model(model)
np.savetxt(textpath, model, fmt="%+.6e %+.6e %.2f")
mykpi = kpi.KPI(fname=textpath, bmax=bmax)
mykpi.package_as_fits(fname=fitspath)

# Plot pupil model.
# f = mykpi.plot_pupil_and_uv(cmap="inferno")
# f = mykpi.plot_pupil_and_uv()
mykpi.plot_pupil_and_uv(
    xymax=3.5, cmap=plt.cm.plasma_r, ssize=9, figsize=(10, 5), marker="o"
)
plt.savefig(textpath[:-3] + "pdf")
plt.show(block=True)
plt.close()
