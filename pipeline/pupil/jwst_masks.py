# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import poppy
import webbpsf
from astropy.io import fits
from webbpsf import constants

# %%
pupil_path = Path(os.getenv("WEBBPSF_PATH")) / "NIRISS/optics/MASKCLEAR.fits.gz"
# pupil_path = Path("/home/vandal/astro/xara/jwst/MASK_CLEAR.fits")
hdul = fits.open(pupil_path)
apert_clear = hdul["PRIMARY"].data
hdr = hdul["PRIMARY"].header
ppscale = hdr["REALSCAL"]  # pupil pixel size in meters

PSZ = apert_clear.shape[0]  # pupil image size

# %%
clearp = webbpsf.optics.NIRISS_CLEARP()
# Rescale clearp to match pix scale
clearp_pup = np.round(clearp.sample(npix=1024, grid_size=1024 * ppscale))

# %%
aper_from_clear = apert_clear * np.flipud(clearp_pup)

# %%
plt.imshow(aper_from_clear, origin="lower")
plt.show(block=True)

# %%
hdr_clearp = fits.getheader("/home/vandal/astro/xara/jwst/MASK_CLEARP.fits")
aper_from_clearp = fits.getdata("/home/vandal/astro/xara/jwst/MASK_CLEARP.fits")

# %%
plt.imshow(aper_from_clearp, origin="lower")
plt.show(block=True)

# %%
diff = aper_from_clear - np.fliplr(aper_from_clearp)
# diff = aper_from_clear - aper_from_clearp
print(np.std(diff))
plt.imshow(diff, origin="lower")
plt.show(block=True)

# %%
primary = webbpsf.optics.WebbPrimaryAperture()
primary_diam = constants.JWST_CIRCUMSCRIBED_DIAMETER
wheelrad = 39.0  # mm
pupil_mag = primary_diam / wheelrad  # 39 mm corresponds to circumscribed diam

# make a flat pupil wavefront
mag = 1  # Careful, this eats memory quickly
w = poppy.Wavefront(npix=1024 * mag, diam=primary_diam)

primask = primary.get_transmission(w)

aper_from_prim = primask * np.flipud(clearp_pup)


# %%
pupil_path_webbpsf = Path(os.getenv("WEBBPSF_PATH")) / "NIRISS/optics/MASKCLEAR.fits.gz"
pupil_path_calwebb = Path("/home/vandal/astro/xara/jwst/MASK_CLEAR.fits")
apert_webbpsf = fits.getdata(pupil_path_webbpsf)
apert_calwebb = fits.getdata(pupil_path_calwebb)
diff = apert_webbpsf - apert_calwebb
plt.imshow(diff, origin="lower")
plt.show(block=True)


# %%
# The masks from Xara and NIRISS MASKCLEAR are the same
pupil_path_webbpsf = Path(os.getenv("WEBBPSF_PATH")) / "NIRISS/optics/MASKCLEAR.fits.gz"
pupil_path_calwebb = Path("/home/vandal/astro/xara/jwst/MASK_CLEAR.fits")
apert_webbpsf = fits.getdata(pupil_path_webbpsf)
apert_calwebb = fits.getdata(pupil_path_calwebb)
diff = apert_webbpsf - apert_calwebb
plt.imshow(diff, origin="lower")
plt.show(block=True)

# %%
apert_analytical = primask
diff = apert_webbpsf - apert_analytical
plt.imshow(diff, origin="lower")
plt.show(block=True)


# %%
# Compare AMI-work and calwebb CLEARs
pupil_path_amiwork = Path("/home/vandal/astro/ami-work/MASK_CLEAR_mag10.fits")
pupil_path_calwebb = Path("/home/vandal/astro/xara/jwst/MASK_CLEAR.fits")
pupil_path_calwebb_cp = Path("/home/vandal/astro/xara/jwst/MASK_CLEARP.fits")
apert_amiwork = fits.getdata(pupil_path_amiwork)
apert_calwebb = fits.getdata(pupil_path_calwebb)
apert_calwebb_cp = fits.getdata(pupil_path_calwebb_cp)
diff = apert_amiwork - apert_calwebb
diff_cp = apert_amiwork - apert_calwebb_cp

# %%
plt.subplot(1, 2, 1)
plt.imshow(diff, origin="lower")
plt.title("New CLEAR (analtytical) vs current (MASKCLEAR)")
plt.subplot(1, 2, 2)
plt.title("New CLEAR (analtytical) vs current CLEARP")
plt.imshow(diff_cp, origin="lower")
plt.show(block=True)

# %%
# apert_analytical = primask
diff = apert_amiwork - apert_analytical
plt.imshow(diff, origin="lower")
plt.show(block=True)
