# %% [markdown]
# _Trying_ do generate NIRISS pupil with webbpsf for kernel phase
#
# First, use NIRISS instrument model, wavefront with CLEARP mask
#
# _NOTE: CLEARP is **not** the same as CLEAR (i.e. the WEBBPSF file)._

# %%
# Load NIRISS instrument model
import webbpsf
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.ndimage import rotate

niriss = webbpsf.NIRISS()

# Set F430M filter so that auto_pupil uses CLEARP
niriss.filter = "F430M"

# %% [markdown]
# Now, use optical elements separately

# %%
import os
from pathlib import Path
# optsys = niriss.get_optical_system()
pupil_path = Path(os.getenv("WEBBPSF_PATH")) / "NIRISS/optics/MASKCLEAR.fits.gz"
hdr = fits.getheader(pupil_path)
clear = fits.getdata(pupil_path)
pscale = hdr["PUPLSCAL"]  # ???: Or REALSCAL


# %%
# clear_pup = np.rot90(np.round(clear), k=2)
clear_pup = np.flipud(np.round(clear))

# %%
plt.imshow(clear_pup, origin="lower")
plt.show()

# %%
clearp = webbpsf.optics.NIRISS_CLEARP()
# Rescale clearp to match pix scale
clearp_pup = np.round(clearp.sample(npix=1024, grid_size=1024 * pscale))

data = clearp_pup * clear_pup

# %%
plt.imshow(data, origin="lower")
plt.show()

# %%
hdul = fits.PrimaryHDU(data=data, header=hdr)
hdul.writeto("CLEARP_CUSTOM.fits", overwrite=True)

# %%
plt.imshow(data - wdata, origin="lower")
plt.show()


# %%
import xaosim

pupil = xaosim.pupil.JWST(1024, pscale=0.0064486953125)


# %%
plt.subplot(1, 3, 1)
plt.imshow(data, origin="lower")
plt.title("WebbPSF")
plt.subplot(1, 3, 2)
plt.imshow(pupil)
plt.title("Current xaosim")
plt.subplot(1, 3, 3)
plt.imshow(pupil_new)
plt.title("This PR")
plt.show()
