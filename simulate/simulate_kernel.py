# %%
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import driver_scene
import matplotlib.pyplot as plt
import numpy as np
import pyami.simcode.utils as U
import webbpsf
from astropy.io import fits


# %%
def get_psf(
    filt: str,
    fov_pixels: int,
    oversample: int,
    filename: Optional[Union[str, Path]] = None,
    pupil_mask: str = "CLEARP",
    source: Optional[str] = "A0V",
    display: bool = False,
):
    """
    Recompute the PSF using webbpsf

    :param _filter: str, the filter to use
    :param filename: str, the output filename
    :param fov_pixels: int, the fov in pixels
    :param oversample: int, the oversampling factor
    :param pupil_mask: str, the pupil mask to use

    :return: None
    """
    # get niriss instance from webb psf
    niriss = webbpsf.NIRISS()
    # set the filter name
    niriss.filter = filt
    # set the pupil mask
    niriss.pupil_mask = pupil_mask

    niriss.pixelscale = U.pixscl

    source = webbpsf.specFromSpectralType(source)

    # run the psf calculation
    psf_fits = niriss.calc_psf(
        outfile=str(filename),
        fov_pixels=fov_pixels,
        oversample=oversample,
        display=display,
        source=source,
    )

    return psf_fits


# %%
def ami_sim_sky(
    fov_pixels: float,
    oversample: int,
    pix_scale: float,
    plot: bool = False,
) -> Tuple[np.ndarray, Dict]:
    """
    Create a simple target scene (one point source at the center)

    :param fov_pixels: the native image size (FOV in pixels)
    :param oversample: the oversampling factor
    :param pix_scale: pixel scale (expected to be very close to 0.0656
                      arcsec/pixel)
    :param ext_flux: extracted flux in e-/s
    :param tot_exp: total exposure time in seconds

    :return: image and header
    """
    # work out the over sample pixel width
    osample_pix_width = fov_pixels * oversample
    # create an image full of zeros
    image = np.zeros([osample_pix_width, osample_pix_width])
    # get the central pixel
    xcen = image.shape[1] // 2
    ycen = image.shape[0] // 2
    # -------------------------------------------------------------------------
    # Add primary
    # -------------------------------------------------------------------------
    # add primary
    # sky data is normalized to 1 and mult by count rate in ami_sim
    COUNT0 = 1.0
    image[xcen, ycen] = COUNT0

    if plot:
        plt.imshow(image)
        plt.show(block=True)
        plt.close()

    # -------------------------------------------------------------------------
    # Write sky file
    # -------------------------------------------------------------------------
    hdict = dict()
    hdict["FOVPX"] = (fov_pixels, "Input FOV pixels")
    hdict["OSAMPLE"] = (oversample, "Input Oversample")
    hdict["PIXSCALE"] = (pix_scale, "Input pixel scale")
    hdict["COUNT0"] = (COUNT0, "Count rate of primary")

    return image, hdict


def sky_add_companion(
    image: np.ndarray,
    hdict: Dict,
    num: int,
    position_angle: float,
    separation: float,
    contrast: float,
    plot: bool = False,
) -> Tuple[np.ndarray, Dict]:
    """
    Add a companion to an image where the primary is assumed to be at the
    center of the image

    :param params: ParamDict, parameter dictionary of constants
    :param image: numpy array (2D) - the image before this companion is added
    :param hdict: ParamDict, the keys for header
    :param num: int, the companion number (must be unique for target)
    :param position_angle: float, the position angle of companion from
                           primary (ra-dec coord system) in degrees
    :param separation: float, separation between primary (assume to be at the
                       center or the detector) and companion in arcsec
    :param contrast: contrast of companion
    :param plot: bool, if True plots

    :return: updated  image and header
    """
    # -------------------------------------------------------------------------
    # set up
    # -------------------------------------------------------------------------
    # get parameters
    pix_scale = float(hdict["PIXSCALE"][0])
    oversample = float(hdict["OSAMPLE"][0])
    count_rate = float(hdict["COUNT0"][0])
    # get the central pixel
    xcen = image.shape[1] // 2
    ycen = image.shape[0] // 2

    # -------------------------------------------------------------------------
    # Calculate pixel position for companion
    # -------------------------------------------------------------------------
    # get the dx and dy in pixels
    dx = np.cos(np.pi * position_angle / 180) * (separation / pix_scale)
    dy = np.sin(np.pi * position_angle / 180) * (separation / pix_scale)
    # get oversampling
    odx = dx * oversample
    ody = dy * oversample
    # get the x and y center for companion
    xcen2 = int(xcen + odx)
    ycen2 = int(ycen + ody)
    # -------------------------------------------------------------------------
    # Add companion to image
    # -------------------------------------------------------------------------
    # test bounds
    image[xcen2, ycen2] = count_rate * contrast

    if plot:
        plt.imshow(image)
        plt.show(block=True)
        plt.close()

    # -------------------------------------------------------------------------
    # Work out true separation and angle (given the rounding to the nearest
    #   pixel)
    # -------------------------------------------------------------------------
    # true separation
    tx2 = (xcen - int(xcen2)) ** 2
    ty2 = (ycen - int(ycen2)) ** 2
    true_sep = np.sqrt(tx2 + ty2) * pix_scale / oversample
    # true position angle
    ay2 = ycen - ycen2
    true_angle = 180 * np.arcsin(ay2 / np.sqrt(tx2 + ty2)) / np.pi + 180

    # -------------------------------------------------------------------------
    # add to hdict
    # -------------------------------------------------------------------------
    # text for comment
    ctxt = "companion {0}".format(num)
    # save input separation
    kw_in_sep = "IN_SEP{0}".format(num)
    hdict[kw_in_sep] = (
        separation,
        "Input separation of scene [arcsec] {0}".format(ctxt),
    )
    # save input position angle
    kw_in_ang = "IN_ANG{0}".format(num)
    hdict[kw_in_ang] = (
        position_angle,
        "Input position angle of scene [deg] {0}".format(ctxt),
    )
    # save true separtion
    kw_t_sep = "T_SEP{0}".format(num)
    hdict[kw_t_sep] = (true_sep, "True separation of scene [arcsec] {0}".format(ctxt))
    # save true position angle
    kw_t_ang = "T_ANG{0}".format(num)
    hdict[kw_t_ang] = (
        true_angle,
        "True position angle of scene [deg] {0}".format(ctxt),
    )
    # save count rate of companion
    kw_count = "COUNT{0}".format(num)
    hdict[kw_count] = (count_rate * contrast, "Count rate of {0}".format(ctxt))

    # -------------------------------------------------------------------------
    # return image and header
    return image, hdict


def ami_sim_save_scene(outpath: str, image: np.ndarray, hdict: dict, overwrite=False):
    """
    Save an ami sim scene to disk (using image and hdict)

    :param params: ParamDict, parameter dictionary of constants
    :param outpath: str, the path to save the file to
    :param image: numpy array 2D, the image to save
    :param hdict: dict, the keys/values/comments to add to header

    :return: None - write to fits file
    """
    # load hdict into header
    header = fits.Header()
    # loop around keys and add to header
    if hdict is not None:
        for key in hdict:
            if len(hdict[key]) == 1:
                header[key] = hdict[key]
            else:
                header[key] = (hdict[key][0], hdict[key][1])

    primary = fits.PrimaryHDU(data=image, header=header)
    hdulist = fits.HDUList([primary])
    hdulist.writeto(outpath, overwrite=overwrite)
    hdulist.close()


# %%
FILTER = "F480M"
FOV = 81
OVERSAMPLE = 11
PSF_FILE = f"psf_{FILTER}.fits"
PA = 0.0
# PA = None
SEP = 0.4
CONTRAST = 1e-3
if PA is None:
    SKY_FILE = f"sky_{FILTER}.fits"
else:
    SKY_FILE = f"sky_{FILTER}_pa_{PA}_sep_{SEP}_contrast_{CONTRAST}.fits"
TARGET = "MYTARGET"
COUNTRATE = 1137034
# NINTS = 6214
NINTS = 200
NGROUPS = 19
APPLY_JITTER = 0
APPLY_DITHER = 0
UNIFORM_FLATFIELD = 1
RANDOM_SEED = 45
INCLUDE_PHOTNOISE = 0
INCLUDE_READNOISE = 0
INCLUDE_DARKCURRENT = 0
INCLUDE_BACKGROUND = 0
TAG = "otherseed"

# %%
psf_fits = get_psf(
    filt=FILTER, oversample=OVERSAMPLE, fov_pixels=FOV, source="A0V", filename=PSF_FILE
)

# %%
img, hdict = ami_sim_sky(FOV, OVERSAMPLE, pix_scale=U.pixscl, plot=False)
if PA is not None:
    img, hdict = sky_add_companion(img, hdict, 1, PA, SEP, CONTRAST, plot=True)
hdict["FILTER"] = (FILTER, "Input filter used")
hdict["TNAME"] = TARGET
ami_sim_save_scene(SKY_FILE, img, hdict, overwrite=True)

# %%
args = []
args += ["--target_dir", "ami_sim_output/"]
args += ["--overwrite", "1"]
args += ["--filter", FILTER]
args += ["--psf", PSF_FILE]
args += ["--oversample", OVERSAMPLE]
args += ["--countrate", COUNTRATE]
args += ["--random_seed", RANDOM_SEED]
args += ["--sky", SKY_FILE]
args += ["--nint", NINTS]
args += ["--tag", TAG]
args += ["--ngroups", NGROUPS]
args += ["--uniform_flatfield", UNIFORM_FLATFIELD]
args += ["--overwrite_flatfield", 1]
args += ["--apply_jitter", APPLY_JITTER]
args += ["--apply_dither", APPLY_DITHER]
args += ["--include_photnoise", INCLUDE_PHOTNOISE]
args += ["--include_readnoise", INCLUDE_READNOISE]
args += ["--include_darkcurrent", INCLUDE_DARKCURRENT]
args += ["--include_background", INCLUDE_BACKGROUND]
args = list(map(lambda x: str(x), args))

driver_scene.main(args)
