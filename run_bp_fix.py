from __future__ import division

# =============================================================================
# IMPORTS
# =============================================================================
import os
import argparse
import time
import glob


import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import interpolate_replace_nans
from poppy import matrixDFT
from scipy.ndimage import median_filter
from jwst.datamodels import dqflags

matplotlib.rcParams.update({'font.size': 14})




# =============================================================================
# CODE FROM ANAND FOLLOWS
# =============================================================================

micron = 1.0e-6
filts = ['F277M', 'F380M', 'F430M', 'F480M', 'F356W', 'F444W']
filtwl_d = {  # pivot wavelengths
    'F277M': 2.776e-6,  # less than Nyquist
    'F380M': 3.828e-6,
    'F430M': 4.286e-6,
    'F480M': 4.817e-6,
    'F356W': 3.595e-6,  # semi-forbidden
    'F444W': 4.435e-6,  # semi-forbidden
}
filthp_d = {  # half power limits
    'F277M': (2.413e-6, 3.142e-6),
    'F380M': (3.726e-6, 3.931e-6),
    'F430M': (4.182e-6, 4.395e-6),
    'F480M': (4.669e-6, 4.971e-6),
    'F356W': (3.141e-6, 4.068e-6),
    'F444W': (3.880e-6, 5.023e-6),
}
WL_OVERSIZEFACTOR = 0.1  # increase filter wl support by this amount to 'oversize' in wl space

pix_arcsec = 0.0656  # nominal isotropic pixel scale - refine later
pix_rad = pix_arcsec * np.pi / (60 * 60 * 180)
pupilfile_nrm = "MASK_NRM.fits"
# pupilfile_clear = "CLEARP_CUSTOM.fits"
pupilfile_clear = "CLEARP_WF.fits"
pupil_files = {
    "NRM": pupilfile_nrm,
    "CLEARP": pupilfile_clear,
}

DIAM = 6.559348  # / Flat-to-flat distance across pupil in V3 axis
PUPLDIAM = 6.603464  # / Full pupil file size, incl padding.
PUPL_CRC = 6.603464  # / Circumscribing diameter for JWST primary


def get_bad_mask(
    dqarr,
    bpix_flags="DO_NOT_USE",
) -> np.ndarray:
    """
    Get bad pixel mask from NIRISS data

    :param fpath: Path to the fits file
    :type fpath: str
    :param dq_ext: Extension with the bad pixel map, defaults to "DQ"
    :type dq_ext: str, optional
    :param bpix_flags: Bad pixel flag(s), defaults to "DO_NOT_USE"
    :type bpix_flags: Union[str, List[str]], optional
    :return: Bad pixel mask for the input file
    :rtype: np.ndarray
    :raises ValueError: ValueError is raised if bpix_flags is unknown
    """

    # Handle "all" special case where any non-zero value is flagged
    if isinstance(bpix_flags, str):
        if bpix_flags.upper() == "ALL":
            return dqarr > 0
        elif bpix_flags in dqflags.pixel:
            bpix_flags = [bpix_flags]
        else:
            raise ValueError(
                "bpix_flags should be 'all', a bad pixel flag (e.g. 'DO_NOT_USE'),"
                " or a list of bad pixel flag strings"
            )

    # Iterate flags and use element-wise OR each time to update map
    bmap = np.zeros(dqarr.shape, dtype=bool)
    for flag_name in bpix_flags:
        bmap = bmap | np.bitwise_and(dqarr, dqflags.pixel[flag_name]).astype(bool)

    return bmap


def amical_bad_pix(image, bad_map, add_bad=None, x_stddev=1):
    """
    Replace bad pixels with values interpolated from their neighbors (interpolation
    is made with a gaussian kernel convolution).

    Function from AMICAL
    """

    if add_bad is None:
        add_bad = []

    if len(add_bad) != 0:
        bad_map = bad_map.copy()  # Don't modify input bad pixel map, use a copy
        for j in range(len(add_bad)):
            bad_map[add_bad[j][1], add_bad[j][0]] = 1

    img_nan = image.copy()
    img_nan[bad_map == 1] = np.nan
    kernel = Gaussian2DKernel(x_stddev=x_stddev)
    fixed_image = interpolate_replace_nans(img_nan, kernel)
    return fixed_image


def create_wavelengths(filtername):
    """
    filtername str: filter name
    Extend filter support slightly past half power points.
    Filter transmissions are quasi-rectangular.
    """
    wl_ctr = filtwl_d[filtername]
    wl_hps = filthp_d[filtername]
    # both positive quantities below - left is lower wl, rite is higher wl
    dleft = (wl_ctr - wl_hps[0]) * (1 + WL_OVERSIZEFACTOR)
    drite = (-wl_ctr + wl_hps[1]) * (1 + WL_OVERSIZEFACTOR)

    return (wl_ctr, wl_ctr - dleft, wl_ctr + drite)


def calcsupport(filtername,
                sqfov_npix,
                pupil="NRM"):
    """
    filtername str: filter name
    calculate psf at low center high wavelengths of filter
    coadd psfs
    perform fft-style transform of image w/dft
    send back absolute value of FT(image) in filter - the CV Vsq array
    """
    wls = create_wavelengths(filtername)
    print(f"      {filtername}: {wls[0] / micron:.3f} to {wls[2] / micron:.3f} micron")

    detimage = np.zeros((sqfov_npix, sqfov_npix), float)
    for wl in wls:
        psf = calcpsf(wl, sqfov_npix, pupil=pupil)
        detimage += psf

    return transform_image(detimage)


def transform_image(image):
    ft = matrixDFT.MatrixFourierTransform()
    ftimage = ft.perform(image, image.shape[0], image.shape[0])  # fake the no-loss fft w/ dft

    return np.abs(ftimage)


def calcpsf(wl, fovnpix, pupil="NRM"):
    """
    input wl: float meters wavelength
    input fovnpix: feld of view (square) in number of pixels
    returns monochromatic unnormalized psf
    """
    reselt = wl / PUPLDIAM  # radian
    nlamD = fovnpix * pix_rad / reselt  # Soummer nlamD FOV in reselts
    # instantiate an mft object:
    ft = matrixDFT.MatrixFourierTransform()


    pupilfile = pupil_files[pupil]
    pup = pyfits.getdata(pupilfile)
    image_field = ft.perform(pup, nlamD, fovnpix)
    image_intensity = (image_field * image_field.conj()).real

    return image_intensity


# =============================================================================
# CODE FROM JENS FOLLOWS
# =============================================================================

def bad_pixels(data,
               median_size,
               median_tres):
    """
    Identify bad pixels by subtracting median-filtered data and searching for
    outliers.
    """

    mfil_data = median_filter(data, size=median_size)
    diff_data = np.abs(data - mfil_data)
    pxdq = diff_data > median_tres * np.median(diff_data)
    pxdq = pxdq.astype('int')

    print('         Identified %.0f bad pixels (%.2f%%)' % (np.sum(pxdq), np.sum(pxdq) / np.prod(pxdq.shape) * 100.))
    print('         %.3f' % np.max(diff_data/np.median(diff_data)))

    return pxdq


def fourier_corr(data,
                 pxdq,
                 fmas):
    """
    Compute and apply the bad pixel corrections based on Section 2.5 of
    Ireland 2013. This function is the core of the bad pixel cleaning code.
    """

    # Get the dimensions.
    ww = np.where(pxdq > 0.5)
    ww_ft = np.where(fmas)

    # Compute the B_Z matrix from Section 2.5 of Ireland 2013. This matrix
    # maps the bad pixels onto their Fourier power in the domain Z, which is
    # the complement of the pupil support.
    B_Z = np.zeros((len(ww[0]), len(ww_ft[0]) * 2))
    xh = data.shape[0] // 2
    yh = data.shape[1] // 2
    xx, yy = np.meshgrid(2. * np.pi * np.arange(yh + 1) / data.shape[1],
                         2. * np.pi * (((np.arange(data.shape[0]) + xh) % data.shape[0]) - xh) / data.shape[0])
    for i in range(len(ww[0])):
        cdft = np.exp(-1j * (ww[0][i] * yy + ww[1][i] * xx))
        B_Z[i, :] = np.append(cdft[ww_ft].real, cdft[ww_ft].imag)

    # Compute the corrections for the bad pixels using the Moore-Penrose pseudo
    # inverse of B_Z (Equation 19 of Ireland 2013).
    B_Z_ct = np.transpose(np.conj(B_Z))
    B_Z_mppinv = np.dot(B_Z_ct, np.linalg.inv(np.dot(B_Z, B_Z_ct)))

    # Apply the corrections for the bad pixels.
    data[ww] = 0.
    data_ft = np.fft.rfft2(data)[ww_ft]
    corr = -np.real(np.dot(np.append(data_ft.real, data_ft.imag), B_Z_mppinv))
    data[ww] += corr

    return data


def fix_bad_pixels(indir,
                   odir,
                   fitsfiles,
                   amical_correction=False,
                   show=None,
                   save=False):
    """
    """

    print('Fixing bad pixels...')


    # These values were determined empirically for NIRISS/AMI and need to be
    # tweaked for any other instrument.
    median_size = 5  # pix
    median_tres = 20. # JK: changed from 28 to 20 in order to capture all bad pixels

    # Create the output directory if it does not exist already.
    if (not os.path.exists(odir)):
        os.makedirs(odir)

    # Go through all FITS files.
    Nfitsfiles = len(fitsfiles)
    for i in range(Nfitsfiles):
        print('   File %.0f of %.0f' % (i + 1, Nfitsfiles))
        print('   %s' % fitsfiles[i])
        # Open the FITS file.
        hdul = pyfits.open(os.path.join(indir, fitsfiles[i]), memmap=False)
        data = hdul['SCI'].data
        pxdq0 = hdul['DQ'].data
        imsz = data.shape
        #print('Im size:', imsz)
        if (not (imsz[1] == 80 and imsz[2] == 80)):
            #raise UserWarning('Expecting 80x80 subarrays')
            # skip file; probably TA exposure
            print('Image dimensions in file #%i not 80x80; skipping' % (i+1))
            continue # proceed to next file
        filt = hdul[0].header['FILTER']
        pupil = hdul[0].header["PUPIL"]
        if pupil == "NRM":
            print("SKIPPING NRM FILE")
            continue
        if (filt not in filts):
            raise UserWarning('Filter ' + filt + ' is not supported')

        # code from Rachel:
        # only correct pixels marked DO_NOT_USE in the DQ array
        totpix = imsz[0] * imsz[1] * imsz[2]
        nrefpix = 4*imsz[1]*imsz[0] # 4-pixel-wide stripe on each frame
        nflagged_all = np.count_nonzero(pxdq0) - nrefpix
        DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
        dqmask = pxdq0 & DO_NOT_USE == DO_NOT_USE
        pxdq = np.where(dqmask, pxdq0, 0)
        nflagged_dnu = np.count_nonzero(pxdq) - nrefpix
        print('The following values do not include the reference pixels')
        print('\t %i pixels flagged in DQ array' % nflagged_all)
        print('\t %i pixels flagged DO_NOT_USE' % nflagged_dnu)
        print('\t %.2f percent of all pixels flagged' % ((nflagged_dnu / totpix) * 100))
        flagged_per_int = [np.count_nonzero(slc) - 4*imsz[1] for slc in dqmask]
        print('\t %f pixels (average) flagged per integration'% np.mean(flagged_per_int))

        if amical_correction:
            dq_arr = get_bad_mask(pxdq0)
            for j in range(imsz[0]):
                data[j] = amical_bad_pix(data[j], dq_arr[j])
        else:
            # These values are taken from the JDox and the SVO Filter Profile
            # Service.
            diam = PUPLDIAM  # m
            gain = 1.61  # e-/ADU
            rdns = 18.32  # e-
            pxsc = pix_arcsec * 1000.  # mas/pix

            # Find the PSF centers and determine the maximum possible frame size.
            ww_max = []
            for j in range(imsz[0]):
                ww_max += [np.unravel_index(np.argmax(median_filter(data[j], size=3)), data[j].shape)] # JK: added median filter to catch PSF center despite hot pixels
            ww_max = np.array(ww_max)
            xh = min(imsz[1] - np.max(ww_max[:, 0]), np.min(ww_max[:, 0]) - 4)  # the bottom 4 rows are reference pixels
            yh = min(imsz[2] - np.max(ww_max[:, 1]), np.min(ww_max[:, 1]) - 0)
            sh = min(xh, yh)
            print('      Cutting all frames to %.0fx%.0f pixels' % (2 * sh, 2 * sh))

            # Compute field-of-view and Fourier sampling.
            fov = 2 * sh * pxsc / 1000.  # arcsec
            fsam = filtwl_d[filt] / (fov / 3600. / 180. * np.pi)  # m/pix
            print('      FOV = %.1f arcsec, Fourier sampling = %.3f m/pix' % (fov, fsam))

            #
            cvis = calcsupport(filt,
                               2 * sh,
                               pupil=pupil,
                               )
            cvis /= np.max(cvis)
            fmas = cvis < 1e-3  # 1e-3 seems to be a reasonable threshold
            fmas_show = fmas.copy()
            fmas = np.fft.fftshift(fmas)[:, :2 * sh // 2 + 1]

            # Compute the pupil mask. This mask defines the region where we are
            # measuring the noise. It looks like 15 lambda/D distance from the PSF
            # is reasonable.
            ramp = np.arange(2 * sh) - 2 * sh // 2
            xx, yy = np.meshgrid(ramp, ramp)
            dist = np.sqrt(xx ** 2 + yy ** 2)
            pmas = dist > 8. * filtwl_d[filt] / diam * 180. / np.pi * 1000. * 3600. / pxsc # JK: changed from 15 to 8 lambda/D
            if (np.sum(pmas) <= 10.):
                raise UserWarning('Subframe too small to estimate noise')

            # Go through all frames.
            for j in range(imsz[0]):
                print('         Frame %.0f of %.0f' % (j + 1, imsz[0]))

                # Now cut out the subframe.
                data_cut = data[j, ww_max[j, 0] - sh:ww_max[j, 0] + sh, ww_max[j, 1] - sh:ww_max[j, 1] + sh].copy()
                data_orig = data_cut.copy()
                pxdq_cut = pxdq[j, ww_max[j, 0] - sh:ww_max[j, 0] + sh, ww_max[j, 1] - sh:ww_max[j, 1] + sh].copy()
                pxdq_cut = pxdq_cut > 0.
                pxdq_orig = pxdq_cut.copy()

                # Correct the bad pixels. This is an iterative process. After each
                # iteration, we check whether new (residual) bad pixels are
                # identified. If so, we re-compute the corrections. If not, we
                # terminate the iteration.
                for k in range(10):

                    # Correct the bad pixels.
                    data_cut = fourier_corr(data_cut,
                                            pxdq_cut,
                                            fmas)
                    if (k == 0):
                        data_temp = data_cut.copy()

                    # Identify residual bad pixels by looking at the high spatial
                    # frequency part of the image.
                    fmas_data = np.real(np.fft.irfft2(np.fft.rfft2(data_cut) * fmas))

                    # Analytically determine the noise (Poisson noise + read noise)
                    # and normalize the high spatial frequency part of the image
                    # by it, then identify residual bad pixels.
                    mfil_data = median_filter(data_cut, size=median_size)
                    nois = np.sqrt(mfil_data / gain + rdns ** 2)
                    fmas_data /= nois
                    temp = bad_pixels(fmas_data,
                                      median_size=median_size,
                                      median_tres=median_tres)

                    # Check which bad pixels are new. Also, compare the
                    # analytically determined noise with the empirically measured
                    # noise.
                    pxdq_new = np.sum(temp[pxdq_cut < 0.5])
                    print('         Iteration %.0f: %.0f new bad pixels, sdev of norm noise = %.3f' % (k + 1, pxdq_new, np.std(fmas_data[pmas])))

                    # If no new bad pixels were identified, terminate the
                    # iteration.
                    if (pxdq_new == 0.):
                        break

                    # If new bad pixels were identified, add them to the bad pixel
                    # map.
                    pxdq_cut = ((pxdq_cut > 0.5) | (temp > 0.5)).astype('int')

                # Put the modified subframes back into the data cube.
                data[j, ww_max[j, 0] - sh:ww_max[j, 0] + sh, ww_max[j, 1] - sh:ww_max[j, 1] + sh] = fourier_corr(data_orig,
                                                                                                                 pxdq_cut,
                                                                                                                 fmas)
                pxdq[j, ww_max[j, 0] - sh:ww_max[j, 0] + sh, ww_max[j, 1] - sh:ww_max[j, 1] + sh] = pxdq_cut

                # Show plot if desired.
                if (show is not None and str(j) in show):
                    print('Making plot...')
                    # Make plot.
                    f, ax = plt.subplots(2, 3, figsize=(3 * 6.4, 2 * 4.8))
                    p00 = ax[0, 0].imshow(np.log10(np.abs(data_orig)), origin='lower')
                    c00 = plt.colorbar(p00, ax=ax[0, 0])
                    c00.set_label('log10(abs(ADU))', labelpad=20., rotation=270.)
                    ax[0, 0].imshow(pmas, origin='lower', cmap='binary', alpha=0.25)
                    ax[0, 0].set_title('Region to estimate noise')
                    p01 = ax[0, 1].imshow(np.log10(np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(data_orig))))),
                                          origin='lower')
                    c01 = plt.colorbar(p01, ax=ax[0, 1])
                    c01.set_label('log10(abs(fourier power))', labelpad=20., rotation=270.)
                    ax[0, 1].imshow(fmas_show, origin='lower', cmap='binary', alpha=0.25)
                    ax[0, 1].set_title('Region Z (complement of pupil support)')
                    p02 = ax[0, 2].imshow(
                        np.log10(np.abs(data[j, ww_max[j, 0] - sh:ww_max[j, 0] + sh, ww_max[j, 1] - sh:ww_max[j, 1] + sh])),
                        origin='lower')
                    c02 = plt.colorbar(p02, ax=ax[0, 2])
                    c02.set_label('log10(abs(ADU))', labelpad=20., rotation=270.)
                    ax[0, 2].set_title('Data after')
                    p10 = ax[1, 0].imshow(pxdq_orig, origin='lower')
                    c10 = plt.colorbar(p10, ax=ax[1, 0])
                    ax[1, 0].set_title('Bad pixels before')
                    fmas_data = np.real(np.fft.irfft2(np.fft.rfft2(data_temp) * fmas))
                    mfil_data = median_filter(data_temp, size=median_size)
                    nois = np.sqrt(mfil_data / gain + rdns ** 2)
                    fmas_data /= nois
                    mfil_fmas_data = median_filter(fmas_data, size=median_size)
                    diff_fmas_data = np.abs(fmas_data - mfil_fmas_data)
                    temp = diff_fmas_data / np.median(diff_fmas_data)
                    p11 = ax[1, 1].imshow(temp, origin='lower', vmax=median_tres)
                    c11 = plt.colorbar(p11, ax=ax[1, 1])
                    ax[1, 1].scatter(*np.where(temp > median_tres)[::-1], s=15, c='red', marker='x')
                    ax[1, 1].set_title('Bad pixels new')
                    p12 = ax[1, 2].imshow(pxdq[j, ww_max[j, 0] - sh:ww_max[j, 0] + sh, ww_max[j, 1] - sh:ww_max[j, 1] + sh],
                                          origin='lower')
                    c12 = plt.colorbar(p12, ax=ax[1, 2])
                    ax[1, 2].set_title('Bad pixels after')
                    plt.tight_layout()
                    if save:
                        outname = os.path.join(odir, fitsfiles[i].replace('.fits','_bpfixed.png'))
                        plt.savefig(outname)
                        print('Saved figure to %s' % outname)
                    else:
                        plt.show()

        # Save the corrected data into the original FITS file.
        hdul['SCI'].data = data
        hdul['DQ'].data = pxdq0 # original dq array
        hdul.writeto(os.path.join(odir, fitsfiles[i]), output_verify='fix', overwrite=True)
        hdul.close()

    return None

def correct_fitsfiles(indir, odir, pattern='*calints.fits', show=None, save=False, use_amical=False):
    toc = time.time()
    # gather fits files to correct
    fitsfiles = sorted(map(os.path.basename, glob.glob(os.path.join(indir, pattern))))
    # do the correction
    fix_bad_pixels(indir=indir,
                   odir=odir,
                   fitsfiles=fitsfiles,
                   amical_correction=use_amical,
                   show=show,
                   save=save)
    tic = time.time()
    print("RUNTIME: %.2f s" % (tic - toc))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir",
                        help="Location of files to be fixed")
    parser.add_argument("--odir",
                        help="Location to save bad-pixel-corrected files")
    parser.add_argument("--use-amical", action="store_true", dest="use_amical",
                        help="Use amical correction")
    parser.add_argument("-p", "--pattern", default="*calints.fits", type=str,
                        help="Pattern of file names that will be processed in indir")
    parser.add_argument("--show", help="Slices to show in plot",
                        nargs='*')
    parser.add_argument("--save",  help="Save plots of corrected images",
                        action="store_true")
    args = parser.parse_args()

    correct_fitsfiles(indir=args.indir,
                      odir=args.odir,
                      pattern=args.pattern,
                      use_amical=args.use_amical,
                      show=args.show,
                      save=args.save)

