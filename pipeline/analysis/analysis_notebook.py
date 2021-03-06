# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xara
from astropy.io import fits
from matplotlib import rcParams, rcParamsDefault
from astropy.nddata import Cutout2D

plt.style.use("tableau-colorblind10")
rcParams.update(rcParamsDefault)
rcParams["font.size"] = 18
rcParams["xtick.top"] = True
rcParams["xtick.direction"] = "in"
rcParams["ytick.right"] = True
rcParams["ytick.direction"] = "in"

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
rcParams["figure.autolayout"] = True
rcParams["savefig.bbox"] = "tight"

# %%
dict_wl = {
    "F277W": [2.776, 0.715],
    "F380M": [3.828, 0.205],
    "F430M": [4.286, 0.202],
    "F480M": [4.817, 0.298],
}


# %%
# phead = fits.getheader("CLEARP_CUSTOM.fits")
# pmask = fits.getdata("CLEARP_CUSTOM.fits")
# # pmask = np.flipud(fits.getdata("CLEARP_CUSTOM.fits"))
# #
# pdiam = phead["PUPLDIAM"]
# psz = pmask.shape[0]
# ppscale = pdiam / psz
# mstep = 0.3
# binary_model = True
# assert ppscale == phead["PUPLSCAL"], "Scale not the same as specified in header"
#
show = True
#
# # %%
# niriss_model = xara.core.create_discrete_model(
#     pmask, ppscale, mstep, binary=binary_model, tmin=0.4
# )
#
# # %%
# # Overplot model with data
# f0 = plt.figure(0)
# f0.clf()
# ax = f0.add_subplot(111)
# ax.imshow(pmask, origin="lower")
# ax.plot(
#     psz / 2 + niriss_model[:, 0] / ppscale, psz / 2 + niriss_model[:, 1] / ppscale, "b."
# )
# f0.set_size_inches(5, 5, forward=True)
# plt.show()

# %%
# kpo_0 = xara.KPO(array=niriss_model)
kpo_0 = xara.KPO(fname="jwst_model.txt")

# %%
kpo_0.save_as_fits("jwst_custom_binary_model.fits")

kpo_0.kpi.plot_pupil_and_uv(
    xymax=3.5, cmap=plt.cm.plasma_r, ssize=9, figsize=(10, 5), marker="o"
)
plt.show()

# %%
# data_dir = Path("ami_sim_output/F480M/")
data_dir = Path("pipeline_calibrated_data_wbadpix_corr_custom/")
# data_dir = Path("pipeline_calibrated_data_wbadpix_corr_amical/")
# data_dir = Path("pipeline_calibrated_data_wbadpix_corr/")
sep = 300
pa = 0.0
cr = 0.1
# files = [f"c_sky_F480M_pa_{pa}_sep_{sep*1e-3:.1f}_contrast_{cr}__psf_F480M__00.fits"]
# files = [f"t_sky_F480M_pa_{pa}_sep_{sep*1e-3:.1f}_contrast_{cr}__psf_F480M__00.fits"]
# files = ["t_sky_F480M__psf_F480M__00.fits"]
files = [f.name for f in data_dir.glob("*calints.fits")]

# %%
from matplotlib.colors import PowerNorm

# %%
for fname in files:
# fname = files[0]
# if True:

    print(f"Processing {fname}")
    plt.clf()
    fpath = data_dir / fname
    cube = fits.getdata(fpath)
    hdr = fits.getheader(fpath)

    # plt.imshow(cube[0], norm=PowerNorm(0.5, vmin=0.0))
    plt.imshow(cube[0])
    plt.show()

    # %%
    pscale = 65.6
    wl = dict_wl[hdr["FILTER"]][0] * 1e-6
    kpo1 = xara.KPO(fname="jwst_custom_binary_model.fits")
    # kpo1 = xara.KPO(fname="jwst_model.txt")

    # # %%
    try:
        target_name = hdr["TARGPROP"]
    except:  # noqa
        target_name = "simbin"
    # cropped_cube = []
    # for img in cube:
    #     isz = img.shape[0]
    #     x0, y0 = xara.core.determine_origin(img, algo="BCEN", verbose=False)
    #     cropped_size = int(np.min(np.floor(isz - np.array([x0, y0])))) * 2
    #     cimg = Cutout2D(img, (x0, y0), cropped_size).data
    #     cropped_cube.append(cimg)
    #
    # cube_full = cube.copy()
    #
    # cube = np.array(cropped_cube)

    # # plt.imshow(cube[0], norm=PowerNorm(0.5, vmin=0.0))
    # plt.imshow(cube[0])
    # plt.show()

    # %%
    kpo1.extract_KPD_single_cube(
        cube, pscale, wl, target=target_name, recenter=True
    )
    # kpo1.extract_KPD_single_frame(
    #     cube[0], pscale, wl, target=target_name, recenter=True
    # )

    # %%
    kp_data = np.array(kpo1.KPDT)[0]
    kp_med = np.median(kp_data, axis=0)
    # Standard error gives very small error bars. Aren't we
    # myerr = np.sqrt(np.var(kp_data, axis=0) / (kpo1.KPDT[0].shape[0] - 1))
    myerr = np.std(kp_data, axis=0)

    # %%
    plt.errorbar(np.arange(len(kp_med)), kp_med, yerr=myerr, fmt="k.")
    if show:
        plt.show()
    plt.savefig(f"plots/{fpath.stem}_kp.pdf")

    # %%
    gsize = 100
    gstep = 10  # in mas
    xx, yy = np.meshgrid(np.arange(gsize) - gsize / 2, np.arange(gsize) - gsize / 2)
    azim = -np.arctan2(xx, yy) * 180.0 / np.pi
    dist = np.hypot(xx, yy) * gstep

    # %%
    mmap = kpo1.kpd_binary_match_map(gsize, gstep, kp_med, norm=True)
    x0, y0 = np.argmax(mmap) % gsize, np.argmax(mmap) // gsize
    print(
        "max colinearity found for sep = %.2f mas and ang = %.2f deg"
        % (dist[y0, x0], azim[y0, x0])
    )

    f1 = plt.figure(figsize=(5, 5))
    ax1 = f1.add_subplot(111)
    ax1.imshow(
        mmap,
        extent=(
            gsize / 2 * gstep,
            -gsize / 2 * gstep,
            -gsize / 2 * gstep,
            gsize / 2 * gstep,
        ),
    )
    ax1.set_xlabel("right ascension (mas)")
    ax1.set_ylabel("declination (mas)")
    ax1.plot([0, 0], [0, 0], "w*", ms=16)
    ax1.set_title("Calibrated signal colinearity map")
    ax1.grid()
    f1.set_tight_layout(True)
    f1.canvas.draw()
    if show:
        plt.show()
    plt.savefig(f"plots/{fpath.stem}_colin.pdf")

    # # %%
    # p0 = [dist[y0, x0], azim[y0, x0], mmap.max()]
    #
    # mfit = kpo1.binary_model_fit(p0)
    # p1 = mfit[0]  # Best fit parameter [sep, pa, cr]
    #
    # # %%
    # cvis_b = xara.core.cvis_binary(kpo1.kpi.UVC[:, 0], kpo1.kpi.UVC[:, 1], wl, p1)
    # ker_theo = kpo1.kpi.KPM.dot(np.angle(cvis_b))
    #
    # # %%
    # fig = plt.figure(figsize=(6, 6))
    # ax = fig.add_subplot(111)
    #
    # ax.errorbar(ker_theo, kp_med, yerr=myerr, fmt="none", ecolor="c")
    # ax.plot(ker_theo, kp_med, "b.")
    # mmax = np.round(np.abs(kp_med).max())
    # ax.plot([-mmax, mmax], [-mmax, mmax], "r")
    # ax.set_ylabel("data kernel-phase")
    # ax.set_xlabel("model kernel-phase")
    # ax.set_title("kernel-phase correlation diagram")
    # # ax.axis("equal")
    # # ax.axis([-11, 11, -11, 11])
    # fig.set_tight_layout(True)
    # if show:
    #     plt.show()
    # plt.savefig(f"plots/{fpath.stem}_fit.pdf")
    #
    # if myerr is not None:
    #     chi2 = np.sum(((kp_med - ker_theo) / myerr) ** 2) / kpo1.kpi.nbkp
    # else:
    #     chi2 = np.sum(((kp_med - ker_theo)) ** 2) / kpo1.kpi.nbkp

    # break
