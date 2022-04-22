# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xara
from astropy.io import fits
from matplotlib import rcParams, rcParamsDefault
from matplotlib.colors import PowerNorm

# %%
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

plt.style.use("tableau-colorblind10")

# %%
dict_wl = {
    "F277W": [2.776, 0.715],
    "F380M": [3.828, 0.205],
    "F430M": [4.286, 0.202],
    "F480M": [4.817, 0.298],
}

# %%
data_dir = Path("ami_sim_output/F480M/")
target_name = "binary"
calib_name = "calibrator"
sep = 0.4
cr = 0.001
tag = "otherseed"
target_file = data_dir / f"t_sky_F480M_pa_0.0_sep_{sep}_contrast_{cr}__psf_F480M_{tag}_00.fits"
calib_file = data_dir / ("c" + target_file.name[1:])
# target_file = data_dir / "t_sky_F480M_pa_0.0_sep_0.3_contrast_0.1__psf_F480M__00.fits"
# calib_file = data_dir / "c_sky_F480M_pa_0.0_sep_0.3_contrast_0.1__psf_F480M__00.fits"
# calib_file = data_dir / "t_sky_F480M__psf_F480M__00.fits"


# %%
tgt_cube = fits.getdata(target_file)[:, :-1, :-1]
hdr = fits.getheader(target_file)
cal_cube = fits.getdata(calib_file)[:, :-1, :-1]

pscale = 65.6
wl = dict_wl[hdr["FILTER"]][0] * 1e-6

# %%
# kpo1 = xara.KPO(fname="jwst_model.txt")
kpo1 = xara.KPO(fname="niriss_clear_pupil.txt")
kpo1.kpi.plot_pupil_and_uv(
    xymax=3.5, cmap=plt.cm.plasma_r, ssize=9, figsize=(10, 5), marker="o"
)
plt.show()

kpo2 = kpo1.copy()

# %%
kpo1.extract_KPD_single_cube(tgt_cube, pscale, wl, target=target_name, recenter=True)
kpo2.extract_KPD_single_cube(cal_cube, pscale, wl, target=calib_name, recenter=True)

# %%
plt.imshow(tgt_cube[0], norm=PowerNorm(0.5))
plt.show()

# %%
plt.imshow(np.mean(tgt_cube, axis=0))
plt.show()

# %%
plt.imshow(np.std(tgt_cube, axis=0))
plt.show()

# %%
plt.imshow(np.std(cal_cube, axis=0))
plt.show()

# %%
data1 = np.array(kpo1.KPDT)[0]
data2 = np.array(kpo2.KPDT)[0]

# %%
kpd1 = np.median(data1, axis=0)
# ekpd1 = np.sqrt(
#     np.var(data1, axis=0) / (kpo1.KPDT[0].shape[0] - 1)
# )
ekpd1 = np.std(data1, axis=0) / np.sqrt(data1.shape[0])
# ekpd1 = np.std(data1, axis=0)
kpd2 = np.median(data2, axis=0)
# ekpd2 = np.sqrt(
#     np.var(data2, axis=0) / (kpo2.KPDT[0].shape[0] - 1)
# )
ekpd2 = np.std(data2, axis=0) / np.sqrt(data1.shape[0])
# ekpd2 = np.std(data2, axis=0)

# %%
kpd = kpd1 - kpd2
ekpd = np.sqrt(ekpd1**2 + ekpd2**2)
# ekpd = np.sqrt(
#     np.var(data1, axis=0) / (kpo1.KPDT[0].shape[0] - 1)
#     + np.var(data2, axis=0) / (kpo2.KPDT[0].shape[0] - 1)
# )
# ekpd = np.sqrt(ekpd ** 2 + 1.2 ** 2)

# %%
x = np.arange(len(kpd1))
plt.plot(x, kpd1, label=target_name)
plt.plot(x, kpd2, label=calib_name)
np.savetxt(f"calibrator_{tag}_flipud.txt",kpd2)
print("Saved calibrator")
plt.plot(x, kpd - np.std(kpd1) * 5, label="Calibrated")
plt.legend()
plt.show()

# %%
x = np.arange(len(kpd1))
plt.errorbar(x, kpd1, yerr=ekpd1, label=target_name, fmt=".")
plt.errorbar(x, kpd2, yerr=ekpd2, label=calib_name, fmt=".")
plt.errorbar(x, kpd - np.std(kpd1) * 5, yerr=ekpd, label="Calibrated", fmt=".")
plt.legend()
plt.show()

# %%
gsize = 130
gstep = 10  # in mas
xx, yy = np.meshgrid(np.arange(gsize) - gsize / 2, np.arange(gsize) - gsize / 2)
azim = -np.arctan2(xx, yy) * 180.0 / np.pi
dist = np.hypot(xx, yy) * gstep

# %%
mmap1 = kpo1.kpd_binary_match_map(gsize, gstep, kpd1, norm=True)
x0, y0 = np.argmax(mmap1) % gsize, np.argmax(mmap1) // gsize
print(
    "max colinearity for raw target found for sep = %.2f mas and ang = %.2f deg"
    % (dist[y0, x0], azim[y0, x0])
)

# %%
mmap2 = kpo2.kpd_binary_match_map(gsize, gstep, kpd2, norm=True)
x0, y0 = np.argmax(mmap2) % gsize, np.argmax(mmap2) // gsize
print(
    "max colinearity for raw calibrator found for sep = %.2f mas and ang = %.2f deg"
    % (dist[y0, x0], azim[y0, x0])
)

# %%
mmap = kpo1.kpd_binary_match_map(gsize, gstep, kpd, norm=True)
x0, y0 = np.argmax(mmap) % gsize, np.argmax(mmap) // gsize
print(
    "max colinearity for raw calibrator found for sep = %.2f mas and ang = %.2f deg"
    % (dist[y0, x0], azim[y0, x0])
)
# %%
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

extent = (
    gsize / 2 * gstep,
    -gsize / 2 * gstep,
    -gsize / 2 * gstep,
    gsize / 2 * gstep,
)

ax1.imshow(
    mmap1,
    extent=extent,
)
ax1.set_xlabel("right ascension (mas)")
ax1.set_ylabel("declination (mas)")
ax1.plot([0, 0], [0, 0], "w*", ms=16)
ax1.set_title("Raw target colinearity map")

ax2.imshow(
    mmap2,
    extent=extent,
)
ax2.set_xlabel("right ascension (mas)")
ax2.set_ylabel("declination (mas)")
ax2.plot([0, 0], [0, 0], "w*", ms=16)
ax2.set_title("Raw calibrator colinearity map")

ax3.imshow(
    mmap,
    extent=extent,
)
ax3.set_xlabel("right ascension (mas)")
ax3.set_ylabel("declination (mas)")
ax3.plot([0, 0], [0, 0], "w*", ms=16)
ax3.set_title("Calibrated colinearity map")

fig.set_tight_layout(True)
plt.show()

# %%
p0 = [dist[y0, x0], azim[y0, x0], mmap.max()]

mfit = kpo1.binary_model_fit(p0, calib=kpo2)
p1 = mfit[0]  # Best fit parameter [sep, pa, cr]

cvis_b = xara.core.cvis_binary(kpo1.kpi.UVC[:, 0], kpo1.kpi.UVC[:, 1], wl, p1)
ker_theo = kpo1.kpi.KPM.dot(np.angle(cvis_b))

# %%
plt.errorbar(ker_theo, kpd, yerr=ekpd, fmt="none", ecolor="c")
plt.plot(ker_theo, kpd, "b.")
mmax = np.round(np.abs(kpd).max())
plt.plot([-mmax, mmax], [-mmax, mmax], "r")
plt.ylabel("data kernel-phase")
plt.xlabel("model kernel-phase")
plt.title("kernel-phase correlation diagram")
plt.axis("equal")
# plt.axis([-11, 11, -11, 11])
fig.set_tight_layout(True)
plt.show()

# %%
plt.errorbar(x, kpd, yerr=np.sqrt(ekpd1**2 + ekpd2**2), label="Calibrated", fmt=".")
plt.plot(x, ker_theo, label="Model")
plt.show()

# %%
chi2 = np.sum(((kpd - ker_theo) / ekpd) ** 2) / kpo1.kpi.nbkp

# %%
print(p1)
print(chi2)
