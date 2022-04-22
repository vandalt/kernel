from pathlib import Path

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm, LogNorm

parent_data_dir = Path("data")
wbad_dir  = parent_data_dir / Path("pipeline_calibrated_data_wbadpix/")
old_corr_dir = parent_data_dir / Path("pipeline_calibrated_data_wbadpix_corr/")
custom_corr_dir = parent_data_dir / Path("pipeline_calibrated_data_wbadpix_corr_custom/")
wf_corr_dir = parent_data_dir / Path("pipeline_calibrated_data_wbadpix_corr_wf/")
amical_corr_dir = parent_data_dir / Path("pipeline_calibrated_data_wbadpix_corr_amical/")

corrections = {
    "No correction": wbad_dir,
    # "NRM correction (bad)": old_corr_dir,
    "Custom mask": custom_corr_dir,
    # "WF mask":wf_corr_dir,
    "AMICAL correction (Gaussian convolution)": amical_corr_dir,
}

files = [f.name for f in custom_corr_dir.glob("*calints.fits")]

# %%
# ref_key = "Custom mask"
ref_key = "No correction"
ref_corr = corrections[ref_key]
ref2_key = "Custom mask"
ref2_corr = corrections[ref2_key]
for fname in files:

    fig, axs = plt.subplots(figsize=(19 * len(corrections) / 5, 1.5 * 8), nrows=3, ncols=len(corrections), sharex="col")

    fref = ref_corr / fname
    ref = fits.getdata(fref)
    fref2 = ref2_corr / fname
    ref2 = fits.getdata(fref2)
    for i, corr in enumerate(corrections):

        fpath = corrections[corr] / fname

        data = fits.getdata(fpath)

        psf = axs[0, i].imshow(data[0], norm=PowerNorm(0.5, vmin=0), cmap="afmhot")
        # psf = axs[0, i].imshow(data[0], norm=LogNorm(), cmap="afmhot")
        axs[0, i].set_title(corr)

        res = axs[1, i].imshow(data[0] - ref[0], norm=PowerNorm(0.1, vmin=0.0), cmap="afmhot")
        # res = axs[1, i].imshow(data[0] - ref[0], cmap="afmhot")
        # res = axs[1, i].imshow(data[0] - ref[0], norm=LogNorm(), cmap="afmhot")
        if i == 0:
            axs[1, i].set_ylabel(f"{corr} - {ref_key}")

        res2 = axs[2, i].imshow(data[0] - ref2[0], norm=PowerNorm(0.1, vmin=0.0), cmap="afmhot")
        # res = axs[1, i].imshow(data[0] - ref2[0], cmap="afmhot")
        # res = axs[1, i].imshow(data[0] - ref2[0], norm=LogNorm(), cmap="afmhot")
        if i == 0:
            axs[2, i].set_ylabel(f"{corr} - {ref2_key}")

        fig.colorbar(psf, ax=axs[0, i], shrink=0.6)
        fig.colorbar(res, ax=axs[1, i], shrink=0.6)
        fig.colorbar(res2, ax=axs[2, i], shrink=0.6)
    plt.tight_layout()
    plt.show()
    break
