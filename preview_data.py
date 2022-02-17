from pathlib import Path

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

wbad_dir  = Path("pipeline_calibrated_data_wbadpix/")
old_corr_dir = Path("pipeline_calibrated_data_wbadpix_corr/")
custom_corr_dir = Path("pipeline_calibrated_data_wbadpix_corr_custom/")
wf_corr_dir = Path("pipeline_calibrated_data_wbadpix_corr_wf/")
amical_corr_dir = Path("pipeline_calibrated_data_wbadpix_corr_amical/")

corrections = {
    "No correction": wbad_dir,
    "NRM correction (bad)": old_corr_dir,
    "Custom mask": custom_corr_dir,
    "WF mask":wf_corr_dir,
    "AMICAL correction (Gaussian convolution)": amical_corr_dir,
}

files = [
    "jw01093008001_01101_00006_nis_calints.fits",
    "jw01093008001_01101_00007_nis_calints.fits",
    "jw01093009001_01101_00006_nis_calints.fits",
    "jw01093009001_01101_00007_nis_calints.fits",
    "jw01093010001_01101_00006_nis_calints.fits",
    "jw01093010001_01101_00007_nis_calints.fits",
    "jw01093011001_01101_00006_nis_calints.fits",
    "jw01093011001_01101_00007_nis_calints.fits",
    "jw01093008001_01101_00006_nis_calints.fits",
    "jw01093008001_01101_00007_nis_calints.fits",
    "jw01093009001_01101_00006_nis_calints.fits",
    "jw01093009001_01101_00007_nis_calints.fits",
    "jw01093010001_01101_00006_nis_calints.fits",
    "jw01093010001_01101_00007_nis_calints.fits",
    "jw01093011001_01101_00006_nis_calints.fits",
    "jw01093011001_01101_00007_nis_calints.fits",
]


# %%
ref_corr = corrections["Custom mask"]
for fname in files:

    fig, axs = plt.subplots(figsize=(19, 8), nrows=2, ncols=len(corrections), sharex="col")

    fref = ref_corr / fname
    ref = fits.getdata(fref)
    for i, corr in enumerate(corrections):

        fpath = corrections[corr] / fname

        data = fits.getdata(fpath)

        psf = axs[0, i].imshow(data[0], norm=PowerNorm(0.5, vmin=0), cmap="afmhot")
        axs[0, i].set_title(corr)

        res = axs[1, i].imshow(data[0] - ref[0], norm=PowerNorm(0.5, vmin=0), cmap="afmhot")
        if i == 0:
            axs[1, i].set_ylabel("Residuals vs. 'Custom mask'")

        fig.colorbar(psf, ax=axs[0, i], shrink=0.6)
        fig.colorbar(res, ax=axs[1, i], shrink=0.6)
    plt.tight_layout()
    plt.show()
    break
