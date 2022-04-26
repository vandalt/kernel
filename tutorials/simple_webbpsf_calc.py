# Copied and adapted from first code cell here: https://webbpsf.readthedocs.io/en/latest/usage.html
import webbpsf
import matplotlib.pyplot as plt

ns = webbpsf.NIRISS()
ns.filter = "F480M"
ns.pupil_mask = "MASK_NRM"
psf = ns.calc_psf(oversample=4, display=True)
plt.show()
