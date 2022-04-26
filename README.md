# Misc kernel phase codes

- `tutorials`: Tutorials from external sources (e.g. xara) go there
- `simulate`: Code to generate synthetic kernel phase data (e.g. with ami_sim)
- `pupil_masks`: Pupil mask fits files (the `MASK_*.fits` are the ones from
  Calwebb-kpi stage 3 made by Rachel and Jens).
- `pipeline`: Draft/test code to test things for xara or jwst pipeline
  - `badpix`: Things related to bad pixel corrections
  - `pupil`: Things related to discrete pupil model generation
  - `analysis`: Kernel phase extraction and fitting with xara or other tools

## Bad pixel correction script
For the badpixel correction script, the method used is the one from
[Ireland 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.433.1718I/abstract)
and [Kammerer et al.
2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.486..639K/abstract).

The original script is from Michael J. Ireland and was edited by Jens Kammerer
and Rachel Cooper. I also added the bad pixel correction function from
[AMICAL](https://github.com/SydneyAstrophotonicInstrumentationLab/AMICAL) to
compare its results with the fourier correction.
