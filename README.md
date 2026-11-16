# Prediction of leptospirosis outbreaks by hydroclimatic covariates: a comparative study of statistical models

[![DOI](https://img.shields.io/badge/DOI-10.1007%2Fs00484--022--02378--z-blue)](https://doi.org/10.1007/s00484-022-02378-z)
![Language](https://img.shields.io/badge/language-R-276DC3)

This repository contains the code to reproduce the results of:

> Llop, M. J., Gómez, A., Llop, P., López, M. S., & Müller, G. V. (2022). **Prediction of leptospirosis outbreaks by hydroclimatic covariates: a comparative study of statistical models.** *International Journal of Biometeorology*. https://doi.org/10.1007/s00484-022-02378-z

---

## 1. Description

Leptospirosis is a climate-sensitive zoonosis that is a major public health problem in northeastern Argentina. This work compares five time series forecasting methods to predict monthly leptospirosis cases using hydroclimatic covariates, with data from 2009 to 2020 for three cities (Santa Fe, Paraná and Rosario).

| Type | Method | Uses covariates |
|---|---|---|
| Parametric | ARIMA | No |
| Parametric | ARIMAX | Yes |
| Nonparametric | Nadaraya-Watson kernel estimator, one kernel (NW-1K) | No |
| Nonparametric | Nadaraya-Watson kernel estimator, two kernels (NW-2K) | No |
| Semiparametric | Semi-functional partial linear regression (SFPLR) | Yes |

**Covariates:** total monthly precipitation, monthly maximum hydrometric level of the Paraná River, and the Oceanic Niño Index (ONI).

**Main finding:** no method improves considerably on the rest and can be recommended as a unique tool for leptospirosis prediction. In general, NW-2K achieves the best performance, while parametric methods are more sensitive to outbreaks and tend to overestimate the number of cases.

The main codes are written in **R**.

## 2. Repository organization

```
.
├── CODES/
│   ├── Analysis.R    # Main script: fits the models and computes predictions and RMSE
│   ├── models.R      # Functions required by the predictors
│   ├── plots.R       # Code to generate the figures
│   └── data          # Dataset (not included, see "Data" section)
├── Bibliografia/     # Reference material
└── README.md
```

| File | Description |
|---|---|
| `Analysis.R` | Main script to reproduce the prediction results reported in the manuscript. |
| `models.R` | All the functions needed to run the predictors: `get.best.arima` (ARIMA/ARIMAX order selection), `sfplr` (SFPLR), `one.kernel` (NW-1K), `two.kernel` (NW-2K), the kernel functions and the `pesos` weights function. |
| `plots.R` | Code to generate the graphs of the manuscript (curves of the variables and predictions). |

## 3. Requirements

R and the following packages:

```r
install.packages(c("forecast", "tseries", "TSA", "readxl", "tidyverse",
                   "reshape2", "Hmisc", "fda.usc",
                   "ggplot2", "ggsci", "cowplot"))
```

| Script | Packages |
|---|---|
| `Analysis.R` | `forecast`, `tseries`, `TSA`, `readxl`, `tidyverse`, `reshape2`, `Hmisc` |
| `models.R` | `fda.usc` |
| `plots.R` | `ggplot2`, `ggsci`, `reshape2`, `cowplot` |

## 4. Data

The data are available upon request (see [Contact](#8-contact)). Once obtained, copy the dataset file into the `CODES` folder.

The dataset must contain **monthly** records from January 2009 to December 2020 (144 rows) with four columns:

| Column | Name used in the code | Description |
|---|---|---|
| 1 | `Cases` | Confirmed leptospirosis cases (response) |
| 2 | `Precip` | Total monthly precipitation (mm) |
| 3 | `Hydro_Level` | Monthly maximum hydrometric level (m) |
| 4 | `ONI` | Oceanic Niño Index |

Data sources: confirmed cases from the Ministry of Health of Santa Fe province and the Epidemiology Division of the Ministry of Health of Entre Ríos province; precipitation from the National Weather Service of Argentina (SMN) and INTA; hydrometric levels from the National Water Institute of Argentina (INA); ONI from NOAA/NWS/CPC.

## 5. Usage

The scripts assume the project folder is located at `~/CODES/`. If you place it elsewhere, update the `source("~/CODES/...")` lines in `Analysis.R` and `plots.R`.

1. **Choose the optimal lags** for each covariate (see [Figures and tables](#6-reproducing-figures-and-tables)) and set them in `Analysis.R`:

   ```r
   lp  = ...   # lag of precipitation
   lhl = ...   # lag of hydrometric level
   lo  = ...   # lag of ONI
   ```

   Use `0` for no lag.

2. **Run `Analysis.R`.** It uses 2009–2019 as the training sample and 2020 as the test sample (when covariates are lagged, the first year is dropped). The script:
   - builds the functional data (one curve per year),
   - computes the predictions of SFPLR, ARIMA, ARIMAX, NW-1K and NW-2K,
   - saves the predictions and the observed values in `results.RData`,
   - computes the RMSE of each method on the test year.

3. **Run `plots.R`** to generate the figures. The prediction plots require `results.RData` generated in the previous step.

## 6. Reproducing figures and tables

### Figures

| Figure | How to reproduce |
|---|---|
| Fig. 1 (study area map) | Not generated by the code in this repository. |
| Fig. 2 | Run `plots.R`, section `# plot of variables`. |
| Fig. 3 | Run `plots.R`, section `# plots of annual curvs of leptospirosis cases`. |
| Fig. 4 | See Convertino et al. (2021) and the [OIF Toolbox](https://github.com/HokundaiNexusLab/net-valid). |
| Figs. 5 and 6 | Run `plots.R`, section `# plots of predictions`, after running the `# Predictions` section of `Analysis.R`. |

### Tables

| Table | How to reproduce |
|---|---|
| Table 1 | See Convertino et al. (2021). |
| Tables 2 and 4 | Run `Analysis.R`, section `# RMSE`, after running the `# Predictions` section. Table 2 uses each city separately; Table 4 uses the three cities merged. |
| Table 3 | See Convertino et al. (2021). |

### External tools used in the paper

- **Transfer entropy (TE), mutual information (MI) and optimal time lags:** [OIF Toolbox](https://github.com/HokundaiNexusLab/net-valid) (Li and Convertino, 2021).
- **Variance-based sensitivity analysis:** [SAFE Toolbox](http://bristol.ac.uk/cabot/resources/safe-toolbox/) (Pianosi et al., 2015).

## 7. Citation

If you use this code, please cite:

```bibtex
@article{Llop2022leptospirosis,
  title   = {Prediction of leptospirosis outbreaks by hydroclimatic covariates: a comparative study of statistical models},
  author  = {Llop, Mar{\'i}a Jos{\'e} and G{\'o}mez, Andrea and Llop, Pamela and L{\'o}pez, Mar{\'i}a Soledad and M{\"u}ller, Gabriela V.},
  journal = {International Journal of Biometeorology},
  year    = {2022},
  doi     = {10.1007/s00484-022-02378-z}
}
```

## 8. Contact

Corresponding author: María José Llop, llopmariajose@gmail.com

Affiliations:
- Facultad de Ingeniería Química, Universidad Nacional del Litoral (UNL), Santa Fe, Argentina
- Consejo Nacional de Investigaciones Científicas y Técnicas (CONICET), Santa Fe, Argentina
- CEVARCAM, Facultad de Ingeniería y Ciencias Hídricas, UNL, Santa Fe, Argentina

The code is also available at https://www.fiq.unl.edu.ar/investigacion/investigacion-reproducible/

## References

- Convertino, M., Reddy, A., Liu, Y., Muñoz-Zanzi, C. (2021). Eco-epidemiological scaling of leptospirosis: vulnerability mapping and early warning forecasts. *Sci Total Environ* 799:149102.
- Li, J., Convertino, M. (2021). Inferring ecosystem networks as information flows. *Sci Rep* 11:7094.
- Pianosi, F., Sarrazin, F., Wagener, T. (2015). A Matlab toolbox for global sensitivity analysis. *Environ Model Softw* 70:80–85.
