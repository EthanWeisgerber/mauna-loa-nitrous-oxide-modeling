# Forecasting Atmospheric Nitrous Oxide at Mauna Loa (1997–2023)

This project performs a comprehensive time series analysis of monthly atmospheric Nitrous Oxide (N₂O) levels recorded at Mauna Loa, Hawaii. We modeled data from 1997–2023, with a special focus on the most recent decade (2014–2023), and forecasted future values through 2026 using a variety of time series methods.

Team Members:  
> - Ethan Weisgerber  
> - Asit Das  
> - Rahul Khatri  

## Project Summary

Nitrous Oxide is a potent greenhouse gas that contributes to both global warming and ozone depletion. Using publicly available NOAA data, we developed additive time series models combining:

- **Polynomial trend modeling**
- **Trigonometric seasonality**
- **AR(17) residual modeling**

We compared multiple models (linear, exponential, polynomial) and evaluated performance using **MSE** and **MAPE**. Our final model achieved a **MAPE of 0.039%**, outperforming baseline models and Holt-Winters smoothing in interpretability.

---

## Files

| File | Description |
|------|-------------|
| `TimeSeries_Group3.pdf` | Full written report including visualizations, explanations, results, and forecasting |
| `mauna_loa_analysis.R` | R script containing all code for preprocessing, modeling, evaluation, and forecasting |

---

## Final Results

| Model | MSE | MAPE |
|-------|-----|------|
| Linear + Seasonality + White Noise | 0.091 | 0.078% |
| Exponential + Seasonality + White Noise | 0.041 | 0.050% |
| Polynomial + Seasonality + White Noise | 0.040 | 0.049% |
| Exponential + Seasonality + AR(17) | 0.026 | 0.040% |
| **Polynomial + Seasonality + AR(17)** | **0.025** | **0.039%** ✅ |

---

## Key Insights

- There is a **clear upward trend** in atmospheric N₂O levels at Mauna Loa.
- Seasonality peaks in **late autumn and early winter**, likely tied to environmental cycles.
- AR(17) modeling eliminated residual autocorrelation and improved forecast smoothness.
- Holt-Winters forecasting performed similarly well but was slightly less interpretable.

---

## Requirements

- R version ≥ 4.0
- R packages: `forecast`, `tseries`, `dplyr`

To install packages:
```r
install.packages(c("forecast", "tseries", "dplyr"))
```

## My Role

I contributed to:
- Trend modeling (linear, polynomial, exponential)
- Decomposition analysis
- Residual analysis and AR(17) modeling
- Forecast validation and performance evaluation
- Writing and visualizing sections of the final report

## Reference Data Source
NOAA Global Monitoring Laboratory:<br>
https://gml.noaa.gov/data/dataset.php?item=mlo-n2o-flask-month
