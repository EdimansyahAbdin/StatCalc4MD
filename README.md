# StatCalc MD

**Sample size calculator for medical researchers and clinicians**

🌐 **Live site:** [edimansyahabdin.github.io/statcalc](https://edimansyahabdin.github.io/statcalc)

---

## About

StatCalc MD is a free, open-source sample size calculator built for medical researchers, clinicians, and biostatistics students. In Mobile App version, every calculation runs locally in your browser — no data is sent anywhere and the app works fully offline after the first load.

---

## Features

- **26 study designs** across 10 categories
- **Ready-to-paste manuscript paragraph** for each design — with embedded textbook citations and exact page numbers
- **Verified R, Stata, and SAS code** to independently reproduce every calculation
- **12 authoritative reference books** including Cohen (1988), Machin et al. (2009), Chow et al. (2008), Stata PSS Manual (2023), SAS PROC POWER, WebPower, semPower, and more
- **Share via WhatsApp** or **export to PDF**
- **Works fully offline in Mobile App version** — no account, no sign-in, no data collection
- **iPhone-optimised** mobile layout

---

## Study designs covered

| Category | Designs |
|----------|---------|
| **Means** | Two independent means, Paired means, One-sample mean |
| **Proportions** | Two independent proportions, One-sample proportion, McNemar (paired), Chi-square, Non-inferiority |
| **ANOVA** | One-way, Two-way, Repeated measures, ANCOVA |
| **Regression** | Pearson correlation, Two correlations, Simple slope, Multiple R², Logistic |
| **Survival** | Log-rank, Cluster RCT log-rank, Cox proportional hazards |
| **Clinical Trials** | Equivalence (TOST) |
| **Diagnostic** | ICC reliability / repeatability |
| **Phase II** | Simon two-stage (oncology) |
| **Survey** | Prevalence / proportion, Mean estimation |
| **SEM / Mediation** | RMSEA-based (semPower), Sobel mediation test |

---

## Reference sources

1. Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.)
2. Aberson, C.L. (2019). *Applied Power Analysis for the Behavioral Sciences* (2nd ed.)
3. Zhang, Z. & Yuan, K.-H. (2018). *Practical Statistical Power Analysis Using WebPower and R*
4. Jak, S. et al. (2020). Analytical power calculations for SEM. *Behavior Research Methods*
5. Machin, D., Campbell, M.J., & Tan, S.B. (2009). *Sample Size Tables for Clinical Studies* (3rd ed.)
6. Chow, S.-C., Shao, J., & Wang, H. (2008). *Sample Size Calculations in Clinical Research* (2nd ed.)
7. Aday, L.A. & Cornelius, L.J. (2006). *Designing and Conducting Health Surveys* (3rd ed.)
8. Valliant, R., Dever, J.A., & Kreuter, F. (2013). *Practical Tools for Designing and Weighting Survey Samples*
9. Fleiss, J.L., Levin, B., & Paik, M.C. (2003). *Statistical Methods for Rates and Proportions* (3rd ed.)
10. Moshagen, M. & Bader, M. (2024). semPower. *Behavior Research Methods* 56:2901–2922
11. StataCorp (2023). *Stata 18 Power, Precision, and Sample-Size Reference Manual*
12. SAS Institute (2015). *SAS/STAT 14.1 User's Guide, Chapter 89: The POWER Procedure*

---

## Install on iPhone (PWA)

1. Open [edimansyahabdin.github.io/statcalc-md](https://edimansyahabdin.github.io/statcalc) in **Safari**
2. Tap the **Share** button (square with arrow)
3. Tap **"Add to Home Screen"**
4. Tap **Add** — the app icon appears on your home screen
5. Opens full-screen like a native app, works offline ✅

---

## Run locally

```bash
git clone https://github.com/EdimansyahAbdin/statcalc-md.git
cd statcalc-md
npm install
npm run dev
# Open http://localhost:3000
```

---

## Tech stack

- **React 18** + **Vite 5**
- Pure JavaScript formula engine (no external stats library)
- Deployed via **GitHub Actions → GitHub Pages**
- PWA with offline service worker

---

## Disclaimer

StatCalc MD is intended for research planning purposes only and does not constitute medical or statistical advice. All sample size calculations should be reviewed by a qualified statistician before use in clinical research protocols or ethics submissions.

---

## License

MIT — free to use, share, and adapt with attribution.
