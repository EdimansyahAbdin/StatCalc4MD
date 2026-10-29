import { useState, useCallback } from "react";

// ─────────────────────────────────────────────
// FORMULA ENGINE — pure deterministic math
// ─────────────────────────────────────────────
const jStat = {
  // inverse normal CDF (Beasley-Springer-Moro)
  normalInv(p) {
    const a = [
      -3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2,
       1.383577518672690e2, -3.066479806614716e1, 2.506628277459239,
    ];
    const b = [
      -5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2,
       6.680131188771972e1, -1.328068155288572e1,
    ];
    const c = [
      -7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838,
      -2.549732539343734, 4.374664141464968, 2.938163982698783,
    ];
    const d = [
      7.784695709041462e-3, 3.224671290700398e-1,
      2.445134137142996, 3.754408661907416,
    ];
    const pLow = 0.02425, pHigh = 1 - pLow;
    let q, r;
    if (p < pLow) {
      q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    } else if (p <= pHigh) {
      q = p - 0.5; r = q * q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
             (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    } else {
      q = Math.sqrt(-2 * Math.log(1 - p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
              ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  },
  // normal CDF
  normalCDF(x) {
    const t = 1 / (1 + 0.2316419 * Math.abs(x));
    const d = 0.3989422820 * Math.exp(-x * x / 2);
    let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.7814779 + t * (-1.8212560 + t * 1.3302744))));
    return x > 0 ? 1 - p : p;
  },
  // chi-square inverse CDF (Wilson-Hilferty approximation)
  chiInv(p, df) {
    const z = this.normalInv(p);
    const h = 2 / (9 * df);
    return df * Math.pow(1 - h + z * Math.sqrt(h), 3);
  },
  // non-central chi-square power
  ncchiPower(lambda, df, alpha) {
    const critVal = this.chiInv(1 - alpha, df);
    // approximation: shift mean by lambda
    const z = (critVal - df - lambda) / Math.sqrt(2 * (df + 2 * lambda));
    return 1 - this.normalCDF(z);
  },
  // regularized lower incomplete gamma P(Gamma(a)≤x) via series + Lentz CF
  gammaInc(a, x) {
    if (x <= 0) return 0;
    if (x < a + 1) {
      let del = 1, sum = 1;
      for (let n = 1; n <= 400; n++) { del *= x / (a + n); sum += del; if (del < 1e-14) break; }
      return Math.exp(-x + a * Math.log(x) - this.lgamma(a + 1)) * sum;
    }
    const FPMIN = 1e-30;
    let b = x + 1 - a, c = 1 / FPMIN, d = 1 / b, h = d;
    for (let i = 1; i <= 200; i++) {
      const an = -i * (i - a); b += 2;
      d = an * d + b; if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d;
      c = b + an / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-12) break;
    }
    return 1 - Math.exp(-x + a * Math.log(x) - this.lgamma(a)) * h;
  },
  // exact chi-squared CDF
  chiCDF(x, df) { return this.gammaInc(df / 2, x / 2); },
  // exact chi-squared quantile via bisection
  chiInvExact(p, df) {
    if (df === 1) { const z = this.normalInv((1 + p) / 2); return z * z; }
    let lo = 0, hi = Math.max(20, 2 * df);
    while (this.chiCDF(hi, df) < p) hi *= 2;
    for (let i = 0; i < 100; i++) {
      const mid = (lo + hi) / 2;
      if (this.chiCDF(mid, df) < p) lo = mid; else hi = mid;
      if (hi - lo < 1e-9) break;
    }
    return (lo + hi) / 2;
  },
  // exact non-central chi-squared power (Poisson mixture of gammaInc)
  ncchiPowerExact(lambda, df, alpha) {
    const chiCrit = this.chiInvExact(1 - alpha, df);
    if (lambda <= 0) return 1 - this.chiCDF(chiCrit, df);
    const halfNcp = lambda / 2;
    let logW = -halfNcp, power = 0;
    for (let j = 0; j <= 300; j++) {
      const w = Math.exp(logW);
      if (w < 1e-15 && j > halfNcp * 3) break;
      power += w * (1 - this.gammaInc((df + 2 * j) / 2, chiCrit / 2));
      logW += Math.log(halfNcp) - Math.log(j + 1);
    }
    return power;
  },
  // non-central chi-squared CDF (Poisson mixture of gammaInc)
  ncchiCDF(x, df, lambda) {
    if (x <= 0) return 0;
    if (lambda <= 0) return this.gammaInc(df / 2, x / 2);
    const halfNcp = lambda / 2;
    let logW = -halfNcp, cdf = 0;
    for (let j = 0; j <= 300; j++) {
      const w = Math.exp(logW);
      if (w < 1e-15 && j > halfNcp * 3) break;
      cdf += w * this.gammaInc((df + 2 * j) / 2, x / 2);
      logW += Math.log(halfNcp) - Math.log(j + 1);
    }
    return cdf;
  },
  // non-central chi-squared quantile via bisection
  ncchiInv(p, df, lambda) {
    if (lambda <= 0) return this.chiInvExact(p, df);
    let lo = 0, hi = df + lambda + 10 * Math.sqrt(df + 2 * lambda);
    while (this.ncchiCDF(hi, df, lambda) < p) hi *= 2;
    for (let i = 0; i < 100; i++) {
      const mid = (lo + hi) / 2;
      if (this.ncchiCDF(mid, df, lambda) < p) lo = mid; else hi = mid;
      if (hi - lo < 1e-9) break;
    }
    return (lo + hi) / 2;
  },
  // inverse t-CDF (Cornish-Fisher series approximation)
  tInv(p, df) {
    const z = this.normalInv(p);
    const g1 = (z * z * z + z) / 4;
    const g2 = (5 * z * z * z * z * z + 16 * z * z * z + 3 * z) / 96;
    return z + g1 / df + g2 / (df * df);
  },
  // non-central t power: corrected normal approximation matching R's power.t.test
  nctPower(ncp, t_crit, df) {
    const delta = ncp - t_crit;
    const phi = Math.exp(-delta * delta / 2) / Math.sqrt(2 * Math.PI);
    return this.normalCDF(delta) - 0.05 * (t_crit * ncp / (2 * df)) * phi;
  },
  // log-gamma (Lanczos approximation)
  lgamma(x) {
    const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
               -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - this.lgamma(1 - x);
    x--; let a = c[0]; const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += c[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  },
  // regularized incomplete beta (Lentz continued fraction)
  betaInc(x, a, b) {
    if (x <= 0) return 0; if (x >= 1) return 1;
    if (x > (a + 1) / (a + b + 2)) return 1 - this.betaInc(1 - x, b, a);
    const front = Math.exp(a * Math.log(x) + b * Math.log(1 - x) - this.lgamma(a) - this.lgamma(b) + this.lgamma(a + b)) / a;
    let f = 1, C = 1, D = 1 - (a + b) * x / (a + 1);
    if (Math.abs(D) < 1e-30) D = 1e-30; D = 1 / D; f = D;
    for (let m = 1; m <= 200; m++) {
      let num = m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m));
      D = 1 + num*D; if (Math.abs(D) < 1e-30) D = 1e-30; D = 1/D;
      C = 1 + num/C; if (Math.abs(C) < 1e-30) C = 1e-30; f *= C*D;
      num = -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1));
      D = 1 + num*D; if (Math.abs(D) < 1e-30) D = 1e-30; D = 1/D;
      C = 1 + num/C; if (Math.abs(C) < 1e-30) C = 1e-30; f *= C*D;
      if (Math.abs(C*D - 1) < 1e-12) break;
    }
    return front * f;
  },
  // F distribution CDF
  fCDF(x, d1, d2) { if (x <= 0) return 0; return this.betaInc(d1*x / (d1*x + d2), d1/2, d2/2); },
  // F distribution quantile (bisection)
  fInv(p, d1, d2) {
    let lo = 0, hi = 1000;
    while (this.fCDF(hi, d1, d2) < p) hi *= 2;
    for (let i = 0; i < 100; i++) {
      const mid = (lo + hi) / 2;
      if (this.fCDF(mid, d1, d2) < p) lo = mid; else hi = mid;
      if (hi - lo < 1e-10) break;
    }
    return (lo + hi) / 2;
  },
  // Non-central F power (exact Poisson mixture of beta distributions)
  ncfPower(ncp, df1, df2, alpha) {
    const Fcrit = this.fInv(1 - alpha, df1, df2);
    const x = df1 * Fcrit / (df1 * Fcrit + df2);
    const halfNcp = ncp / 2;
    let logW = -halfNcp, power = 0;
    for (let j = 0; j <= 300; j++) {
      const w = Math.exp(logW);
      if (w < 1e-15 && j > 5) break;
      power += w * (1 - this.betaInc(x, df1/2 + j, df2/2));
      logW += Math.log(halfNcp) - Math.log(j + 1);
    }
    return power;
  },
};

const za = (alpha, sides = 2) => jStat.normalInv(1 - alpha / sides);
const zb = (power) => jStat.normalInv(power);

const FORMULAS = {
  two_means: {
    label: "Two independent means",
    group: "Means",
    params: ["effectMode", "mu1", "mu2", "sd", "d", "alpha", "power", "sides", "attrition", "ratio"],
    defaults: { effectMode: "means", mu1: 10, mu2: 8, sd: 5, d: 0.5, alpha: 0.05, power: 0.80, sides: 2, attrition: 0, ratio: 1 },
    labels: { effectMode: "Effect size input", mu1: "Mean (group 1)", mu2: "Mean (group 2)", sd: "Common SD", d: "Cohen's d", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %", ratio: "Allocation ratio (n2/n1)" },
    units: { mu1: "", mu2: "", sd: "", d: "", alpha: "", power: "", sides: "", attrition: "%", ratio: "" },
    hints: { mu1: "Expected mean in group 1", mu2: "Expected mean in group 2 (the difference that matters clinically)", sd: "Pooled standard deviation (from pilot data or literature)", d: "Standardised difference (|μ₁−μ₂|/σ): small=0.20, medium=0.50, large=0.80", attrition: "Expected dropout rate — n will be inflated to compensate" },
    visibleParams(p) {
      const tail = ["alpha", "power", "sides", "attrition", "ratio"];
      if (p.effectMode === "d") return ["effectMode", "d", ...tail];
      return ["effectMode", "mu1", "mu2", "sd", ...tail];
    },
    calculate(p) {
      const d = p.effectMode === "d" ? Math.abs(p.d) : Math.abs(p.mu1 - p.mu2) / p.sd;
      if (d <= 0) throw new Error("Effect size must be > 0");
      let n1 = Math.max(2, Math.ceil(Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) * (1 + 1 / p.ratio) / (d * d)));
      for (let iter = 0; iter < 300; iter++) {
        const n2t = Math.max(2, Math.ceil(n1 * p.ratio));
        const df = n1 + n2t - 2;
        const ncp = d * Math.sqrt(n1 * n2t / (n1 + n2t));
        const t_crit = jStat.tInv(1 - p.alpha / p.sides, df);
        if (jStat.nctPower(ncp, t_crit, df) >= p.power) break;
        n1++;
      }
      const n1_adj = Math.ceil(n1 / (1 - p.attrition / 100));
      const n2 = Math.ceil(n1_adj * p.ratio);
      return { n1: n1_adj, n2, total: n1_adj + n2, effectSize: d, effectLabel: "Cohen's d" };
    },
  },
  paired_means: {
    label: "Paired / matched means",
    group: "Means",
    params: ["mu_diff", "sd_diff", "alpha", "power", "sides", "attrition"],
    defaults: { mu_diff: 5, sd_diff: 10, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { mu_diff: "Expected mean difference", sd_diff: "SD of differences", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { mu_diff: "Mean of (post − pre) or (treatment − control)", sd_diff: "SD of paired differences — often √2·(1−ρ)·σ" },
    calculate(p) {
      const d = Math.abs(p.mu_diff) / p.sd_diff;
      // paired t-test: df = n-1, ncp = d*sqrt(n)
      let n = Math.max(2, Math.ceil(Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) / (d * d)));
      for (let iter = 0; iter < 300; iter++) {
        const df = n - 1;
        const ncp = d * Math.sqrt(n);
        const t_crit = jStat.tInv(1 - p.alpha / p.sides, df);
        if (jStat.nctPower(ncp, t_crit, df) >= p.power) break;
        n++;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: n_adj, total: n_adj, effectSize: d, effectLabel: "Cohen's d (paired)" };
    },
  },
  one_mean: {
    label: "One-sample mean",
    group: "Means",
    params: ["mu0", "mu1", "sd", "alpha", "power", "sides", "attrition"],
    defaults: { mu0: 0, mu1: 3, sd: 8, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { mu0: "Null value (μ₀)", mu1: "Expected mean (μ₁)", sd: "SD", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { mu0: "The hypothesised population value under H₀", mu1: "The true mean you expect to observe" },
    calculate(p) {
      const d = Math.abs(p.mu1 - p.mu0) / p.sd;
      // one-sample t-test: df = n-1, ncp = d*sqrt(n)
      let n = Math.max(2, Math.ceil(Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) / (d * d)));
      for (let iter = 0; iter < 300; iter++) {
        const df = n - 1;
        const ncp = d * Math.sqrt(n);
        const t_crit = jStat.tInv(1 - p.alpha / p.sides, df);
        if (jStat.nctPower(ncp, t_crit, df) >= p.power) break;
        n++;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj, effectSize: d, effectLabel: "Cohen's d" };
    },
  },
  two_proportions: {
    label: "Two independent proportions",
    group: "Proportions",
    params: ["p1", "p2", "alpha", "power", "sides", "attrition", "ratio"],
    defaults: { p1: 0.30, p2: 0.50, alpha: 0.05, power: 0.80, sides: 2, attrition: 0, ratio: 1 },
    labels: { p1: "Proportion (group 1)", p2: "Proportion (group 2)", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %", ratio: "Allocation ratio (n2/n1)" },
    hints: { p1: "Event rate in control group (from literature)", p2: "Expected event rate in intervention group" },
    calculate(p) {
      const pbar = (p.p1 + p.ratio * p.p2) / (1 + p.ratio);
      const num = Math.pow(
        za(p.alpha, p.sides) * Math.sqrt((1 + 1 / p.ratio) * pbar * (1 - pbar)) +
        zb(p.power) * Math.sqrt(p.p1 * (1 - p.p1) + p.p2 * (1 - p.p2) / p.ratio),
        2
      );
      const denom = Math.pow(p.p1 - p.p2, 2);
      const n1 = Math.ceil(num / denom / (1 - p.attrition / 100));
      const n2 = Math.ceil(n1 * p.ratio);
      const h = 2 * (Math.asin(Math.sqrt(p.p2)) - Math.asin(Math.sqrt(p.p1)));
      return { n1, n2, total: n1 + n2, effectSize: Math.abs(h), effectLabel: "Cohen's h" };
    },
  },
  one_proportion: {
    label: "One-sample proportion",
    group: "Proportions",
    params: ["p0", "p1", "alpha", "power", "sides", "attrition"],
    defaults: { p0: 0.50, p1: 0.65, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { p0: "Null proportion (p₀)", p1: "Expected proportion (p₁)", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { p0: "The reference or historical proportion under H₀", p1: "The true proportion you expect to detect" },
    calculate(p) {
      const h = Math.abs(2 * (Math.asin(Math.sqrt(p.p1)) - Math.asin(Math.sqrt(p.p0))));
      const n = Math.ceil(Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) / (h * h));
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj, effectSize: h, effectLabel: "Cohen's h" };
    },
  },
  one_way_anova: {
    label: "One-way ANOVA",
    group: "ANOVA",
    params: ["f", "k", "alpha", "power", "attrition"],
    defaults: { f: 0.25, k: 3, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { f: "Effect size (Cohen's f)", k: "Number of groups", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { f: "Small=0.10, Medium=0.25, Large=0.40 (Cohen 1988, Ch.8)", k: "Total number of independent groups" },
    calculate(p) {
      const df1 = p.k - 1;
      let n_per = 2;
      for (let iter = 0; iter < 1000; iter++) {
        const N = n_per * p.k;
        const df2 = N - p.k;
        const lambda = p.f * p.f * N;
        if (jStat.ncfPower(lambda, df1, df2, p.alpha) >= p.power) break;
        n_per++;
      }
      const n_adj = Math.ceil(n_per / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj * p.k, perGroup: n_adj, effectSize: p.f, effectLabel: "Cohen's f" };
    },
  },
  correlation: {
    label: "Pearson correlation",
    group: "Regression",
    params: ["rho", "alpha", "power", "sides", "attrition"],
    defaults: { rho: 0.30, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { rho: "Expected correlation (ρ)", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { rho: "Small=0.10, Medium=0.30, Large=0.50 (Cohen 1988, Ch.3)" },
    calculate(p) {
      const z_rho = 0.5 * Math.log((1 + p.rho) / (1 - p.rho));
      const n = Math.ceil(Math.pow((za(p.alpha, p.sides) + zb(p.power)) / z_rho, 2) + 2.2);
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj, effectSize: p.rho, effectLabel: "ρ (Pearson)" };
    },
  },
  regression_r2: {
    label: "Multiple regression (R²)",
    group: "Regression",
    params: ["f2", "u", "alpha", "power", "attrition"],
    defaults: { f2: 0.15, u: 3, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { f2: "Effect size (Cohen's f²)", u: "Number of predictors", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { f2: "Small=0.02, Medium=0.15, Large=0.35 (Cohen 1988, Ch.9)", u: "Number of predictors in the model" },
    calculate(p) {
      const zApprox = jStat.normalInv(1 - p.alpha / 2) + jStat.normalInv(p.power);
      let n = Math.max(p.u + 5, Math.floor(zApprox * zApprox / p.f2));
      while (true) {
        const df2 = n - p.u - 1;
        if (df2 > 0 && jStat.ncfPower(p.f2 * n, p.u, df2, p.alpha) >= p.power) break;
        n++;
        if (n > 100000) break;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj, effectSize: p.f2, effectLabel: "Cohen's f²" };
    },
  },
  logrank: {
    label: "Log-rank test (survival)",
    group: "Survival",
    params: ["hr", "p_event", "alpha", "power", "sides", "attrition", "ratio"],
    defaults: { hr: 0.70, p_event: 0.60, alpha: 0.05, power: 0.80, sides: 2, attrition: 0, ratio: 1 },
    labels: { hr: "Hazard ratio (HR)", p_event: "Overall event probability", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %", ratio: "Allocation ratio (n2/n1)" },
    hints: { hr: "HR < 1 = intervention better; minimum clinically important HR to detect", p_event: "Proportion of participants expected to experience the event" },
    calculate(p) {
      // Schoenfeld (1981) formula
      const events = Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) / Math.pow(Math.log(p.hr), 2) * Math.pow(1 + p.ratio, 2) / p.ratio;
      const e = Math.ceil(events);
      const n_total_raw = e / p.p_event;
      const n1 = Math.ceil(n_total_raw / (1 + p.ratio) / (1 - p.attrition / 100));
      const n2 = Math.ceil(n1 * p.ratio);
      return { n1, n2, total: n1 + n2, events: e, effectSize: p.hr, effectLabel: "Hazard ratio" };
    },
  },
  survey_proportion: {
    label: "Survey — prevalence / proportion",
    group: "Survey",
    params: ["p", "d", "alpha", "deff", "response_rate"],
    defaults: { p: 0.50, d: 0.05, alpha: 0.05, deff: 1.0, response_rate: 100 },
    labels: { p: "Expected proportion", d: "Margin of error (±d)", alpha: "Alpha (α)", deff: "Design effect (DEFF)", response_rate: "Expected response rate %" },
    hints: { p: "Use 0.50 for maximum sample size when unknown", d: "Desired precision — e.g. 0.05 = ±5%", deff: "1.0 for simple random sample; >1 for cluster sampling", response_rate: "Inflate n to account for non-response" },
    calculate(p) {
      const z = za(p.alpha, 2);
      const n_raw = (z * z * p.p * (1 - p.p)) / (p.d * p.d);
      const n_deff = Math.ceil(n_raw * p.deff);
      const n_final = Math.ceil(n_deff / (p.response_rate / 100));
      return { n1: n_final, n2: null, total: n_final, effectSize: null, effectLabel: null };
    },
  },
  sem_rmsea: {
    label: "SEM — RMSEA-based power",
    group: "SEM",
    params: ["rmsea0", "rmsea1", "df", "alpha", "power"],
    defaults: { rmsea0: 0.05, rmsea1: 0.08, df: 10, alpha: 0.05, power: 0.80 },
    labels: { rmsea0: "RMSEA under H₀", rmsea1: "RMSEA under H₁", df: "Model degrees of freedom", alpha: "Alpha (α)", power: "Power (1-β)" },
    hints: { rmsea0: "Typically 0.05 (close fit threshold)", rmsea1: "Expected true RMSEA — 0.08 = mediocre fit", df: "df = observed statistics − free parameters" },
    calculate(p) {
      let n = 50;
      for (let iter = 0; iter < 2000; iter++) {
        const lambda0 = p.rmsea0 * p.rmsea0 * p.df * n;
        const lambda1 = p.rmsea1 * p.rmsea1 * p.df * n;
        const crit = jStat.ncchiInv(1 - p.alpha, p.df, lambda0);
        const pow = 1 - jStat.ncchiCDF(crit, p.df, lambda1);
        if (pow >= p.power) break;
        n++;
      }
      return { n1: n, n2: null, total: n, effectSize: p.rmsea1, effectLabel: "RMSEA" };
    },
  },

  // ── NEW DESIGNS ─────────────────────────────────────────────────────────────

  // Repeated measures ANOVA (O'Brien & Kaiser 1985; Aberson 2019 Ch.6; Zhang Ch.8; Stata PSS-2 p.393)
  repeated_measures: {
    label: "Repeated measures ANOVA",
    group: "ANOVA",
    params: ["f", "k", "rho", "alpha", "power", "attrition"],
    defaults: { f: 0.25, k: 3, rho: 0.50, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { f: "Effect size (Cohen's f)", k: "Number of time points", rho: "Correlation between measurements (ρ)", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { f: "Small=0.10, Medium=0.25, Large=0.40 — same as one-way ANOVA", k: "Number of repeated measurements per subject (e.g. 3 time points)", rho: "Average correlation between repeated measures — typically 0.30–0.70", attrition: "Expected dropout across the study period" },
    calculate(p) {
      // O'Brien & Kaiser (1985): variance reduction = 1 - rho for compound symmetry
      // Effective n inflated by 1/(1-rho) relative to independent groups ANOVA
      // Iterate: for each candidate n, compute noncentrality and check power
      const varianceReduction = 1 - p.rho;
      let n = 3;
      for (let iter = 0; iter < 2000; iter++) {
        // lambda = f² * n * k / varianceReduction  (within-subject NCP)
        const lambda = (p.f * p.f * n * p.k) / varianceReduction;
        const df_effect = p.k - 1;
        const pow = jStat.ncchiPowerExact(lambda, df_effect, p.alpha);
        if (pow >= p.power) break;
        n++;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return {
        n1: n_adj,
        n2: null,
        total: n_adj,
        perGroup: n_adj,
        effectSize: p.f,
        effectLabel: "Cohen's f",
        note: `n = ${n_adj} subjects, each measured ${p.k} times (total observations = ${n_adj * p.k})`,
      };
    },
  },

  // Two-way ANOVA (Aberson 2019 Ch.5; Zhang Ch.6; Stata PSS-2 p.372)
  two_way_anova: {
    label: "Two-way ANOVA",
    group: "ANOVA",
    params: ["fa", "fb", "fab", "a", "b", "alpha", "power", "attrition"],
    defaults: { fa: 0.25, fb: 0.25, fab: 0.10, a: 2, b: 3, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { fa: "Effect size — Factor A (f)", fb: "Effect size — Factor B (f)", fab: "Effect size — Interaction A×B (f)", a: "Levels of Factor A", b: "Levels of Factor B", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { fa: "Cohen's f for main effect A — small=0.10, medium=0.25, large=0.40", fb: "Cohen's f for main effect B", fab: "Interaction effect — often smaller; start with 0.10", a: "Number of levels for factor A (e.g. 2 groups)", b: "Number of levels for factor B (e.g. 3 time points or doses)" },
    calculate(p) {
      const k = p.a * p.b; // total cells
      // Calculate n needed for each effect separately, take maximum
      const calcN = (f, df_effect) => {
        let n = 3;
        for (let iter = 0; iter < 2000; iter++) {
          const N = n * k;
          const lambda = f * f * N;
          const pow = jStat.ncchiPowerExact(lambda, df_effect, p.alpha);
          if (pow >= p.power) break;
          n++;
        }
        return n;
      };
      const nA  = calcN(p.fa,  p.a - 1);
      const nB  = calcN(p.fb,  p.b - 1);
      const nAB = calcN(p.fab, (p.a - 1) * (p.b - 1));
      const n = Math.max(nA, nB, nAB);
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return {
        n1: n_adj,
        n2: null,
        total: n_adj * k,
        perGroup: n_adj,
        effectSize: Math.max(p.fa, p.fb, p.fab),
        effectLabel: "Max Cohen's f",
        note: `${n_adj} per cell × ${k} cells = ${n_adj * k} total. Powered for: A (f=${p.fa}), B (f=${p.fb}), A×B (f=${p.fab}).`,
      };
    },
  },

  // ANCOVA (Cohen 1988 Ch.9; Aberson 2019 Ch.9; Chow 2008 p.39)
  ancova: {
    label: "ANCOVA (covariate-adjusted)",
    group: "ANOVA",
    params: ["f", "k", "rho_cov", "alpha", "power", "attrition"],
    defaults: { f: 0.25, k: 2, rho_cov: 0.50, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { f: "Effect size (Cohen's f)", k: "Number of groups", rho_cov: "Outcome–covariate correlation (ρ)", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { f: "Target effect size — ANCOVA increases power so smaller n needed vs. ANOVA", k: "Number of independent groups", rho_cov: "Correlation between outcome and covariate (e.g. baseline). Higher ρ = greater power gain. Typically 0.30–0.70 for baseline as covariate." },
    calculate(p) {
      const varianceReduction = 1 - p.rho_cov * p.rho_cov;
      let n = 3;
      for (let iter = 0; iter < 2000; iter++) {
        const N = n * p.k;
        // ANCOVA lambda: adjust f upward by 1/sqrt(1-rho²)
        const lambda = (p.f * p.f * N) / varianceReduction;
        const pow = jStat.ncchiPowerExact(lambda, p.k - 1, p.alpha);
        if (pow >= p.power) break;
        n++;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return {
        n1: n_adj,
        n2: null,
        total: n_adj * p.k,
        perGroup: n_adj,
        effectSize: p.f,
        effectLabel: "Cohen's f",
        note: `Variance reduction factor: ${varianceReduction.toFixed(3)} (ρ²=${(p.rho_cov**2).toFixed(3)}). Compare: ANOVA would need ~${Math.ceil(n_adj / varianceReduction)} per group.`,
      };
    },
  },

  // Paired proportions / McNemar (Machin 2009 Ch.7; Chow 2008 Ch.7; Stata PSS-2 p.273)
  paired_proportions: {
    label: "Paired proportions (McNemar)",
    group: "Proportions",
    params: ["p01", "p10", "alpha", "power", "sides", "attrition"],
    defaults: { p01: 0.10, p10: 0.25, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { p01: "Discordant proportion p₀₁ (−/+)", p10: "Discordant proportion p₁₀ (+/−)", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { p01: "Pairs where condition 1 negative, condition 2 positive", p10: "Pairs where condition 1 positive, condition 2 negative — usually larger if intervention works", attrition: "Inflate pairs to allow for dropout" },
    calculate(p) {
      const pdis = p.p01 + p.p10;
      const diff = Math.abs(p.p10 - p.p01);
      // Connett et al. (1987) formula for McNemar
      const z_a = za(p.alpha, p.sides);
      const z_b = zb(p.power);
      const n_raw = Math.pow(z_a * Math.sqrt(pdis) + z_b * Math.sqrt(pdis - diff * diff), 2) / (diff * diff);
      const n = Math.ceil(n_raw / (1 - p.attrition / 100));
      return { n1: n, n2: n, total: n, effectSize: null, effectLabel: null, note: `n = ${n} matched pairs required` };
    },
  },

  // Chi-square (Cohen 1988 Ch.7; Aberson 2019 Ch.2; Chow 2008 Ch.7)
  chi_square: {
    label: "Chi-square (goodness-of-fit / independence)",
    group: "Proportions",
    params: ["w", "df", "alpha", "power", "attrition"],
    defaults: { w: 0.30, df: 1, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { w: "Effect size (Cohen's w)", df: "Degrees of freedom", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { w: "Small=0.10, Medium=0.30, Large=0.50 (Cohen 1988, Ch.7)", df: "(rows−1)×(cols−1) for contingency table; k−1 for goodness-of-fit" },
    calculate(p) {
      const zApprox = jStat.normalInv(1 - p.alpha) + jStat.normalInv(p.power);
      let n = Math.max(5, Math.floor(zApprox * zApprox / (p.w * p.w)));
      while (true) {
        if (jStat.ncchiPowerExact(p.w * p.w * n, p.df, p.alpha) >= p.power) break;
        n++;
        if (n > 200000) break;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj, effectSize: p.w, effectLabel: "Cohen's w" };
    },
  },

  // Two correlations (Cohen 1988 Ch.4; Aberson 2019 Ch.4; Stata PSS-2 p.337)
  two_correlations: {
    label: "Comparing two correlations",
    group: "Regression",
    params: ["rho1", "rho2", "alpha", "power", "sides", "attrition"],
    defaults: { rho1: 0.30, rho2: 0.55, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { rho1: "Correlation group 1 (ρ₁)", rho2: "Correlation group 2 (ρ₂)", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { rho1: "Expected correlation in group 1 (e.g. control)", rho2: "Expected correlation in group 2 (e.g. intervention)" },
    calculate(p) {
      const z1 = 0.5 * Math.log((1 + p.rho1) / (1 - p.rho1));
      const z2 = 0.5 * Math.log((1 + p.rho2) / (1 - p.rho2));
      // n per group (each group needs 3 more for z-transform correction)
      const n_raw = Math.pow((za(p.alpha, p.sides) + zb(p.power)) / (z1 - z2), 2) + 3;
      const n = Math.ceil(n_raw / (1 - p.attrition / 100));
      return { n1: n, n2: n, total: n * 2, effectSize: Math.abs(p.rho1 - p.rho2), effectLabel: "|ρ₁ − ρ₂|" };
    },
  },

  // Simple regression slope (Cohen 1988 Ch.9; Stata PSS-2 p.421)
  simple_slope: {
    label: "Simple linear regression (slope)",
    group: "Regression",
    params: ["rho", "alpha", "power", "sides", "attrition"],
    defaults: { rho: 0.30, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { rho: "Expected correlation (ρ = standardised slope)", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { rho: "For simple regression, power equals that of testing the correlation. Small=0.10, Medium=0.30, Large=0.50." },
    calculate(p) {
      const z_rho = 0.5 * Math.log((1 + p.rho) / (1 - p.rho));
      const n = Math.ceil(Math.pow((za(p.alpha, p.sides) + zb(p.power)) / z_rho, 2) + 2.2);
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj, effectSize: p.rho, effectLabel: "ρ (= β standardised)" };
    },
  },

  // Logistic regression (Zhang 2018 Ch.10; Aberson 2019 Ch.9; SAS PROC POWER LOGISTIC)
  logistic_regression: {
    label: "Logistic regression",
    group: "Regression",
    params: ["or", "p0", "r2", "alpha", "power", "sides", "attrition"],
    defaults: { or: 2.0, p0: 0.30, r2: 0.10, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { or: "Odds ratio (OR) for predictor", p0: "Baseline event probability", r2: "R² with other covariates", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { or: "OR for the predictor of interest — e.g. 2.0 = doubled odds per unit increase", p0: "Expected event rate in the sample (between 0.10 and 0.90 for good sensitivity)", r2: "Variance explained by other covariates in the model. Set to 0 if no covariates." },
    calculate(p) {
      // Hsieh, Bloch & Larsen (1998) Stat Med approximation
      const p1 = (p.or * p.p0) / (1 - p.p0 + p.or * p.p0);
      const pbar = (p.p0 + p1) / 2;
      const v = pbar * (1 - pbar);
      const n = Math.ceil(
        Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) /
        (v * Math.pow(Math.log(p.or), 2) * (1 - p.r2)) /
        (1 - p.attrition / 100)
      );
      return { n1: n, n2: null, total: n, effectSize: p.or, effectLabel: "Odds ratio" };
    },
  },

  // Non-inferiority proportions (Machin 2009 Ch.9; Chow 2008 Ch.5; Stata PSS-2 p.232 margin option)
  non_inferiority_prop: {
    label: "Non-inferiority — proportions",
    group: "Proportions",
    params: ["p1", "p2", "delta", "alpha", "power", "attrition"],
    defaults: { p1: 0.70, p2: 0.70, delta: 0.10, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { p1: "Reference proportion (control)", p2: "Expected test proportion", delta: "Non-inferiority margin (δ)", alpha: "Alpha (α) — one-sided", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { delta: "Maximum acceptable inferiority — e.g. 0.10 means test may be up to 10% worse", p2: "Set equal to p1 for conservative estimate; higher if test expected better", alpha: "Typically 0.025 (one-sided) for NI trials per ICH E9" },
    calculate(p) {
      const z_a = jStat.normalInv(1 - p.alpha); // one-sided
      const z_b = zb(p.power);
      const pbar = (p.p1 + p.p2) / 2;
      const n = Math.ceil(
        Math.pow(z_a * Math.sqrt(2 * pbar * (1 - pbar)) + z_b * Math.sqrt(p.p1 * (1 - p.p1) + p.p2 * (1 - p.p2)), 2) /
        Math.pow(p.delta - (p.p1 - p.p2), 2) /
        (1 - p.attrition / 100)
      );
      return { n1: n, n2: n, total: n * 2, effectSize: p.delta, effectLabel: "NI margin (δ)" };
    },
  },

  // Equivalence means TOST (Machin 2009 Ch.9; Chow 2008 Ch.3; Stata PSS-2 p.144 margin option)
  equivalence_means: {
    label: "Equivalence (TOST) — means",
    group: "Clinical Trials",
    params: ["delta", "mu_diff", "sd", "alpha", "power", "attrition"],
    defaults: { delta: 5, mu_diff: 0, sd: 10, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { delta: "Equivalence margin (±δ)", mu_diff: "Expected true difference (μ₁−μ₂)", sd: "Common SD", alpha: "Alpha (α) — one-sided", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { delta: "Maximum acceptable difference to claim equivalence — set from clinical judgement", mu_diff: "Usually 0 for conservative estimate; positive if test expected slightly better", alpha: "Typically 0.05 one-sided (= 0.10 two-sided TOST) per regulatory guidance" },
    calculate(p) {
      const z_a = jStat.normalInv(1 - p.alpha);
      const z_b = zb(p.power);
      const n = Math.ceil(
        Math.pow(z_a + z_b, 2) * 2 * p.sd * p.sd /
        Math.pow(p.delta - Math.abs(p.mu_diff), 2) /
        (1 - p.attrition / 100)
      );
      return { n1: n, n2: n, total: n * 2, effectSize: (p.delta / p.sd).toFixed(4), effectLabel: "δ / SD" };
    },
  },

  // Log-rank cluster RCT (Zhang 2018 Ch.12; Stata PSS-2 p.601)
  logrank_cluster: {
    label: "Cluster RCT — log-rank (survival)",
    group: "Survival",
    params: ["hr", "p_event", "icc", "m", "alpha", "power", "sides", "attrition"],
    defaults: { hr: 0.70, p_event: 0.60, icc: 0.05, m: 20, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { hr: "Hazard ratio (HR)", p_event: "Overall event probability", icc: "ICC (ρ)", m: "Cluster size (m)", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { hr: "HR < 1 = intervention reduces hazard", icc: "Intra-cluster correlation for the survival event", m: "Average number of individuals per cluster" },
    calculate(p) {
      // Individual-level events needed (Schoenfeld)
      const e_ind = Math.ceil(
        Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) /
        Math.pow(Math.log(p.hr), 2) * 4
      );
      const deff = 1 + (p.m - 1) * p.icc;
      const n_ind = Math.ceil((e_ind / p.p_event) * deff / (1 - p.attrition / 100));
      const k = Math.ceil(n_ind / p.m); // clusters per arm
      return {
        n1: k, n2: k, total: k * 2 * p.m,
        events: e_ind,
        effectSize: p.hr, effectLabel: "Hazard ratio",
        note: `k = ${k} clusters/arm × ${p.m} subjects = ${k * p.m}/arm. DEFF = ${deff.toFixed(2)}.`,
      };
    },
  },

  // Cox proportional hazards (Schoenfeld 1983; Chow 2008 Ch.10; Stata PSS-2 p.519)
  cox_regression: {
    label: "Cox proportional hazards",
    group: "Survival",
    params: ["hr", "p_event", "r2", "alpha", "power", "sides", "attrition"],
    defaults: { hr: 0.70, p_event: 0.60, r2: 0.10, alpha: 0.05, power: 0.80, sides: 2, attrition: 0 },
    labels: { hr: "Hazard ratio for predictor", p_event: "Overall event probability", r2: "R² with other covariates", alpha: "Alpha (α)", power: "Power (1-β)", sides: "Sides", attrition: "Attrition %" },
    hints: { hr: "HR for the key predictor of interest", p_event: "Expected proportion experiencing the event", r2: "Variance explained by other covariates (Schoenfeld 1983). Set to 0 if none." },
    calculate(p) {
      const events = Math.ceil(
        Math.pow(za(p.alpha, p.sides) + zb(p.power), 2) /
        (Math.pow(Math.log(p.hr), 2) * (1 - p.r2))
      );
      const n = Math.ceil((events / p.p_event) / (1 - p.attrition / 100));
      return { n1: n, n2: null, total: n, events, effectSize: p.hr, effectLabel: "Hazard ratio" };
    },
  },

  // ICC reliability (Machin 2009 Ch.14; Stata — no native command)
  icc_reliability: {
    label: "ICC — reliability / repeatability",
    group: "Diagnostic",
    params: ["icc0", "icc1", "k", "alpha", "power", "attrition"],
    defaults: { icc0: 0.60, icc1: 0.80, k: 2, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { icc0: "ICC under H₀ (minimum acceptable)", icc1: "Expected ICC (ρ₁)", k: "Raters / replications per subject", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { icc0: "Below this value, reliability is considered poor — set from clinical standards", icc1: "Expected true ICC in your measurement study", k: "Number of repeated measurements or raters per subject (usually 2 or 3)" },
    calculate(p) {
      // Donner & Eliasziw (1987) approximation
      const f0 = (1 + (p.k - 1) * p.icc0) / (1 - p.icc0);
      const f1 = (1 + (p.k - 1) * p.icc1) / (1 - p.icc1);
      const lnF0 = Math.log(f0);
      const lnF1 = Math.log(f1);
      const se = Math.sqrt(2 / (p.k - 1));
      const z_diff = (lnF1 - lnF0) / se;
      let n = 5;
      for (let iter = 0; iter < 2000; iter++) {
        const se_n = Math.sqrt(2 * p.k / ((p.k - 1) * n));
        const z_stat = (lnF1 - lnF0) / se_n;
        const pow = 1 - jStat.normalCDF(za(p.alpha, 2) - z_stat);
        if (pow >= p.power) break;
        n++;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return {
        n1: n_adj, n2: null, total: n_adj,
        effectSize: (p.icc1 - p.icc0).toFixed(3), effectLabel: "ΔICC",
        note: `${n_adj} subjects × ${p.k} raters = ${n_adj * p.k} total observations`,
      };
    },
  },

  // Simon two-stage phase II (Machin 2009 Ch.16; clinfun R package)
  simon_two_stage: {
    label: "Simon two-stage (Phase II oncology)",
    group: "Phase II",
    params: ["p0", "p1", "alpha", "beta"],
    defaults: { p0: 0.20, p1: 0.40, alpha: 0.10, beta: 0.10 },
    labels: { p0: "Unacceptable response rate (p₀)", p1: "Target response rate (p₁)", alpha: "Type I error (α)", beta: "Type II error (β = 1−power)" },
    hints: { p0: "Response rate that would make the drug unworthy of further study", p1: "Response rate that justifies further development", alpha: "Usually 0.05–0.10 for phase II", beta: "Usually 0.10–0.20 (power = 1−β = 80–90%)" },
    calculate(p) {
      const binomPMF = (k, n, prob) => {
        let c = 1;
        for (let j = 0; j < k; j++) c = c * (n - j) / (j + 1);
        return c * Math.pow(prob, k) * Math.pow(1 - prob, n - k);
      };
      const binomCDF = (k, n, prob) => {
        let s = 0;
        for (let i = 0; i <= k; i++) s += binomPMF(i, n, prob);
        return s;
      };
      let best = { en: Infinity, n1: 0, r1: 0, n: 0, r: 0, alpha_act: 1, power_act: 0 };
      for (let n1 = 4; n1 <= 35; n1++) {
        for (let r1 = 0; r1 < n1; r1++) {
          const pet0 = binomCDF(r1, n1, p.p0);
          for (let n2 = 1; n2 <= 55; n2++) {
            const ntot = n1 + n2;
            for (let r = r1 + 1; r <= ntot; r++) {
              let alpha_act = 0, power_act = 0;
              for (let x = r1 + 1; x <= n1; x++) {
                const px0 = binomPMF(x, n1, p.p0);
                const px1 = binomPMF(x, n1, p.p1);
                const p2crit = r - x;
                alpha_act += px0 * (1 - binomCDF(p2crit, n2, p.p0));
                power_act += px1 * (1 - binomCDF(p2crit, n2, p.p1));
              }
              if (alpha_act > p.alpha || (1 - power_act) > p.beta) continue;
              const en = n1 + (1 - pet0) * n2;
              if (en < best.en) {
                best = { en, n1, r1, n: ntot, r, alpha_act: +alpha_act.toFixed(4), power_act: +power_act.toFixed(4) };
              }
            }
          }
        }
      }
      if (best.en === Infinity) {
        return { n1: null, n2: null, total: null, effectSize: null, effectLabel: null, note: "No feasible design found — try relaxing α or β." };
      }
      return {
        n1: best.n1, n2: best.n - best.n1, total: best.n,
        effectSize: (p.p1 - p.p0).toFixed(3), effectLabel: "Δp (p₁ − p₀)",
        note: `Stage 1: n=${best.n1}, stop early if ≤${best.r1} responses. Stage 2: total n=${best.n}, reject drug if ≤${best.r} total responses. E[n|H₀]=${best.en.toFixed(1)}, actual α=${best.alpha_act}, actual power=${best.power_act}.`,
      };
    },
  },

  // Survey mean (Valliant 2013 Ch.3; Aday 2006 Ch.7; Stata PSS-3 p.702)
  survey_mean: {
    label: "Survey — estimating a mean",
    group: "Survey",
    params: ["sd", "d", "alpha", "deff", "response_rate"],
    defaults: { sd: 15, d: 2, alpha: 0.05, deff: 1.0, response_rate: 100 },
    labels: { sd: "Population SD", d: "Margin of error (±d)", alpha: "Alpha (α)", deff: "Design effect (DEFF)", response_rate: "Expected response rate %" },
    hints: { sd: "Estimated from pilot data or published literature", d: "Acceptable half-width of CI on the mean scale", deff: "1.0 = simple random sample; multiply for cluster designs", response_rate: "Adjust for expected non-response" },
    calculate(p) {
      const z = za(p.alpha, 2);
      const n = Math.ceil(Math.ceil(z * z * p.sd * p.sd / (p.d * p.d)) * p.deff / (p.response_rate / 100));
      return { n1: n, n2: null, total: n, effectSize: null, effectLabel: null };
    },
  },

  // Mediation Sobel (Zhang 2018 Ch.14)
  mediation_sobel: {
    label: "Simple mediation (Sobel test)",
    group: "SEM",
    params: ["a", "b", "sd_a", "sd_b", "alpha", "power", "attrition"],
    defaults: { a: 0.40, b: 0.30, sd_a: 1.0, sd_b: 1.0, alpha: 0.05, power: 0.80, attrition: 0 },
    labels: { a: "Path a coefficient (X → M)", b: "Path b coefficient (M → Y)", sd_a: "SD of predictor X", sd_b: "SD of mediator M", alpha: "Alpha (α)", power: "Power (1-β)", attrition: "Attrition %" },
    hints: { a: "Unstandardised coefficient from X to mediator M", b: "Unstandardised coefficient from M to outcome Y (controlling X)", sd_a: "Standard deviation of the predictor variable", sd_b: "Standard deviation of the mediator variable" },
    calculate(p) {
      const ab = Math.abs(p.a * p.b);
      // Sobel SE = sqrt(b²·sa² + a²·sb²)  (per unit, before scaling by 1/sqrt(n))
      const se_unit = Math.sqrt(p.b * p.b * p.sd_a * p.sd_a + p.a * p.a * p.sd_b * p.sd_b);
      let n = 20;
      for (let iter = 0; iter < 3000; iter++) {
        const se = se_unit / Math.sqrt(n);
        const pow = 1 - jStat.normalCDF(za(p.alpha, 2) - ab / se);
        if (pow >= p.power) break;
        n++;
      }
      const n_adj = Math.ceil(n / (1 - p.attrition / 100));
      return { n1: n_adj, n2: null, total: n_adj, effectSize: ab.toFixed(4), effectLabel: "Indirect effect (ab)" };
    },
  },
};

// ─────────────────────────────────────────────
// CITATION DATABASE
// ─────────────────────────────────────────────
const CITATIONS = {
  two_means: [
    { ref: "Cohen (1988)", detail: "Ch. 2, p. 19–66", note: "t-test for means, effect size d" },
    { ref: "Aberson (2019)", detail: "Ch. 3, p. 34–53", note: "Applied R-based examples" },
    { ref: "Machin et al. (2009)", detail: "Ch. 5, p. 47–57", note: "Clinical sample size tables" },
    { ref: "Chow et al. (2008)", detail: "Ch. 3, p. 47–63", note: "FDA-aligned regulatory context" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 4, p. 49–68", note: "WebPower implementation" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 144", note: "power twomeans" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, TWOSAMPLEMEANS", note: "proc power; twosamplemeans" },
  ],
  paired_means: [
    { ref: "Cohen (1988)", detail: "Ch. 2, p. 48–51", note: "Paired t-test" },
    { ref: "Machin et al. (2009)", detail: "Ch. 7, p. 67–83", note: "Paired continuous data" },
    { ref: "Chow et al. (2008)", detail: "Ch. 4, p. 64–67", note: "Crossover / paired design" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 182", note: "power pairedmeans" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, PAIREDMEANS", note: "proc power; pairedmeans" },
  ],
  one_mean: [
    { ref: "Cohen (1988)", detail: "Ch. 2, p. 19–26", note: "One-sample t-test" },
    { ref: "Machin et al. (2009)", detail: "Ch. 5, p. 57", note: "One-sample table" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 114", note: "power onemean" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, ONESAMPLEMEANS", note: "proc power; onesamplemeans" },
  ],
  two_proportions: [
    { ref: "Cohen (1988)", detail: "Ch. 7, p. 179–213", note: "Arcsine transformation, effect size h" },
    { ref: "Fleiss, Levin & Paik (2003)", detail: "Ch. 4, p. 85–86", note: "Equation 4.19 — gold standard formula" },
    { ref: "Aberson (2019)", detail: "Ch. 2, p. 28–33", note: "Independent proportions with R" },
    { ref: "Machin et al. (2009)", detail: "Ch. 3, p. 30–41", note: "Clinical table with OR" },
    { ref: "Chow et al. (2008)", detail: "Ch. 5, p. 82–100", note: "Superiority & NI framing" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 3, p. 38–47", note: "wp.prop() implementation" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 232", note: "power twoproportions" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, TWOSAMPLEFREQ", note: "proc power; twosamplefreq" },
  ],
  one_proportion: [
    { ref: "Aberson (2019)", detail: "Ch. 2, p. 29–31", note: "One-sample proportion" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 3, p. 29–37", note: "wp.prop() single sample" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 198", note: "power oneproportion" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, ONESAMPLEFREQ", note: "proc power; onesamplefreq" },
  ],
  one_way_anova: [
    { ref: "Cohen (1988)", detail: "Ch. 8, p. 273–406", note: "One-way ANOVA, effect size f" },
    { ref: "Aberson (2019)", detail: "Ch. 5, p. 69–87", note: "ANOVA power with contrasts" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 5, p. 69–88", note: "wp.anova() implementation" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 351", note: "power oneway" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, ONEWAYANOVA", note: "proc power; onewayanova" },
  ],
  correlation: [
    { ref: "Cohen (1988)", detail: "Ch. 3, p. 75–107", note: "Pearson r, Fisher's z-transform" },
    { ref: "Aberson (2019)", detail: "Ch. 4, p. 54–68", note: "Correlation power with R" },
    { ref: "Machin et al. (2009)", detail: "Ch. 12, p. 153–157", note: "Correlation coefficient table" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 2, p. 17–28", note: "wp.correlation() implementation" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 324", note: "power onecorrelation" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, ONECORR", note: "proc power; onecorr" },
  ],
  regression_r2: [
    { ref: "Cohen (1988)", detail: "Ch. 9, p. 407–468", note: "Multiple regression R², effect size f²" },
    { ref: "Aberson (2019)", detail: "Ch. 7, p. 99–119", note: "Regression power with R" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 9, p. 151–162", note: "wp.regression() implementation" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 435", note: "power rsquared" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, MULTREG", note: "proc power; multreg" },
  ],
  logrank: [
    { ref: "Schoenfeld (1981)", detail: "Biometrika 68(1):316–319", note: "Original log-rank formula" },
    { ref: "Lachin & Foulkes (1986)", detail: "Biometrics 42:93–109", note: "Extended survival formula" },
    { ref: "Machin et al. (2009)", detail: "Ch. 8, p. 84–101", note: "Survival tables — log-rank" },
    { ref: "Chow et al. (2008)", detail: "Ch. 10, p. 177–181", note: "Log-rank regulatory context" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 575", note: "power logrank" },
  ],
  survey_proportion: [
    { ref: "Valliant, Dever & Kreuter (2013)", detail: "Ch. 3–5", note: "Survey design with DEFF" },
    { ref: "Aday & Cornelius (2006)", detail: "Ch. 7, p. 154–193", note: "Health survey sample size" },
    { ref: "Lwanga & Lemeshow / WHO (1991)", detail: "Section 2", note: "Classic WHO formula" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-3, p. 619", note: "ciwidth (precision approach)" },
  ],
  sem_rmsea: [
    { ref: "MacCallum, Browne & Sugawara (1996)", detail: "Psychol Methods 1(2):130–149", note: "RMSEA-based power — foundational paper" },
    { ref: "Jak et al. (2020)", detail: "Behav Res Methods, power4SEM app", note: "χ² and RMSEA methods, tutorial" },
    { ref: "Moshagen & Bader (2024)", detail: "Behav Res Methods 56:2901–2922", note: "semPower 2 — most current R package" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 15, p. 249–270", note: "wp.sem.rmsea() implementation" },
    { ref: "Jobst, Bader & Moshagen (2023)", detail: "Psychol Methods 28:207–221", note: "SEM power tutorial" },
  ],
  repeated_measures: [
    { ref: "Aberson (2019)", detail: "Ch. 6, p. 88–97", note: "Within-subjects ANOVA power in R" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 8, p. 117–136", note: "wp.rmanova() WebPower implementation" },
    { ref: "Machin et al. (2009)", detail: "Ch. 6, p. 58–66", note: "Repeated measures multiplying factor" },
    { ref: "Chow et al. (2008)", detail: "Ch. 3, p. 302", note: "Repeated measures ANOVA — regulatory context" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 393", note: "power repeated" },
  ],
  two_way_anova: [
    { ref: "Aberson (2019)", detail: "Ch. 5, p. 78–87", note: "Two-factor ANOVA with interactions in R" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 6, p. 85–100", note: "wp.anova() two-way" },
    { ref: "Cohen (1988)", detail: "Ch. 8, p. 273–406", note: "Factorial ANOVA effect sizes" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 372", note: "power twoway" },
  ],
  ancova: [
    { ref: "Cohen (1988)", detail: "Ch. 9, p. 407–468", note: "ANCOVA as special case of regression" },
    { ref: "Aberson (2019)", detail: "Ch. 9, p. 120–140", note: "ANCOVA power analysis in R" },
    { ref: "Chow et al. (2008)", detail: "Ch. 3, p. 39", note: "Adjustment for covariates — regulatory context" },
    { ref: "Machin et al. (2009)", detail: "Ch. 5, p. 47–57", note: "Baseline adjustment in clinical trials" },
  ],
  paired_proportions: [
    { ref: "Machin et al. (2009)", detail: "Ch. 7, p. 82", note: "Sample sizes for paired binary data" },
    { ref: "Chow et al. (2008)", detail: "Ch. 7, p. 143–155", note: "McNemar's test — regulatory context" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 273", note: "power pairedproportions" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, PAIREDFREQ", note: "proc power; pairedfreq" },
  ],
  chi_square: [
    { ref: "Cohen (1988)", detail: "Ch. 7, p. 215–271", note: "Chi-square tests, effect size w — foundational" },
    { ref: "Aberson (2019)", detail: "Ch. 2, p. 18–33", note: "Chi-square with pwr.chisq.test() in R" },
    { ref: "Chow et al. (2008)", detail: "Ch. 7, p. 143–155", note: "Contingency tables — clinical context" },
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 3, p. 29–38", note: "wp.chisq() WebPower implementation" },
  ],
  two_correlations: [
    { ref: "Cohen (1988)", detail: "Ch. 4, p. 109–131", note: "Comparing two independent correlations" },
    { ref: "Aberson (2019)", detail: "Ch. 4, p. 57–68", note: "Correlation difference power in R" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 337", note: "power twocorrelations" },
  ],
  simple_slope: [
    { ref: "Cohen (1988)", detail: "Ch. 9, p. 407–468", note: "Simple regression as special case of multiple regression" },
    { ref: "Aberson (2019)", detail: "Ch. 7, p. 99–119", note: "Regression slope power in R" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 421", note: "power oneslope" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, MULTREG", note: "proc power; multreg (u=1)" },
  ],
  logistic_regression: [
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 10, p. 163–172", note: "wp.logistic() WebPower implementation" },
    { ref: "Aberson (2019)", detail: "Ch. 9, p. 120–140", note: "Logistic regression power in R" },
    { ref: "SAS PROC POWER (2015)", detail: "Ch. 89, LOGISTIC", note: "proc power; logistic — native SAS command" },
    { ref: "Hsieh, Bloch & Larsen (1998)", detail: "Stat Med 17(13):1503–1506", note: "Underlying approximation formula" },
  ],
  non_inferiority_prop: [
    { ref: "Machin et al. (2009)", detail: "Ch. 9, p. 120–141", note: "NI and equivalence — sample size tables" },
    { ref: "Chow et al. (2008)", detail: "Ch. 5, p. 87–100", note: "Non-inferiority proportions — FDA context" },
    { ref: "Fleiss, Levin & Paik (2003)", detail: "Ch. 4, p. 85–86", note: "Proportions formula basis" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 232 (margin option)", note: "power twoproportions, margin()" },
  ],
  equivalence_means: [
    { ref: "Machin et al. (2009)", detail: "Ch. 9, p. 102–121", note: "Equivalence tables for means — TOST" },
    { ref: "Chow et al. (2008)", detail: "Ch. 3, p. 59–63", note: "TOST for means — bioequivalence context" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 144 (margin option)", note: "power twomeans, margin()" },
  ],
  logrank_cluster: [
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 12, p. 183–210", note: "wp.crt2arm() cluster survival" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 601", note: "power logrank, cluster" },
    { ref: "Schoenfeld (1981)", detail: "Biometrika 68(1):316–319", note: "Underlying log-rank events formula" },
  ],
  cox_regression: [
    { ref: "Schoenfeld (1983)", detail: "Biometrics 39:499–503", note: "Cox regression power — foundational formula" },
    { ref: "Chow et al. (2008)", detail: "Ch. 10, p. 172–177", note: "Cox model — regulatory context" },
    { ref: "Machin et al. (2009)", detail: "Ch. 8, p. 97–101", note: "Cox model in clinical trials" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-2, p. 519", note: "power cox" },
  ],
  icc_reliability: [
    { ref: "Machin et al. (2009)", detail: "Ch. 14, p. 193–194", note: "ICC sample size tables" },
    { ref: "Donner & Eliasziw (1987)", detail: "Stat Med 6(4):431–441", note: "ICC power formula" },
  ],
  simon_two_stage: [
    { ref: "Machin et al. (2009)", detail: "Ch. 16, p. 207–235", note: "Simon designs — tables, software and examples" },
    { ref: "Simon (1989)", detail: "Control Clin Trials 10(1):1–10", note: "Original Simon two-stage paper" },
  ],
  survey_mean: [
    { ref: "Valliant, Dever & Kreuter (2013)", detail: "Ch. 3, p. 75–120", note: "Mean estimation in complex surveys" },
    { ref: "Aday & Cornelius (2006)", detail: "Ch. 7, p. 154–193", note: "Health survey sample size for means" },
    { ref: "Stata PSS Manual (2023)", detail: "PSS-3, p. 702", note: "ciwidth onemean" },
  ],
  mediation_sobel: [
    { ref: "Zhang & Yuan (2018)", detail: "Ch. 14, p. 239–248", note: "wp.mediation() — Sobel test implementation" },
    { ref: "Fritz & MacKinnon (2007)", detail: "Psychol Sci 18(3):233–239", note: "Required sample sizes for mediation" },
  ],
};

// ─────────────────────────────────────────────
// R / STATA / SAS CODE GENERATORS
// ─────────────────────────────────────────────
const CODE = {
  two_means: (p, r) => ({
    r: `# Two independent means — power.t.test (base R)
power.t.test(
  n       = NULL,        # solve for n
  delta   = 1,           # set delta=1 and sd=1/d to express as Cohen's d
  sd      = ${(1 / r.effectSize).toFixed(4)},
  sig.level = ${p.alpha},
  power   = ${p.power},
  type    = "two.sample",
  alternative = "${p.sides === 2 ? "two.sided" : "one.sided"}"
)
# Result: n = ${r.n1} per group (Cohen's d = ${r.effectSize.toFixed(4)})

# Alternatively using pwr package:
# install.packages("pwr")
library(pwr)
pwr.t.test(d = ${r.effectSize.toFixed(4)}, sig.level = ${p.alpha},
           power = ${p.power}, type = "two.sample")`,
    stata: `* Two independent means (Cohen's d = ${r.effectSize.toFixed(4)})
power twomeans 0 ${r.effectSize.toFixed(4)}, ///
    sd(1) ///
    alpha(${p.alpha}) ///
    power(${p.power}) ///
    ${p.sides === 1 ? "onesided" : ""}
* Reference: Stata PSS Manual (2023), PSS-2 p.144`,
    sas: `/* Two independent means (Cohen's d = ${r.effectSize.toFixed(4)}) */
proc power;
  twosamplemeans
    meandiff = ${r.effectSize.toFixed(4)}
    stddev   = 1
    alpha    = ${p.alpha}
    power    = ${p.power}
    sides    = ${p.sides}
    ntotal   = .;  /* solve for N */
run;
/* Reference: SAS/STAT 14.1, Ch.89 TWOSAMPLEMEANS */`,
  }),
  paired_means: (p, r) => ({
    r: `# Paired means — power.t.test (base R)
power.t.test(
  n       = NULL,
  delta   = ${p.mu_diff},
  sd      = ${p.sd_diff},
  sig.level = ${p.alpha},
  power   = ${p.power},
  type    = "paired",
  alternative = "${p.sides === 2 ? "two.sided" : "one.sided"}"
)
# Result: n = ${r.total} pairs`,
    stata: `* Paired means
power pairedmeans ${p.mu_diff}, ///
    sddiff(${p.sd_diff}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.182`,
    sas: `/* Paired means */
proc power;
  pairedmeans
    meandiff = ${p.mu_diff}
    std      = ${p.sd_diff}
    alpha    = ${p.alpha}
    power    = ${p.power}
    npairs   = .;
run;`,
  }),
  one_mean: (p, r) => ({
    r: `# One-sample mean
power.t.test(
  n       = NULL,
  delta   = ${Math.abs(p.mu1 - p.mu0)},
  sd      = ${p.sd},
  sig.level = ${p.alpha},
  power   = ${p.power},
  type    = "one.sample"
)
# Result: n = ${r.total}`,
    stata: `* One-sample mean
power onemean ${p.mu0} ${p.mu1}, ///
    sd(${p.sd}) alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.114`,
    sas: `/* One-sample mean */
proc power;
  onesamplemeans
    mean  = ${p.mu1}
    nullm = ${p.mu0}
    std   = ${p.sd}
    alpha = ${p.alpha}
    power = ${p.power}
    ntotal = .;
run;`,
  }),
  two_proportions: (p, r) => ({
    r: `# Two independent proportions
# Method 1: base R
power.prop.test(
  n  = NULL,
  p1 = ${p.p1},
  p2 = ${p.p2},
  sig.level = ${p.alpha},
  power     = ${p.power},
  alternative = "${p.sides === 2 ? "two.sided" : "one.sided"}"
)

# Method 2: pwr package (arcsine transformation — Cohen 1988)
library(pwr)
pwr.2p.test(
  h         = ES.h(${p.p1}, ${p.p2}),
  sig.level = ${p.alpha},
  power     = ${p.power}
)
# Result: n = ${r.n1} per group`,
    stata: `* Two independent proportions
power twoproportions ${p.p1} ${p.p2}, ///
    alpha(${p.alpha}) power(${p.power}) ///
    ${p.sides === 1 ? "onesided" : ""}
* Reference: Stata PSS Manual (2023), PSS-2 p.232`,
    sas: `/* Two proportions */
proc power;
  twosamplefreq
    refproportion = ${p.p1}
    proportion    = ${p.p2}
    alpha         = ${p.alpha}
    power         = ${p.power}
    sides         = ${p.sides}
    ntotal        = .;
run;
/* Reference: SAS/STAT 14.1, Ch.89 TWOSAMPLEFREQ */`,
  }),
  one_proportion: (p, r) => ({
    r: `# One-sample proportion
library(pwr)
pwr.p.test(
  h         = ES.h(${p.p1}, ${p.p0}),
  sig.level = ${p.alpha},
  power     = ${p.power}
)
# Result: n = ${r.total}`,
    stata: `* One-sample proportion
power oneproportion ${p.p0} ${p.p1}, ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.198`,
    sas: `/* One-sample proportion */
proc power;
  onesamplefreq
    nullproportion = ${p.p0}
    proportion     = ${p.p1}
    alpha          = ${p.alpha}
    power          = ${p.power}
    ntotal         = .;
run;`,
  }),
  one_way_anova: (p, r) => ({
    r: `# One-way ANOVA
library(pwr)
pwr.anova.test(
  k         = ${p.k},
  f         = ${p.f},
  sig.level = ${p.alpha},
  power     = ${p.power}
)
# Result: n = ${r.perGroup} per group (${r.total} total)

# WebPower alternative:
library(WebPower)
wp.anova(k = ${p.k}, f = ${p.f}, alpha = ${p.alpha}, power = ${p.power})`,
    stata: `* One-way ANOVA
power oneway, ///
    k(${p.k}) f(${p.f}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.351`,
    sas: `/* One-way ANOVA */
proc power;
  onewayanova
    groupmeans = . /* specify group means or */
    stddev     = . /* and SD */
    alpha      = ${p.alpha}
    power      = ${p.power}
    npergroup  = .;
run;`,
  }),
  correlation: (p, r) => ({
    r: `# Pearson correlation (Fisher's z-transformation)
library(pwr)
pwr.r.test(
  r         = ${p.rho},
  sig.level = ${p.alpha},
  power     = ${p.power},
  alternative = "${p.sides === 2 ? "two.sided" : "one.sided"}"
)
# Result: n = ${r.total}

# WebPower:
library(WebPower)
wp.correlation(r = ${p.rho}, alpha = ${p.alpha}, power = ${p.power})`,
    stata: `* Pearson correlation
power onecorrelation ${p.rho}, ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.324`,
    sas: `/* Pearson correlation */
proc power;
  onecorr
    corr  = ${p.rho}
    alpha = ${p.alpha}
    power = ${p.power}
    ntotal = .;
run;
/* Reference: SAS/STAT 14.1, Ch.89 ONECORR */`,
  }),
  regression_r2: (p, r) => ({
    r: `# Multiple regression — R² test
library(pwr)
pwr.f2.test(
  u         = ${p.u},          # numerator df = predictors
  f2        = ${p.f2},         # Cohen's f² (R²/(1-R²))
  sig.level = ${p.alpha},
  power     = ${p.power}
)
# Result: n = ${r.total} (v = n - u - 1 = ${r.total - p.u - 1})

# WebPower:
library(WebPower)
wp.regression(p1 = ${p.u}, f2 = ${p.f2}, alpha = ${p.alpha}, power = ${p.power})`,
    stata: `* Multiple regression R²
power rsquared, ///
    u(${p.u}) f2(${p.f2}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.435`,
    sas: `/* Multiple regression */
proc power;
  multreg
    numpredictors = ${p.u}
    rsquare       = ${(p.f2 / (1 + p.f2)).toFixed(4)}
    alpha         = ${p.alpha}
    power         = ${p.power}
    ntotal        = .;
run;`,
  }),
  logrank: (p, r) => ({
    r: `# Log-rank test (Schoenfeld 1981 formula)
# install.packages("gsDesign")
library(gsDesign)
nSurv(
  lambdaC = -log(1 - ${p.p_event}),  # approximate hazard
  hr      = ${p.hr},
  sided   = ${p.sides},
  alpha   = ${p.alpha},
  beta    = ${1 - p.power}
)

# Simpler: powerSurvEpi package
library(powerSurvEpi)
ssizeCT.default(
  power    = ${p.power},
  k        = ${p.ratio},
  pE       = ${p.p_event},
  pC       = ${p.p_event},
  RR       = ${p.hr},
  alpha    = ${p.alpha}
)
# Required events: ${r.events}, Total n: ${r.total}`,
    stata: `* Log-rank test
power logrank, ///
    hratio(${p.hr}) ///
    eventprob(${p.p_event}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.575`,
    sas: `/* Log-rank / survival */
proc power;
  twosamplesurvival
    hazardratio = ${p.hr}
    accrualtime = .  /* specify recruitment period */
    followuptime = . /* specify follow-up */
    alpha       = ${p.alpha}
    power       = ${p.power}
    ntotal      = .;
run;`,
  }),
  survey_proportion: (p, r) => ({
    r: `# Survey proportion — simple random sample
# Formula: n = z² * p * (1-p) / d²  × DEFF
z     <- qnorm(1 - ${p.alpha}/2)
p     <- ${p.p}
d     <- ${p.d}
deff  <- ${p.deff}
rr    <- ${p.response_rate / 100}

n_raw   <- (z^2 * p * (1-p)) / d^2
n_deff  <- ceiling(n_raw * deff)
n_final <- ceiling(n_deff / rr)
cat("Required n:", n_final)

# PracTools package (Valliant et al. 2013):
library(PracTools)
nProp(CV0 = ${p.d}/${p.p}, pU = ${p.p})`,
    stata: `* Survey proportion (precision approach)
ciwidth oneproportion, ///
    p(${p.p}) width(${p.d * 2}) ///
    level(${(1 - p.alpha) * 100})
* Then multiply n by DEFF = ${p.deff}
* Reference: Stata PSS Manual (2023), PSS-3 p.619`,
    sas: `/* Survey — proportion */
/* SAS does not have direct survey power proc;
   use formula approach: */
data survey_n;
  z    = probit(1 - ${p.alpha}/2);
  p    = ${p.p};
  d    = ${p.d};
  deff = ${p.deff};
  rr   = ${p.response_rate / 100};
  n    = ceil(ceil((z**2 * p*(1-p)) / d**2) * deff / rr);
  put "Required n = " n;
run;`,
  }),
  sem_rmsea: (p, r) => ({
    r: `# SEM power — RMSEA-based (MacCallum et al. 1996)

# Method 1: semPower package (Moshagen & Bader 2024)
# install.packages("semPower")
library(semPower)
ap <- semPower.aPriori(
  effect        = ${p.rmsea1},
  effect.measure = "RMSEA",
  alpha         = ${p.alpha},
  power         = ${p.power},
  df            = ${p.df}
)
summary(ap)

# Method 2: WebPower (Zhang & Yuan 2018)
library(WebPower)
wp.sem.rmsea(
  df    = ${p.df},
  rmsea0 = ${p.rmsea0},
  rmsea1 = ${p.rmsea1},
  alpha  = ${p.alpha},
  power  = ${p.power}
)
# Result: N = ${r.total}

# Online Shiny app: https://sempower.shinyapps.io/sempower`,
    stata: `* SEM RMSEA power — no native Stata command
* Use R semPower or WebPower packages
* See: Moshagen & Bader (2024) Behav Res Methods 56:2901-2922
* Online: https://sempower.shinyapps.io/sempower`,
    sas: `/* SEM RMSEA power — no native SAS PROC POWER command.
   Use R package semPower (Moshagen & Bader 2024) or
   WebPower (Zhang & Yuan 2018).
   Online: https://sempower.shinyapps.io/sempower */`,
  }),

  // ── NEW DESIGNS ──────────────────────────────────────────────────────────
  repeated_measures: (p, r) => ({
    r: `# Repeated measures ANOVA
library(WebPower)
wp.rmanova(n = NULL, ng = 1, nm = ${p.k}, f = ${p.f},
           cor.rm = ${p.rho}, alpha = ${p.alpha}, power = ${p.power}, type = 1)
# Result: n = ${r.total} subjects (each measured ${p.k} times)

# Alternatively using pwr (approximate):
library(pwr)
pwr.anova.test(k = ${p.k}, f = ${p.f}, sig.level = ${p.alpha}, power = ${p.power})
# Then divide n by variance reduction factor (1 - rho) = ${(1 - p.rho).toFixed(3)}`,
    stata: `* Repeated measures ANOVA
power repeated, k(${p.k}) f(${p.f}) corr(${p.rho}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.393`,
    sas: `/* Repeated measures — approximate via GLM power */
proc power;
  onewayanova
    groupmeans = 0 | 1   /* adjust to match f = ${p.f} */
    stddev     = 1
    corr       = ${p.rho}
    alpha      = ${p.alpha}
    power      = ${p.power}
    npergroup  = .;
run;`,
  }),

  two_way_anova: (p, r) => ({
    r: `# Two-way ANOVA (${p.a} x ${p.b} factorial)
# install.packages("pwr2ppl")
library(pwr2ppl)
# Factor A (f = ${p.fa}):
pwr.2way(a = ${p.a}, b = ${p.b}, alpha = ${p.alpha}, power = ${p.power},
         f.A = ${p.fa}, f.B = ${p.fb})
# Interaction (f = ${p.fab}):
# n per cell = ${r.perGroup} (total = ${r.total})

# Verify each effect with pwr:
library(pwr)
pwr.anova.test(k = ${p.a * p.b}, f = ${Math.max(p.fa, p.fb, p.fab)},
               sig.level = ${p.alpha}, power = ${p.power})`,
    stata: `* Two-way ANOVA
power twoway, a(${p.a}) b(${p.b}) ///
    f_a(${p.fa}) f_b(${p.fb}) f_ab(${p.fab}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.372`,
    sas: `/* Two-way ANOVA — use PROC POWER with cellmeans */
proc power;
  twosamplefreq /* approximate: specify cell means for exact */
    alpha = ${p.alpha} power = ${p.power};
run;
/* Recommended: use R pwr2ppl package for factorial designs */`,
  }),

  ancova: (p, r) => ({
    r: `# ANCOVA — covariate-adjusted group comparison
# Variance reduction factor: 1 - rho^2 = ${(1 - p.rho_cov ** 2).toFixed(3)}
library(pwr)
pwr.anova.test(k = ${p.k},
               f = ${p.f} / sqrt(${(1 - p.rho_cov ** 2).toFixed(3)}),
               sig.level = ${p.alpha}, power = ${p.power})
# Adjusted f = ${(p.f / Math.sqrt(1 - p.rho_cov ** 2)).toFixed(4)}
# Result: n = ${r.perGroup} per group (${r.total} total)`,
    stata: `* ANCOVA — adjust ANOVA power for covariate correlation
* rho = ${p.rho_cov}, variance reduction = ${(1 - p.rho_cov ** 2).toFixed(3)}
power oneway, k(${p.k}) f(${(p.f / Math.sqrt(1 - p.rho_cov ** 2)).toFixed(4)}) ///
    alpha(${p.alpha}) power(${p.power})
* Effective f adjusted for covariate (Cohen 1988, Ch.9)`,
    sas: `/* ANCOVA power — adjust for covariate */
proc power;
  onewayanova
    stddev     = 1
    alpha      = ${p.alpha}
    power      = ${p.power}
    /* Use effective f = original f / sqrt(1 - rho^2) */
    npergroup  = .;
run;`,
  }),

  paired_proportions: (p, r) => ({
    r: `# McNemar test for paired proportions
# Discordant proportions: p01 = ${p.p01}, p10 = ${p.p10}
library(WebPower)
wp.prop(p1 = ${p.p01}, p2 = ${p.p10}, n = NULL,
        alpha = ${p.alpha}, power = ${p.power}, type = "paired")
# Result: n = ${r.total} pairs

# Alternatively (base calculation):
n <- ${r.total}  # matched pairs required
cat("McNemar n =", n, "pairs")`,
    stata: `* McNemar test — paired proportions
power pairedproportions ${p.p01} ${p.p10}, ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.273`,
    sas: `/* McNemar test — paired proportions */
proc power;
  pairedfreq
    discproportion = ${p.p10}
    odds           = ${(p.p10 / p.p01).toFixed(4)}
    alpha          = ${p.alpha}
    power          = ${p.power}
    npairs         = .;
run;
/* Reference: SAS/STAT 14.1, Ch.89 PAIREDFREQ */`,
  }),

  chi_square: (p, r) => ({
    r: `# Chi-square test (df = ${p.df}, w = ${p.w})
library(pwr)
pwr.chisq.test(w = ${p.w}, df = ${p.df},
               sig.level = ${p.alpha}, power = ${p.power})
# Result: n = ${r.total}

# WebPower:
library(WebPower)
wp.chisq(w = ${p.w}, df = ${p.df}, alpha = ${p.alpha}, power = ${p.power})`,
    stata: `* Chi-square — use proportions approach
* For contingency table (w = ${p.w}, df = ${p.df}):
power twoproportions, power(${p.power}) alpha(${p.alpha})
* Note: specify row/column proportions for exact chi-square`,
    sas: `/* Chi-square test */
/* SAS PROC POWER does not have a direct chi-square statement.
   Use the formula: n = (z_a + z_b)^2 / w^2 */
data _null_;
  w = ${p.w}; df = ${p.df};
  alpha = ${p.alpha}; power = ${p.power};
  z_a = probit(1 - alpha/2);
  z_b = probit(power);
  n = ceil((z_a + z_b)**2 / w**2);
  put "Required n = " n;
run;`,
  }),

  two_correlations: (p, r) => ({
    r: `# Comparing two independent correlations (Fisher z-test)
library(pwr)
# Group 1 (rho1 = ${p.rho1}):
pwr.r.test(r = ${p.rho1}, sig.level = ${p.alpha}, power = ${p.power})
# Group 2 (rho2 = ${p.rho2}):
pwr.r.test(r = ${p.rho2}, sig.level = ${p.alpha}, power = ${p.power})
# For two-sample comparison: n = ${r.n1} per group (${r.total} total)

# Precise two-sample test via cocor package:
library(cocor)
# n per group needed: ${r.n1}`,
    stata: `* Comparing two independent correlations
power twocorrelations ${p.rho1} ${p.rho2}, ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.337`,
    sas: `/* Two independent correlations — Fisher z-test */
data _null_;
  r1 = ${p.rho1}; r2 = ${p.rho2};
  z1 = 0.5 * log((1+r1)/(1-r1));
  z2 = 0.5 * log((1+r2)/(1-r2));
  za = probit(1 - ${p.alpha}/2);
  zb = probit(${p.power});
  n  = ceil((za + zb)**2 / (z1-z2)**2 + 3);
  put "n per group = " n;
run;`,
  }),

  simple_slope: (p, r) => ({
    r: `# Simple linear regression — slope test
# (Equivalent to correlation power for rho = ${p.rho})
library(pwr)
pwr.f2.test(u = 1, f2 = ${(p.rho ** 2 / (1 - p.rho ** 2)).toFixed(4)},
            sig.level = ${p.alpha}, power = ${p.power})
# Result: n = ${r.total}

# Or via correlation:
pwr.r.test(r = ${p.rho}, sig.level = ${p.alpha}, power = ${p.power})`,
    stata: `* Simple regression slope test
power oneslope, r2(${(p.rho ** 2).toFixed(4)}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.421`,
    sas: `/* Simple linear regression slope */
proc power;
  multreg
    numpredictors = 1
    rsquare       = ${(p.rho ** 2).toFixed(4)}
    alpha         = ${p.alpha}
    power         = ${p.power}
    ntotal        = .;
run;`,
  }),

  logistic_regression: (p, r) => ({
    r: `# Logistic regression (Hsieh, Bloch & Larsen 1998)
library(WebPower)
wp.logistic(n = NULL, p0 = ${p.p0}, p1 = NULL,
            alpha = ${p.alpha}, power = ${p.power},
            alternative = "${p.sides === 2 ? "two.sided" : "one.sided"}", OR = ${p.or})
# Result: n = ${r.total}

# Alternatively:
library(pwr)
# f2 for logistic (approximate):
f2 <- ${p.or > 0 ? ((Math.log(p.or) ** 2 * p.p0 * (1 - p.p0)) / 4).toFixed(4) : "NA"}
pwr.f2.test(u = 1, f2 = f2, sig.level = ${p.alpha}, power = ${p.power})`,
    stata: `* Logistic regression — no native Stata power logistic
* Use R: library(WebPower); wp.logistic(OR = ${p.or}, p0 = ${p.p0})
* Or SAS PROC POWER with LOGISTIC statement`,
    sas: `/* Logistic regression */
proc power;
  logistic
    alpha          = ${p.alpha}
    power          = ${p.power}
    /* responseprob and oddsratio: */
    responseprob   = ${p.p0}
    oddsratio      = ${p.or}
    ntotal         = .;
run;
/* Reference: SAS/STAT 14.1, Ch.89 LOGISTIC */`,
  }),

  non_inferiority_prop: (p, r) => ({
    r: `# Non-inferiority trial — proportions
# delta = ${p.delta}, p1 = ${p.p1}, p2 = ${p.p2}
library(TrialSize)
n.NI.Prop(alpha = ${p.alpha}, beta = ${1 - p.power},
           p1 = ${p.p1}, p2 = ${p.p2}, delta = ${p.delta})
# Result: n = ${r.n1} per group (${r.total} total)

# Manual formula:
za1 <- qnorm(1 - ${p.alpha})  # one-sided
zb  <- qnorm(${p.power})
pbar <- (${p.p1} + ${p.p2}) / 2
n <- ceiling((za1 * sqrt(2 * pbar * (1-pbar)) +
              zb * sqrt(${p.p1}*(1-${p.p1}) + ${p.p2}*(1-${p.p2})))^2 /
             (${p.delta} - (${p.p1} - ${p.p2}))^2)
cat("n per group =", n)`,
    stata: `* Non-inferiority — proportions
power twoproportions ${p.p1} ${p.p2}, ///
    margin(${p.delta}) alpha(${p.alpha}) power(${p.power}) onesided
* Reference: Stata PSS Manual (2023), PSS-2 p.232 (margin option)`,
    sas: `/* Non-inferiority proportions */
proc power;
  twosamplefreq
    refproportion = ${p.p1}
    proportion    = ${p.p2}
    margin        = ${p.delta}
    alpha         = ${p.alpha}
    power         = ${p.power}
    sides         = 1
    ntotal        = .;
run;`,
  }),

  equivalence_means: (p, r) => ({
    r: `# Equivalence (TOST) — means
# Margin: ±delta = ±${p.delta}, SD = ${p.sd}
library(TOSTER)
powerTOSTtwo(alpha = ${p.alpha},
             statistical_power = ${p.power},
             low_eqbound_d  = ${(-p.delta / p.sd).toFixed(4)},
             high_eqbound_d = ${(p.delta / p.sd).toFixed(4)})
# Result: n = ${r.n1} per group (${r.total} total)

# Manual formula (one-sided TOST):
za <- qnorm(1 - ${p.alpha})
zb <- qnorm(${p.power})
n  <- ceiling(2 * ${p.sd}^2 * (za + zb)^2 / (${p.delta} - ${Math.abs(p.mu_diff)})^2)
cat("n per group =", n)`,
    stata: `* Equivalence — means (TOST)
power twomeans 0 ${p.mu_diff || 0}, ///
    sd(${p.sd}) margin(${p.delta}) ///
    alpha(${p.alpha}) power(${p.power}) onesided
* Reference: Stata PSS Manual (2023), PSS-2 p.144 (margin option)`,
    sas: `/* Equivalence (TOST) — means */
proc power;
  twosamplemeans
    lower  = ${-p.delta}
    upper  = ${p.delta}
    stddev = ${p.sd}
    alpha  = ${p.alpha}
    power  = ${p.power}
    ntotal = .;
run;`,
  }),

  logrank_cluster: (p, r) => ({
    r: `# Cluster RCT — log-rank (survival)
# HR = ${p.hr}, ICC = ${p.icc}, cluster size = ${p.m}
library(WebPower)
wp.crt2arm(n = ${p.m}, J = NULL, icc = ${p.icc},
           power = ${p.power}, alpha = ${p.alpha})
# Then apply to survival: events needed = ${r.events}
# DEFF = 1 + (m-1)*ICC = ${(1 + (p.m - 1) * p.icc).toFixed(3)}
# Total n = ${r.total} (${r.n1} clusters per arm)`,
    stata: `* Cluster RCT — log-rank survival
power logrank, hratio(${p.hr}) eventprob(${p.p_event}) ///
    icc(${p.icc}) m(${p.m}) ///
    alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.601`,
    sas: `/* Cluster log-rank: no direct SAS procedure.
   Step 1: Calculate individual log-rank n
   Step 2: Multiply by DEFF = ${(1 + (p.m - 1) * p.icc).toFixed(3)}
   Use R WebPower: wp.crt2arm() for cluster adjustment */`,
  }),

  cox_regression: (p, r) => ({
    r: `# Cox proportional hazards (Schoenfeld 1983)
library(powerSurvEpi)
ssizeCox(power = ${p.power}, k = 1, psi = ${p.p_event},
         rho2 = ${p.r2}, alpha = ${p.alpha},
         alternative = "${p.sides === 2 ? "two.sided" : "one.sided"}")
# Required events = ${r.events}, Total n = ${r.total}

# Or using gsDesign:
library(gsDesign)
nSurv(hr = ${p.hr}, sided = ${p.sides},
      alpha = ${p.alpha}, beta = ${1 - p.power},
      ratio = 1)`,
    stata: `* Cox proportional hazards
power cox, hratio(${p.hr}) eventprob(${p.p_event}) ///
    r2(${p.r2}) alpha(${p.alpha}) power(${p.power})
* Reference: Stata PSS Manual (2023), PSS-2 p.519`,
    sas: `/* Cox proportional hazards */
proc power;
  twosamplesurvival
    hazardratio  = ${p.hr}
    alpha        = ${p.alpha}
    power        = ${p.power}
    /* specify accrualtime and followuptime */
    ntotal       = .;
run;`,
  }),

  icc_reliability: (p, r) => ({
    r: `# ICC reliability study (Donner & Eliasziw 1987)
# icc0 = ${p.icc0}, icc1 = ${p.icc1}, k = ${p.k} raters
library(ICC.Sample.Size)
pwr.icc(rho0 = ${p.icc0}, rho1 = ${p.icc1},
        k = ${p.k}, alpha = ${p.alpha}, power = ${p.power})
# Result: n = ${r.total} subjects

# Verify manually:
f0 <- (1 + (${p.k}-1)*${p.icc0}) / (1-${p.icc0})
f1 <- (1 + (${p.k}-1)*${p.icc1}) / (1-${p.icc1})
# Use Fisher z-transform of F-ratios to find n`,
    stata: `* ICC reliability — no native Stata command
* Use R: library(ICC.Sample.Size)
* pwr.icc(rho0=${p.icc0}, rho1=${p.icc1}, k=${p.k},
*         alpha=${p.alpha}, power=${p.power})
* Result: n = ${r.total} subjects × ${p.k} raters`,
    sas: `/* ICC reliability — no SAS PROC POWER statement.
   Use R package ICC.Sample.Size:
   pwr.icc(rho0=${p.icc0}, rho1=${p.icc1},
           k=${p.k}, alpha=${p.alpha}, power=${p.power})
   Result: n = ${r.total} subjects */`,
  }),

  simon_two_stage: (p, r) => ({
    r: `# Simon two-stage design (Phase II oncology)
library(clinfun)
ph2simon(pu = ${p.p0}, pa = ${p.p1},
         ep1 = ${p.alpha}, ep2 = ${p.beta})
# Optimal design:
# Stage 1: n = ${r.n1}, stop if <= ${r.note?.match(/≤(\d+) responses/)?.[1] || "?"} responses
# Total:   n = ${r.total}
# E[n|H0] = ${r.note?.match(/E\[n\|H₀\]=([\d.]+)/)?.[1] || "?"}`,
    stata: `* Simon two-stage — no native Stata command
* Use R: library(clinfun)
* ph2simon(pu = ${p.p0}, pa = ${p.p1}, ep1 = ${p.alpha}, ep2 = ${p.beta})`,
    sas: `/* Simon two-stage — no SAS PROC POWER statement.
   Use R package clinfun:
   ph2simon(pu=${p.p0}, pa=${p.p1},
            ep1=${p.alpha}, ep2=${p.beta}) */`,
  }),

  survey_mean: (p, r) => ({
    r: `# Survey — estimating a mean
# SD = ${p.sd}, margin of error = ±${p.d}
z    <- qnorm(1 - ${p.alpha}/2)
n_raw  <- ceiling(z^2 * ${p.sd}^2 / ${p.d}^2)
n_deff <- ceiling(n_raw * ${p.deff})
n_final <- ceiling(n_deff / ${p.response_rate / 100})
cat("Required n:", n_final)   # = ${r.total}

# PracTools package (Valliant et al. 2013):
library(PracTools)
nCont(CV0 = ${p.d}/${p.sd}, S2 = ${p.sd ** 2}, N = Inf)`,
    stata: `* Survey — estimating a mean
ciwidth onemean, width(${p.d * 2}) sd(${p.sd}) ///
    level(${(1 - p.alpha) * 100})
* Then multiply by DEFF = ${p.deff}
* Reference: Stata PSS Manual (2023), PSS-3 p.702`,
    sas: `/* Survey — mean estimation */
data _null_;
  z    = probit(${1 - p.alpha / 2});
  sd   = ${p.sd}; d = ${p.d};
  deff = ${p.deff}; rr = ${p.response_rate / 100};
  n    = ceil(ceil(z**2 * sd**2 / d**2) * deff / rr);
  put "Required n = " n;
run;`,
  }),

  mediation_sobel: (p, r) => ({
    r: `# Simple mediation — Sobel test
library(WebPower)
wp.mediation(n = NULL,
             a = ${p.a}, b = ${p.b},
             varx = ${p.sd_a ** 2}, varm = ${p.sd_b ** 2},
             alpha = ${p.alpha}, power = ${p.power})
# Result: n = ${r.total}

# Fritz & MacKinnon (2007) tables also available for
# standardised a*b product: ${(p.a * p.b).toFixed(4)}`,
    stata: `* Mediation (Sobel) — no native Stata command
* Use R: library(WebPower)
* wp.mediation(a=${p.a}, b=${p.b}, varx=${p.sd_a ** 2},
*              varm=${p.sd_b ** 2}, alpha=${p.alpha}, power=${p.power})
* Result: n = ${r.total}`,
    sas: `/* Mediation (Sobel) — no SAS PROC POWER statement.
   Use R WebPower:
   wp.mediation(a=${p.a}, b=${p.b},
                varx=${p.sd_a ** 2}, varm=${p.sd_b ** 2},
                alpha=${p.alpha}, power=${p.power})
   Result: n = ${r.total} */`,
  }),
};

// ─────────────────────────────────────────────
// MANUSCRIPT / SUMMARY GENERATOR — all 26 designs
// ─────────────────────────────────────────────
function generateParagraph(formulaKey, params, result) {
  if (!result || result.error || !result.total) return null;
  const cites = CITATIONS[formulaKey] || [];
  const c0 = cites[0] ? `${cites[0].ref} (${cites[0].detail})` : "";
  const c1 = cites[1] ? `${cites[1].ref} (${cites[1].detail})` : c0;
  const att = params.attrition > 0
    ? ` To account for an anticipated attrition rate of ${params.attrition}%, the sample size was inflated accordingly, yielding a final target of ${result.total}.`
    : "";
  const sides = params.sides === 1 ? "one-sided" : "two-sided";
  const es = (v, d = 3) => typeof v === "number" ? v.toFixed(d) : v;

  const paragraphs = {
    two_means: (() => {
      const effectDesc = params.effectMode === "d"
        ? `Cohen's d = ${es(result.effectSize)}`
        : `a mean of ${params.mu1} in group 1 and ${params.mu2} in group 2 with a pooled SD of ${params.sd} (Cohen's d = ${es(result.effectSize)})`;
      return `Sample size was calculated for a two-sample independent t-test comparing means between two parallel groups. Assuming ${effectDesc}, a ${sides} significance level of ${params.alpha}, and ${params.power * 100}% power, a minimum of ${result.n1} participants per group was required (${c1}).${att} The total sample size is ${result.total}.`;
    })(),

    paired_means:
      `Sample size was determined for a paired-samples t-test. Assuming a mean difference of ${params.mu_diff} with a standard deviation of differences of ${params.sd_diff} (Cohen's d = ${es(result.effectSize)}), a ${sides} alpha of ${params.alpha}, and ${params.power * 100}% power, a minimum of ${result.total} pairs was required (${c0}).${att}`,

    one_mean:
      `Sample size was calculated for a one-sample t-test comparing the sample mean against a null value of ${params.mu0}. Assuming a true mean of ${params.mu1}, standard deviation of ${params.sd} (Cohen's d = ${es(result.effectSize)}), a ${sides} alpha of ${params.alpha}, and ${params.power * 100}% power, a minimum of ${result.total} participants was required (${c0}).${att}`,

    two_proportions:
      `Sample size was calculated using the formula for comparing two independent proportions (${c1}). Assuming a proportion of ${(params.p1 * 100).toFixed(1)}% in the control group and ${(params.p2 * 100).toFixed(1)}% in the intervention group (Cohen's h = ${es(result.effectSize)}), with a ${sides} significance level of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.n1} participants per group was required.${att} The total sample size is ${result.total}.`,

    one_proportion:
      `Sample size was calculated for a one-sample proportion test, comparing an expected proportion of ${(params.p1 * 100).toFixed(1)}% against a null hypothesis of ${(params.p0 * 100).toFixed(1)}% (Cohen's h = ${es(result.effectSize)}). With a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.total} participants was required (${c0}).${att}`,

    one_way_anova:
      `Sample size was determined for a one-way analysis of variance (ANOVA) with ${params.k} independent groups. Assuming a medium effect size of Cohen's f = ${params.f} (${c0}), a significance level of ${params.alpha}, and ${params.power * 100}% power, a minimum of ${result.perGroup} participants per group was required, giving a total sample size of ${result.total}.${att}`,

    correlation:
      `Sample size was calculated for testing the significance of a Pearson correlation coefficient using Fisher's z-transformation. Assuming a population correlation of ρ = ${params.rho} (${c0}), with a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.total} participants was required.${att}`,

    regression_r2:
      `Sample size was determined for a multiple regression analysis with ${params.u} predictors, using Cohen's f² = ${params.f2} as the effect size measure (${c0}). With a significance level of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.total} participants was required, yielding ${result.total - params.u - 1} residual degrees of freedom.${att}`,

    logrank:
      `Sample size was calculated for a log-rank test comparing survival curves between two arms, using the formula described by ${cites[0]?.ref} (${cites[0]?.detail}). Assuming a hazard ratio of ${params.hr} and an overall event probability of ${(params.p_event * 100).toFixed(0)}%, with a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.events} events were required. This corresponds to a total sample size of ${result.total} (${result.n1} per arm).${att}`,

    survey_proportion:
      `The required sample size for this cross-sectional survey was estimated assuming a population proportion of ${(params.p * 100).toFixed(0)}%, a margin of error of ±${(params.d * 100).toFixed(0)}%, and a 95% confidence level (${c0}). A design effect (DEFF) of ${params.deff} was applied to account for the sampling design. After adjusting for an expected response rate of ${params.response_rate}%, a final sample of ${result.total} participants was required.`,

    sem_rmsea:
      `Statistical power for the structural equation model was evaluated using the RMSEA-based method described by MacCallum, Browne, and Sugawara (1996, Psychological Methods 1(2):130–149), as implemented in the semPower package (${c1}). The model had ${params.df} degrees of freedom. Power was calculated to reject a null hypothesis of RMSEA ≤ ${params.rmsea0} when the true population RMSEA was ${params.rmsea1}, with a ${sides} alpha of ${params.alpha} and desired power of ${params.power * 100}%. A minimum sample size of N = ${result.total} was required.`,

    repeated_measures:
      `Sample size was determined for a one-way repeated measures analysis of variance with ${params.k} measurement time points. Assuming a medium effect size of Cohen's f = ${params.f} and an average within-subject correlation of ρ = ${params.rho} between measurements (${c0}), with a significance level of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.total} subjects was required.${att} Each subject was measured ${params.k} times, yielding ${result.total * params.k} total observations.`,

    two_way_anova:
      `Sample size was calculated for a ${params.a} × ${params.b} factorial analysis of variance. Effect sizes (Cohen's f) were set at f = ${params.fa} for factor A, f = ${params.fb} for factor B, and f = ${params.fab} for the A × B interaction (${c0}). With a significance level of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.perGroup} participants per cell was required, giving a total sample size of ${result.total} across ${params.a * params.b} cells.${att}`,

    ancova:
      `Sample size was calculated for a one-way analysis of covariance (ANCOVA) with ${params.k} groups and one baseline covariate. Assuming Cohen's f = ${params.f} for the adjusted group effect and a correlation of ρ = ${params.rho_cov} between the outcome and the covariate (variance reduction factor = ${(1 - params.rho_cov ** 2).toFixed(3)}), with alpha = ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.perGroup} participants per group (${result.total} total) was required (${c0}).${att}`,

    paired_proportions:
      `Sample size was determined for a McNemar test comparing paired binary proportions. With discordant pair proportions of ${(params.p01 * 100).toFixed(1)}% (negative/positive) and ${(params.p10 * 100).toFixed(1)}% (positive/negative), a ${sides} alpha of ${params.alpha}, and ${params.power * 100}% power, a minimum of ${result.total} matched pairs was required (${c0}).${att}`,

    chi_square:
      `Sample size was determined for a chi-square test with ${params.df} degree(s) of freedom. Using Cohen's w = ${params.w} as the effect size (small = 0.10, medium = 0.30, large = 0.50; ${c0}), with a significance level of ${params.alpha} and ${params.power * 100}% power, a minimum total sample of ${result.total} was required.${att}`,

    two_correlations:
      `Sample size was calculated to detect a difference between two independent Pearson correlation coefficients using Fisher's z-transformation. Assuming ρ₁ = ${params.rho1} in group 1 and ρ₂ = ${params.rho2} in group 2 (${c0}), with a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.n1} participants per group (${result.total} total) was required.${att}`,

    simple_slope:
      `Sample size was calculated for a slope test in simple linear regression. Assuming a population correlation of ρ = ${params.rho} between the predictor and outcome (equivalent to the standardised regression coefficient; ${c0}), with a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.total} participants was required.${att}`,

    logistic_regression:
      `Sample size was estimated for a logistic regression analysis using the approximation method of Hsieh, Bloch & Larsen (1998, Statistics in Medicine 17:1503–1506), as implemented in WebPower (${c0}). Assuming an odds ratio of ${params.or} for the primary predictor, a baseline event probability of ${(params.p0 * 100).toFixed(1)}%, and R² = ${params.r2} for other covariates in the model, with a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.total} participants was required.${att}`,

    non_inferiority_prop:
      `Sample size was calculated for a non-inferiority trial comparing two proportions. The non-inferiority margin was set at δ = ${params.delta}. Assuming a reference proportion of ${(params.p1 * 100).toFixed(1)}% and an expected test proportion of ${(params.p2 * 100).toFixed(1)}%, with a one-sided alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.n1} participants per group (${result.total} total) was required (${c0}).${att}`,

    equivalence_means:
      `Sample size was calculated for an equivalence trial using the two one-sided tests (TOST) procedure. The equivalence margin was set at ±δ = ${params.delta}. Assuming a standard deviation of ${params.sd} and an expected true difference of ${params.mu_diff}, with a one-sided alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.n1} participants per group (${result.total} total) was required (${c0}).${att}`,

    logrank_cluster:
      `Sample size was calculated for a cluster randomised controlled trial with a time-to-event primary outcome. Using the clustered log-rank approach (${c0}), and assuming a hazard ratio of ${params.hr}, overall event probability of ${(params.p_event * 100).toFixed(0)}%, intra-cluster correlation of ${params.icc}, and average cluster size of ${params.m} (design effect = ${result.note?.match(/DEFF = ([\d.]+)/)?.[1] || "—"}), with a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.n1} clusters per arm was required. The total sample size is ${result.total}.${att}`,

    cox_regression:
      `Sample size was estimated for a Cox proportional hazards regression analysis using the formula of Schoenfeld (1983, Biometrics 39:499–503). Assuming a hazard ratio of ${params.hr} for the primary predictor, an overall event probability of ${(params.p_event * 100).toFixed(0)}%, and R² = ${params.r2} for other covariates, with a ${sides} alpha of ${params.alpha} and ${params.power * 100}% power, a minimum of ${result.events} events (total n = ${result.total}) was required (${c0}).${att}`,

    icc_reliability:
      `Sample size was calculated for a reliability study using the intra-class correlation coefficient (ICC). Setting a minimum acceptable ICC of ${params.icc0} (H₀) and an expected true ICC of ${params.icc1} (H₁), with ${params.k} repeated measurements per subject, a ${sides} alpha of ${params.alpha}, and ${params.power * 100}% power, a minimum of ${result.total} subjects was required (${c0}).${att}`,

    simon_two_stage:
      result.total
        ? `A Simon two-stage design was used for this single-arm phase II trial (${c0}; ${cites[1]?.ref}, ${cites[1]?.detail}). With an unacceptable response rate of p₀ = ${(params.p0 * 100).toFixed(0)}% and a target response rate of p₁ = ${(params.p1 * 100).toFixed(0)}%, using a one-sided alpha of ${params.alpha} and beta of ${params.beta}, the optimal design requires ${result.n1} patients in Stage 1: the drug is abandoned if ${result.note?.match(/≤(\d+) responses/)?.[1] || "—"} or fewer responses are observed. If Stage 1 proceeds, a total of ${result.total} patients are enrolled. The expected sample size under H₀ is ${result.note?.match(/E\[n\|H₀\]=([\d.]+)/)?.[1] || "—"}.`
        : "No feasible Simon two-stage design found with the specified parameters. Consider relaxing alpha or beta.",

    survey_mean:
      `The required sample size for estimating the population mean was calculated assuming a population standard deviation of ${params.sd} and a desired margin of error of ±${params.d} at the 95% confidence level (${c0}). A design effect of ${params.deff} was applied to account for the sampling design. After adjusting for an expected response rate of ${params.response_rate}%, a final sample of ${result.total} was required.`,

    mediation_sobel:
      `Sample size was estimated for a simple mediation analysis using the Sobel test (${c0}). Assuming path coefficients of a = ${params.a} (X → M) and b = ${params.b} (M → Y), with corresponding standard deviations of ${params.sd_a} and ${params.sd_b}, an alpha of ${params.alpha}, and ${params.power * 100}% power, a minimum of ${result.total} participants was required.${att}`,
  };

  return paragraphs[formulaKey] || null;
}

// ─────────────────────────────────────────────
// LANDING PAGE — first page / about screen
// ─────────────────────────────────────────────
// ─────────────────────────────────────────────
// OVERVIEW DATA — used by landing page
// ─────────────────────────────────────────────
const OVERVIEW_DATA = {
  designs: [
    { group: "Means",          color: "#3b82f6", icon: "〰", count: 3,  items: ["Two independent means", "Paired / matched means", "One-sample mean"] },
    { group: "Proportions",    color: "#8b5cf6", icon: "⊡", count: 5,  items: ["Two independent proportions", "One-sample proportion", "Paired proportions (McNemar)", "Chi-square test", "Non-inferiority — proportions"] },
    { group: "ANOVA",          color: "#06b6d4", icon: "⊞", count: 4,  items: ["One-way ANOVA", "Two-way ANOVA", "Repeated measures ANOVA", "ANCOVA (covariate-adjusted)"] },
    { group: "Regression",     color: "#10b981", icon: "↗", count: 5,  items: ["Pearson correlation", "Two correlations (compare)", "Simple linear regression", "Multiple regression (R²)", "Logistic regression"] },
    { group: "Survival",       color: "#f59e0b", icon: "⏳", count: 3,  items: ["Log-rank test", "Cluster RCT — log-rank", "Cox proportional hazards"] },
    { group: "Clinical Trials",color: "#dc2626", icon: "⚖", count: 1,  items: ["Equivalence (TOST) — means"] },
    { group: "Diagnostic",     color: "#7c3aed", icon: "🔬", count: 1,  items: ["ICC — reliability / repeatability"] },
    { group: "Phase II",       color: "#0891b2", icon: "🧪", count: 1,  items: ["Simon two-stage (Phase II oncology)"] },
    { group: "Survey",         color: "#ef4444", icon: "📋", count: 2,  items: ["Survey — prevalence / proportion", "Survey — estimating a mean"] },
    { group: "SEM",            color: "#ec4899", icon: "⬡", count: 2,  items: ["SEM — RMSEA-based power", "Simple mediation (Sobel test)"] },
  ],
  references: [
    { short: "Cohen (1988)",            full: "Statistical Power Analysis for the Behavioral Sciences, 2nd ed.",          covers: "t-tests, ANOVA, correlation, regression, chi-square" },
    { short: "Aberson (2019)",          full: "Applied Power Analysis for the Behavioral Sciences, 2nd ed.",              covers: "R-based examples for all major designs" },
    { short: "Zhang & Yuan (2018)",     full: "Practical Statistical Power Analysis Using WebPower and R",                covers: "Cluster RCT, mediation, SEM, logistic, Poisson" },
    { short: "Jak et al. (2020)",       full: "Analytical power calculations for SEM: A tutorial and Shiny app",          covers: "SEM χ² and RMSEA power — power4SEM" },
    { short: "Machin et al. (2009)",    full: "Sample Size Tables for Clinical Studies, 3rd ed.",                         covers: "Clinical tables, diagnostics, phase II, ICC, observer agreement" },
    { short: "Chow et al. (2008)",      full: "Sample Size Calculations in Clinical Research, 2nd ed.",                   covers: "FDA/regulatory, equivalence, crossover, survival" },
    { short: "Aday & Cornelius (2006)", full: "Designing and Conducting Health Surveys, 3rd ed.",                         covers: "Health survey sample size, stratification" },
    { short: "Valliant et al. (2013)",  full: "Practical Tools for Designing and Weighting Survey Samples",               covers: "Complex survey design, DEFF, PracTools R package" },
    { short: "Fleiss et al. (2003)",    full: "Statistical Methods for Rates and Proportions, 3rd ed.",                   covers: "Proportions, OR, case-control — Eq. 4.19" },
    { short: "Moshagen & Bader (2024)", full: "semPower: General Power Analysis for SEM. Behav Res Methods 56:2901–2922", covers: "Most current R package for SEM power" },
    { short: "Stata PSS Manual (2023)", full: "Stata 18 Power, Precision, and Sample-Size Reference Manual",              covers: "All power & ciwidth commands with exact page numbers" },
    { short: "SAS PROC POWER (2015)",   full: "SAS/STAT 14.1 User's Guide, Chapter 89: The POWER Procedure",             covers: "SAS syntax for all major designs incl. logistic" },
  ],
  outputs: [
    { icon: "📊", label: "Interpretation",  desc: "Plain-language explanation of the result, effect size benchmarks, and what each assumption means" },
    { icon: "📝", label: "Summary",         desc: "Ready-to-paste manuscript paragraph with embedded citations — copy, share via WhatsApp, or export to PDF" },
    { icon: "💻", label: "Code",            desc: "Verified R, Stata, and SAS commands to independently reproduce the calculation" },
    { icon: "📚", label: "Citations",       desc: "Authoritative textbook references with chapter and page numbers — primary and backup citations per design" },
  ],
};


function LandingPage({ onEnter }) {
  const [refOpen, setRefOpen] = useState(false);
  const totalDesigns = OVERVIEW_DATA.designs.reduce((s, d) => s + d.count, 0);

  return (
    <div style={{
      fontFamily: "'Inter', system-ui, sans-serif",
      background: "#ffffff",
      minHeight: "100vh", overflowY: "auto", color: "#1e3a5f",
    }}>
      {/* Hero */}
      <div style={{ textAlign: "center", padding: "52px 24px 36px" }}>
        <div style={{
          width: 64, height: 64, borderRadius: 18,
          background: "#eef2ff", border: "1px solid #c7d2fe",
          display: "flex", alignItems: "center", justifyContent: "center",
          fontSize: 30, fontWeight: 700, margin: "0 auto 20px",
          boxShadow: "0 8px 32px rgba(99,102,241,0.35)",
        }}>Σ</div>

        <h1 style={{ fontSize: 32, fontWeight: 800, letterSpacing: "-0.03em", margin: "0 0 10px", lineHeight: 1.15 }}>
          StatCalc <span style={{ color: "#818cf8" }}>MD</span>
        </h1>
        <p style={{ fontSize: 16, color: "#1e3a5f", margin: "0 0 28px", fontWeight: 500 }}>
          Sample size calculator for medical researchers & clinicians
        </p>

        {/* Stat pills */}
        <div style={{ display: "flex", flexWrap: "wrap", gap: 8, justifyContent: "center", marginBottom: 36 }}>
          {[
            { val: totalDesigns, label: "study designs" },
            { val: "12", label: "reference books" },
            { val: "R · Stata · SAS", label: "code output" },
            { val: "26", label: "manuscript paragraphs" },
            { val: "WhatsApp + PDF", label: "share & export" },
          ].map(({ val, label }) => (
            <div key={label} style={{
              display: "flex", alignItems: "center", gap: 6,
              background: "#eef2ff", border: "1px solid #c7d2fe",
              borderRadius: 24, padding: "5px 14px",
            }}>
              <span style={{ fontSize: 13, fontWeight: 700, color: "#1e3a5f" }}>{val}</span>
              <span style={{ fontSize: 11, color: "#4f46e5" }}>{label}</span>
            </div>
          ))}
        </div>

        {/* Description */}
        <p style={{ fontSize: 14, color: "#1e3a5f", lineHeight: 1.75, maxWidth: 620, margin: "0 auto 40px" }}>
          Every calculation runs locally in your browser — no server, no data sent anywhere. For each design
          you get a <strong style={{ color: "#1e3a5f" }}>plain-language interpretation</strong>, a{" "}
          <strong style={{ color: "#1e3a5f" }}>ready-to-paste manuscript paragraph</strong> with embedded textbook
          citations, <strong style={{ color: "#1e3a5f" }}>verified R / Stata / SAS code</strong> to reproduce the
          result independently, and the exact book chapter and page number so a reviewer can check the formula themselves.
        </p>

        {/* Enter button */}
        <button
          onClick={onEnter}
          style={{
            padding: "14px 48px", borderRadius: 14,
            background: "linear-gradient(135deg, #4f46e5, #6366f1)",
            color: "#fff", border: "none", fontSize: 16, fontWeight: 700,
            cursor: "pointer", letterSpacing: "0.01em",
            boxShadow: "0 8px 32px rgba(99,102,241,0.50)",
            transition: "transform 0.15s, box-shadow 0.15s",
          }}
          onMouseEnter={e => { e.target.style.transform = "translateY(-2px)"; e.target.style.boxShadow = "0 12px 40px rgba(99,102,241,0.60)"; }}
          onMouseLeave={e => { e.target.style.transform = "translateY(0)"; e.target.style.boxShadow = "0 8px 32px rgba(99,102,241,0.50)"; }}
        >
          Enter Calculator →
        </button>
        <p style={{ fontSize: 11, color: "#64748b", marginTop: 10 }}>No sign-up required · Runs offline</p>
      </div>

      {/* What each result gives you */}
      <div style={{ maxWidth: 860, margin: "0 auto", padding: "0 20px 32px" }}>
        <div style={{ fontSize: 10, fontWeight: 700, color: "#4f46e5", letterSpacing: "0.09em", marginBottom: 12, textAlign: "center" }}>
          WHAT EACH RESULT GIVES YOU
        </div>
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(190px, 1fr))", gap: 10 }}>
          {OVERVIEW_DATA.outputs.map(o => (
            <div key={o.label} style={{
              background: "#f8fafc", border: "1px solid #e2e8f0",
              borderRadius: 12, padding: "14px 16px",
            }}>
              <div style={{ fontSize: 22, marginBottom: 6 }}>{o.icon}</div>
              <div style={{ fontSize: 13, fontWeight: 700, color: "#1e3a5f", marginBottom: 5 }}>{o.label}</div>
              <div style={{ fontSize: 11.5, color: "#64748b", lineHeight: 1.55 }}>{o.desc}</div>
            </div>
          ))}
        </div>
      </div>

      {/* Supported designs */}
      <div style={{ maxWidth: 860, margin: "0 auto", padding: "0 20px 32px" }}>
        <div style={{ fontSize: 10, fontWeight: 700, color: "#4f46e5", letterSpacing: "0.09em", marginBottom: 12, textAlign: "center" }}>
          {totalDesigns} SUPPORTED STUDY DESIGNS ACROSS 10 CATEGORIES
        </div>
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(195px, 1fr))", gap: 9 }}>
          {OVERVIEW_DATA.designs.map(d => (
            <div key={d.group} style={{
              background: "#f8fafc", border: `1px solid ${d.color}35`,
              borderRadius: 11, padding: "11px 14px",
            }}>
              <div style={{ display: "flex", alignItems: "center", gap: 7, marginBottom: 7 }}>
                <span style={{ fontSize: 15 }}>{d.icon}</span>
                <span style={{ fontSize: 12, fontWeight: 700, color: d.color }}>{d.group}</span>
                <span style={{
                  marginLeft: "auto", fontSize: 10, fontWeight: 700,
                  background: d.color + "25", color: d.color,
                  borderRadius: 10, padding: "1px 7px",
                }}>{d.count}</span>
              </div>
              {d.items.map(item => (
                <div key={item} style={{ fontSize: 10.5, color: "#475569", lineHeight: 1.65, paddingLeft: 2 }}>
                  · {item}
                </div>
              ))}
            </div>
          ))}
        </div>
      </div>

      {/* References */}
      <div style={{ maxWidth: 860, margin: "0 auto", padding: "0 20px 48px" }}>
        <button
          onClick={() => setRefOpen(r => !r)}
          style={{
            display: "flex", alignItems: "center", gap: 7, width: "100%", justifyContent: "center",
            background: "#f1f5f9", border: "1px solid #e2e8f0",
            borderRadius: 10, padding: "10px 20px", color: "#4f46e5",
            fontSize: 11, fontWeight: 700, letterSpacing: "0.07em", cursor: "pointer", marginBottom: 12,
          }}
        >
          <span>{refOpen ? "▼" : "▶"}</span>
          <span>12 AUTHORITATIVE REFERENCE SOURCES {refOpen ? "— click to collapse" : "— click to expand"}</span>
        </button>

        {refOpen && (
          <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(380px, 1fr))", gap: 7 }}>
            {OVERVIEW_DATA.references.map((ref, i) => (
              <div key={i} style={{
                background: "#f8fafc", border: "1px solid #e2e8f0",
                borderRadius: 9, padding: "10px 13px", display: "flex", gap: 10, alignItems: "flex-start",
              }}>
                <span style={{
                  minWidth: 20, height: 20, borderRadius: "50%",
                  background: "#4338ca", color: "#fff",
                  fontSize: 9.5, fontWeight: 700, display: "flex",
                  alignItems: "center", justifyContent: "center", flexShrink: 0, marginTop: 2,
                }}>{i + 1}</span>
                <div>
                  <div style={{ fontSize: 12, fontWeight: 700, color: "#1e3a5f" }}>{ref.short}</div>
                  <div style={{ fontSize: 10.5, color: "#4f46e5", lineHeight: 1.45, marginTop: 2 }}>{ref.full}</div>
                  <div style={{ fontSize: 10, color: "#475569", marginTop: 3 }}>Covers: {ref.covers}</div>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Bottom enter CTA */}
      <div style={{ textAlign: "center", padding: "0 20px 52px" }}>
        <button
          onClick={onEnter}
          style={{
            padding: "12px 40px", borderRadius: 12,
            background: "#eef2ff", border: "1px solid #c7d2fe",
            color: "#1e3a5f", fontSize: 14, fontWeight: 700, cursor: "pointer",
          }}
        >
          Enter Calculator →
        </button>
      </div>
    </div>
  );
}


function SummaryTab({ manuscript, formula, params, result }) {
  const [copied, setCopied] = useState(false);

  const shareWhatsApp = () => {
    if (!manuscript) return;
    const text = encodeURIComponent(
      `📊 Sample Size Calculation — ${formula.label}\n\n` +
      manuscript +
      `\n\n🔢 Result: N = ${result.total}` +
      (result.n2 ? ` (${result.n1} per group)` : "") +
      `\n⚙ α = ${params.alpha} · Power = ${params.power * 100}%` +
      `\n\n[Calculated with StatCalc MD]`
    );
    window.open(`https://wa.me/?text=${text}`, "_blank");
  };

  const exportPDF = () => {
    if (!manuscript) return;
    const content = `
      <html><head><style>
        body { font-family: Georgia, serif; max-width: 700px; margin: 40px auto; color: #1e293b; line-height: 1.8; }
        h1 { font-size: 20px; color: #312e81; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; }
        h2 { font-size: 14px; color: #64748b; font-weight: 600; margin-top: 24px; }
        .paragraph { font-size: 13px; background: #f8fafc; border-left: 4px solid #6366f1; padding: 14px 18px; border-radius: 6px; font-style: italic; }
        .grid { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 12px; margin: 16px 0; }
        .card { background: #f1f5f9; border-radius: 8px; padding: 10px 14px; }
        .card-label { font-size: 10px; text-transform: uppercase; color: #94a3b8; font-weight: 700; }
        .card-val { font-size: 22px; font-weight: 700; color: #4f46e5; }
        .warning { font-size: 11px; color: #94a3b8; margin-top: 20px; border-top: 1px solid #e2e8f0; padding-top: 12px; }
        .meta { font-size: 12px; color: #64748b; }
      </style></head><body>
        <h1>Sample Size Calculation Report</h1>
        <p class="meta"><strong>Design:</strong> ${formula.label} &nbsp;|&nbsp; <strong>Date:</strong> ${new Date().toLocaleDateString()}</p>

        <h2>KEY RESULTS</h2>
        <div class="grid">
          <div class="card"><div class="card-label">Total N</div><div class="card-val">${result.total?.toLocaleString()}</div></div>
          ${result.n2 != null ? `<div class="card"><div class="card-label">Per Group</div><div class="card-val">${result.n1}</div></div>` : ""}
          ${result.events ? `<div class="card"><div class="card-label">Events</div><div class="card-val">${result.events}</div></div>` : ""}
          ${result.effectSize != null ? `<div class="card"><div class="card-label">${result.effectLabel}</div><div class="card-val">${typeof result.effectSize === "number" ? result.effectSize.toFixed(3) : result.effectSize}</div></div>` : ""}
          <div class="card"><div class="card-label">Alpha (α)</div><div class="card-val">${params.alpha}</div></div>
          <div class="card"><div class="card-label">Power</div><div class="card-val">${params.power * 100}%</div></div>
        </div>

        <h2>MANUSCRIPT / PROTOCOL PARAGRAPH</h2>
        <div class="paragraph">${manuscript}</div>

        ${result.note ? `<p class="meta" style="margin-top:12px">ℹ ${result.note}</p>` : ""}

        <p class="warning">⚠ Review and adapt this paragraph to your specific study context before submission. All parameter values must reflect your actual study design and should be justified with reference to pilot data, literature, or clinical judgement.</p>
        <p class="warning" style="margin-top:8px">Generated with StatCalc MD · ${new Date().toLocaleString()}</p>
      </body></html>
    `;
    const w = window.open("", "_blank");
    w.document.write(content);
    w.document.close();
    w.focus();
    setTimeout(() => { w.print(); }, 500);
  };

  if (!manuscript) {
    return (
      <div style={{ background: "#fff", borderRadius: 12, padding: 24, border: "1px solid #e2e8f0", textAlign: "center", color: "#94a3b8" }}>
        <div style={{ fontSize: 28, marginBottom: 8 }}>📝</div>
        <div style={{ fontSize: 13, fontWeight: 500 }}>No summary available for this design yet.</div>
        <div style={{ fontSize: 11.5, marginTop: 4 }}>Calculate first, then check back here.</div>
      </div>
    );
  }

  return (
    <div>
      {/* Manuscript box */}
      <div style={{ background: "#fff", borderRadius: 12, padding: 16, border: "1px solid #e2e8f0", marginBottom: 10 }}>
        <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", marginBottom: 12 }}>
          <div style={{ fontSize: 13, fontWeight: 600, color: "#1e293b" }}>Manuscript / Protocol paragraph</div>
          <button
            onClick={() => { navigator.clipboard.writeText(manuscript); setCopied(true); setTimeout(() => setCopied(false), 2000); }}
            style={{ fontSize: 11, padding: "4px 10px", borderRadius: 6, border: "1px solid #e2e8f0", background: "#f8fafc", cursor: "pointer", color: "#4f46e5", fontWeight: 600 }}
          >
            {copied ? "✓ Copied" : "Copy text"}
          </button>
        </div>

        <div style={{
          fontSize: 13, lineHeight: 1.85, color: "#334155",
          padding: 14, background: "#f8fafc", borderRadius: 8,
          borderLeft: "3px solid #6366f1", fontStyle: "italic", marginBottom: 14,
        }}>
          {manuscript}
        </div>

        {/* Share / Export buttons */}
        <div style={{ display: "flex", gap: 8, flexWrap: "wrap" }}>
          <button
            onClick={shareWhatsApp}
            style={{
              display: "flex", alignItems: "center", gap: 7,
              padding: "9px 16px", borderRadius: 9, border: "none",
              background: "#25D366", color: "#fff",
              fontSize: 12.5, fontWeight: 700, cursor: "pointer",
              boxShadow: "0 2px 8px rgba(37,211,102,0.30)",
            }}
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
              <path d="M17.472 14.382c-.297-.149-1.758-.867-2.03-.967-.273-.099-.471-.148-.67.15-.197.297-.767.966-.94 1.164-.173.199-.347.223-.644.075-.297-.15-1.255-.463-2.39-1.475-.883-.788-1.48-1.761-1.653-2.059-.173-.297-.018-.458.13-.606.134-.133.298-.347.446-.52.149-.174.198-.298.298-.497.099-.198.05-.371-.025-.52-.075-.149-.669-1.612-.916-2.207-.242-.579-.487-.5-.669-.51-.173-.008-.371-.01-.57-.01-.198 0-.52.074-.792.372-.272.297-1.04 1.016-1.04 2.479 0 1.462 1.065 2.875 1.213 3.074.149.198 2.096 3.2 5.077 4.487.709.306 1.262.489 1.694.625.712.227 1.36.195 1.871.118.571-.085 1.758-.719 2.006-1.413.248-.694.248-1.289.173-1.413-.074-.124-.272-.198-.57-.347m-5.421 7.403h-.004a9.87 9.87 0 01-5.031-1.378l-.361-.214-3.741.982.998-3.648-.235-.374a9.86 9.86 0 01-1.51-5.26c.001-5.45 4.436-9.884 9.888-9.884 2.64 0 5.122 1.03 6.988 2.898a9.825 9.825 0 012.893 6.994c-.003 5.45-4.437 9.884-9.885 9.884m8.413-18.297A11.815 11.815 0 0012.05 0C5.495 0 .16 5.335.157 11.892c0 2.096.547 4.142 1.588 5.945L.057 24l6.305-1.654a11.882 11.882 0 005.683 1.448h.005c6.554 0 11.89-5.335 11.893-11.893a11.821 11.821 0 00-3.48-8.413Z"/>
            </svg>
            Share via WhatsApp
          </button>

          <button
            onClick={exportPDF}
            style={{
              display: "flex", alignItems: "center", gap: 7,
              padding: "9px 16px", borderRadius: 9, border: "1.5px solid #e2e8f0",
              background: "#fff", color: "#374151",
              fontSize: 12.5, fontWeight: 700, cursor: "pointer",
            }}
          >
            <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
              <polyline points="14 2 14 8 20 8"/>
              <line x1="16" y1="13" x2="8" y2="13"/>
              <line x1="16" y1="17" x2="8" y2="17"/>
              <polyline points="10 9 9 9 8 9"/>
            </svg>
            Export / Print PDF
          </button>
        </div>
      </div>

      <div style={{ fontSize: 11, color: "#94a3b8", padding: "0 3px" }}>
        ⚠ Review before submission. Adapt to your study context. All parameters must match your actual design.
      </div>
    </div>
  );
}


const GROUPS = ["Means", "Proportions", "ANOVA", "Regression", "Survival", "Clinical Trials", "Diagnostic", "Phase II", "Survey", "SEM"];
const TAG_COLORS = {
  Means: "#3b82f6", Proportions: "#8b5cf6", ANOVA: "#06b6d4",
  Regression: "#10b981", Survival: "#f59e0b", Survey: "#ef4444", SEM: "#ec4899",
  "Clinical Trials": "#dc2626", Diagnostic: "#7c3aed", "Phase II": "#0891b2",
};

function Pill({ label, color }) {
  return (
    <span style={{
      display: "inline-block", fontSize: 10, fontWeight: 600, padding: "2px 8px",
      borderRadius: 20, background: color + "20", color, border: `1px solid ${color}40`,
      letterSpacing: "0.04em", textTransform: "uppercase",
    }}>{label}</span>
  );
}

function ParamInput({ name, label, value, hint, unit, onChange, isPercent, isSides, isEffectMode }) {
  if (isEffectMode) {
    const opts = [{ v: "means", l: "Means & SD" }, { v: "d", l: "Cohen's d" }];
    return (
      <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
        <label style={{ fontSize: 11, fontWeight: 600, color: "#64748b", letterSpacing: "0.03em" }}>{label}</label>
        <div style={{ display: "flex", gap: 6 }}>
          {opts.map(o => (
            <button key={o.v} onClick={() => onChange(name, o.v)} style={{
              flex: 1, padding: "7px 0", borderRadius: 8, border: "1.5px solid",
              borderColor: value === o.v ? "#6366f1" : "#e2e8f0",
              background: value === o.v ? "#eef2ff" : "#fff",
              color: value === o.v ? "#4f46e5" : "#475569",
              fontWeight: 600, fontSize: 11, cursor: "pointer",
            }}>{o.l}</button>
          ))}
        </div>
      </div>
    );
  }
  if (isSides) {
    return (
      <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
        <label style={{ fontSize: 11, fontWeight: 600, color: "#64748b", letterSpacing: "0.03em" }}>{label}</label>
        <div style={{ display: "flex", gap: 6 }}>
          {[1, 2].map(s => (
            <button key={s} onClick={() => onChange(name, s)} style={{
              flex: 1, padding: "7px 0", borderRadius: 8, border: "1.5px solid",
              borderColor: value === s ? "#6366f1" : "#e2e8f0",
              background: value === s ? "#eef2ff" : "#fff",
              color: value === s ? "#4f46e5" : "#475569",
              fontWeight: 600, fontSize: 12, cursor: "pointer",
            }}>{s === 2 ? "Two-sided" : "One-sided"}</button>
          ))}
        </div>
      </div>
    );
  }
  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
      <label style={{ fontSize: 11, fontWeight: 600, color: "#64748b", letterSpacing: "0.03em" }}>
        {label} {unit && <span style={{ color: "#94a3b8" }}>{unit}</span>}
      </label>
      <div style={{ position: "relative" }}>
        <input
          type="number" step="any" value={value}
          onChange={e => onChange(name, parseFloat(e.target.value) || 0)}
          style={{
            width: "100%", boxSizing: "border-box", padding: "8px 10px",
            borderRadius: 8, border: "1.5px solid #e2e8f0",
            fontSize: 13, fontFamily: "'SF Mono', 'Fira Code', monospace",
            color: "#1e293b", outline: "none", background: "#fafafa",
            transition: "border-color 0.15s",
          }}
          onFocus={e => e.target.style.borderColor = "#6366f1"}
          onBlur={e => e.target.style.borderColor = "#e2e8f0"}
        />
      </div>
      {hint && <p style={{ margin: 0, fontSize: 10.5, color: "#94a3b8", lineHeight: 1.4 }}>{hint}</p>}
    </div>
  );
}

function CodeBlock({ code, lang }) {
  const [copied, setCopied] = useState(false);
  const copy = () => {
    navigator.clipboard.writeText(code);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };
  const colors = { r: "#2563eb", stata: "#0891b2", sas: "#7c3aed" };
  return (
    <div style={{ borderRadius: 10, overflow: "hidden", border: "1px solid #e2e8f0" }}>
      <div style={{
        display: "flex", alignItems: "center", justifyContent: "space-between",
        padding: "6px 12px", background: colors[lang] || "#374151",
      }}>
        <span style={{ fontSize: 11, fontWeight: 700, color: "#fff", letterSpacing: "0.06em", textTransform: "uppercase" }}>{lang}</span>
        <button onClick={copy} style={{
          fontSize: 11, color: "#fff", background: "rgba(255,255,255,0.2)",
          border: "none", borderRadius: 5, padding: "2px 8px", cursor: "pointer",
        }}>{copied ? "✓ Copied" : "Copy"}</button>
      </div>
      <pre style={{
        margin: 0, padding: "12px 14px", background: "#0f172a",
        fontSize: 11.5, lineHeight: 1.7, color: "#e2e8f0",
        fontFamily: "'SF Mono', 'Fira Code', monospace", overflowX: "auto", whiteSpace: "pre",
      }}>{code}</pre>
    </div>
  );
}

function CitationCard({ cite, idx }) {
  return (
    <div style={{
      display: "flex", gap: 10, padding: "8px 10px",
      background: idx % 2 === 0 ? "#f8fafc" : "#fff",
      borderRadius: 6, alignItems: "flex-start",
    }}>
      <span style={{
        minWidth: 20, height: 20, borderRadius: "50%",
        background: "#e0e7ff", color: "#4338ca",
        fontSize: 10, fontWeight: 700, display: "flex", alignItems: "center",
        justifyContent: "center", flexShrink: 0,
      }}>{idx + 1}</span>
      <div>
        <div style={{ fontSize: 12, fontWeight: 600, color: "#1e293b" }}>{cite.ref}</div>
        <div style={{ fontSize: 11, color: "#6366f1", fontWeight: 500 }}>{cite.detail}</div>
        <div style={{ fontSize: 11, color: "#64748b" }}>{cite.note}</div>
      </div>
    </div>
  );
}

// ─────────────────────────────────────────────
// MAIN APP
// ─────────────────────────────────────────────
export default function StatCalcMD() {
  const [page, setPage] = useState("landing");
  const [mobileScreen, setMobileScreen] = useState("design");
  const [selectedKey, setSelectedKey] = useState("two_proportions");
  const [params, setParams] = useState({ ...FORMULAS["two_proportions"].defaults });
  const [result, setResult] = useState(null);
  const [activeTab, setActiveTab] = useState("results");
  const [activeCode, setActiveCode] = useState("r");
  const [copied, setCopied] = useState(false);
  const [activeGroup, setActiveGroup] = useState("All");

  const formula = FORMULAS[selectedKey];

  const selectFormula = useCallback((key) => {
    setSelectedKey(key);
    setParams({ ...FORMULAS[key].defaults });
    setResult(null);
    setActiveTab("results");
  }, []);

  const updateParam = useCallback((name, val) => {
    setParams(prev => ({ ...prev, [name]: val }));
  }, []);

  const calculate = useCallback(() => {
    try {
      const r = formula.calculate(params);
      setResult(r);
      setActiveTab("results");
    } catch (e) {
      setResult({ error: e.message });
    }
  }, [formula, params]);

  const goToParams = useCallback((key) => {
    selectFormula(key);
    setMobileScreen("params");
  }, [selectFormula]);

  const handleCalculate = useCallback(() => {
    try {
      const r = formula.calculate(params);
      setResult(r);
      setActiveTab("results");
      setMobileScreen("results");
    } catch (e) {
      setResult({ error: e.message });
      setMobileScreen("results");
    }
  }, [formula, params]);

  // Show landing page (after all hooks are declared — Rules of Hooks)
  if (page === "landing") {
    return <LandingPage onEnter={() => setPage("app")} />;
  }

  const manuscript = result && !result.error
    ? generateParagraph(selectedKey, params, result)
    : null;
  const code = result && !result.error && CODE[selectedKey]
    ? CODE[selectedKey](params, result)
    : null;
  const citations = CITATIONS[selectedKey] || [];

  const filteredKeys = Object.keys(FORMULAS).filter(k =>
    activeGroup === "All" || FORMULAS[k].group === activeGroup
  );

  const tabStyle = (tab) => ({
    padding: "7px 14px", borderRadius: 8, border: "none", cursor: "pointer",
    fontSize: 12, fontWeight: 600, transition: "all 0.15s",
    background: activeTab === tab ? "#4f46e5" : "transparent",
    color: activeTab === tab ? "#fff" : "#64748b",
  });

  return (
    <div style={{
      fontFamily: "'Inter', system-ui, sans-serif",
      background: "#f1f5f9",
      height: "100dvh",
      display: "flex",
      flexDirection: "column",
      overflow: "hidden",
      maxWidth: 430,
      margin: "0 auto",
      position: "relative",
    }}>

      {/* ── TOP NAV BAR ── */}
      <div style={{
        background: "#ffffff",
        borderBottom: "1px solid #e2e8f0",
        padding: "12px 16px 10px",
        display: "flex", alignItems: "center", justifyContent: "space-between",
        flexShrink: 0,
        paddingTop: "max(12px, env(safe-area-inset-top, 12px))",
      }}>
        <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
          {mobileScreen !== "design" && (
            <button
              onClick={() => setMobileScreen(mobileScreen === "results" ? "params" : "design")}
              style={{ background: "none", border: "none", color: "#1e3a5f", fontSize: 20, cursor: "pointer", padding: "0 4px 0 0", lineHeight: 1 }}
            >‹</button>
          )}
          <div style={{ width: 26, height: 26, borderRadius: 7, background: "#eef2ff", display: "flex", alignItems: "center", justifyContent: "center", fontSize: 13, fontWeight: 700, color: "#4f46e5" }}>Σ</div>
          <span style={{ fontSize: 15, fontWeight: 700, color: "#1e3a5f", letterSpacing: "-0.01em" }}>StatCalc MD</span>
        </div>
        <button
          onClick={() => setPage("landing")}
          style={{ background: "#eef2ff", border: "1px solid #c7d2fe", borderRadius: 7, padding: "4px 11px", color: "#4f46e5", fontSize: 11.5, fontWeight: 600, cursor: "pointer" }}
        >About</button>
      </div>

      {/* ── SCREEN CONTENT (scrollable) ── */}
      <div style={{ flex: 1, overflowY: "auto", WebkitOverflowScrolling: "touch" }}>

        {/* ── SCREEN 1: DESIGN PICKER ── */}
        {mobileScreen === "design" && (
          <div style={{ padding: "12px 14px 100px" }}>
            <div style={{ fontSize: 11, fontWeight: 700, color: "#94a3b8", letterSpacing: "0.07em", marginBottom: 10 }}>SELECT STUDY DESIGN</div>

            {/* Category filter chips */}
            <div style={{ display: "flex", gap: 6, overflowX: "auto", paddingBottom: 10, marginBottom: 4, WebkitOverflowScrolling: "touch" }}>
              {["All", ...GROUPS].map(g => (
                <button key={g} onClick={() => setActiveGroup(g)} style={{
                  flexShrink: 0, padding: "5px 12px", borderRadius: 16, border: "1px solid",
                  borderColor: activeGroup === g ? (TAG_COLORS[g] || "#6366f1") : "#e2e8f0",
                  background: activeGroup === g ? ((TAG_COLORS[g] || "#6366f1") + "18") : "#fff",
                  color: activeGroup === g ? (TAG_COLORS[g] || "#4f46e5") : "#64748b",
                  fontSize: 12, fontWeight: 600, cursor: "pointer", whiteSpace: "nowrap",
                }}>{g}</button>
              ))}
            </div>

            {/* Formula cards */}
            {GROUPS.filter(g => activeGroup === "All" || g === activeGroup).map(group => {
              const keys = filteredKeys.filter(k => FORMULAS[k].group === group);
              if (!keys.length) return null;
              const color = TAG_COLORS[group] || "#6366f1";
              return (
                <div key={group} style={{ marginBottom: 16 }}>
                  <div style={{ fontSize: 10.5, fontWeight: 700, color, letterSpacing: "0.07em", textTransform: "uppercase", marginBottom: 6, paddingLeft: 2 }}>{group}</div>
                  <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
                    {keys.map(key => (
                      <button key={key} onClick={() => goToParams(key)} style={{
                        display: "flex", alignItems: "center", justifyContent: "space-between",
                        width: "100%", padding: "13px 14px", borderRadius: 12, border: "1.5px solid",
                        borderColor: selectedKey === key ? color : "#e8edf2",
                        background: selectedKey === key ? color + "0d" : "#fff",
                        cursor: "pointer", textAlign: "left",
                        boxShadow: "0 1px 3px rgba(0,0,0,0.05)",
                      }}>
                        <span style={{ fontSize: 13.5, fontWeight: selectedKey === key ? 600 : 500, color: selectedKey === key ? color : "#1e293b" }}>
                          {FORMULAS[key].label}
                        </span>
                        <span style={{ fontSize: 18, color: selectedKey === key ? color : "#cbd5e1", flexShrink: 0 }}>›</span>
                      </button>
                    ))}
                  </div>
                </div>
              );
            })}
          </div>
        )}

        {/* ── SCREEN 2: PARAMETERS ── */}
        {mobileScreen === "params" && (
          <div style={{ padding: "16px 14px 100px" }}>
            {/* Design header */}
            <div style={{ display: "flex", alignItems: "center", gap: 8, marginBottom: 18 }}>
              <span style={{
                fontSize: 10.5, fontWeight: 700, padding: "3px 9px", borderRadius: 16,
                background: (TAG_COLORS[formula.group] || "#6366f1") + "18",
                color: TAG_COLORS[formula.group] || "#6366f1",
              }}>{formula.group}</span>
              <span style={{ fontSize: 14, fontWeight: 700, color: "#1e293b" }}>{formula.label}</span>
            </div>

            {/* Parameter inputs */}
            <div style={{ display: "flex", flexDirection: "column", gap: 14 }}>
              {(formula.visibleParams ? formula.visibleParams(params) : formula.params).map(name => (
                <ParamInput
                  key={name} name={name}
                  label={formula.labels[name]}
                  value={params[name]}
                  hint={formula.hints?.[name]}
                  unit={formula.units?.[name]}
                  onChange={updateParam}
                  isSides={name === "sides"}
                  isEffectMode={name === "effectMode"}
                />
              ))}
            </div>

            {/* Calculate button */}
            <button onClick={handleCalculate} style={{
              width: "100%", marginTop: 24, padding: "15px 0",
              background: "linear-gradient(135deg, #4f46e5, #6366f1)",
              color: "#fff", border: "none", borderRadius: 14,
              fontSize: 15, fontWeight: 700, cursor: "pointer",
              letterSpacing: "0.01em",
              boxShadow: "0 6px 20px rgba(99,102,241,0.40)",
            }}>
              Calculate Sample Size →
            </button>
          </div>
        )}

        {/* ── SCREEN 3: RESULTS ── */}
        {mobileScreen === "results" && (
          <div style={{ padding: "14px 14px 100px" }}>
            {!result ? (
              <div style={{ textAlign: "center", padding: "60px 20px", color: "#94a3b8" }}>
                <div style={{ fontSize: 40, marginBottom: 12 }}>Σ</div>
                <div style={{ fontSize: 14, fontWeight: 500 }}>No result yet</div>
                <div style={{ fontSize: 12, marginTop: 6 }}>Go to Parameters and tap Calculate</div>
                <button onClick={() => setMobileScreen("params")} style={{ marginTop: 16, padding: "10px 24px", borderRadius: 10, background: "#4f46e5", color: "#fff", border: "none", fontSize: 13, fontWeight: 700, cursor: "pointer" }}>
                  Set Parameters →
                </button>
              </div>
            ) : result.error ? (
              <div style={{ padding: 20, color: "#ef4444", background: "#fff", borderRadius: 12, border: "1px solid #fecaca" }}>
                ⚠ Error: {result.error}
              </div>
            ) : (
              <>
                {/* Design label */}
                <div style={{ fontSize: 12, color: "#64748b", marginBottom: 10, fontWeight: 500 }}>
                  {formula.label}
                </div>

                {/* Result cards — 2-column grid */}
                <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 8, marginBottom: 12 }}>
                  {result.n2 != null ? (
                    <>
                      <div style={cardStyle("#4f46e5")}>
                        <div style={cardLabel}>Group 1 (n)</div>
                        <div style={cardVal("#4f46e5")}>{result.n1?.toLocaleString()}</div>
                      </div>
                      <div style={cardStyle("#4f46e5")}>
                        <div style={cardLabel}>Group 2 (n)</div>
                        <div style={cardVal("#4f46e5")}>{result.n2?.toLocaleString()}</div>
                      </div>
                      <div style={{ ...cardStyle("#0891b2"), gridColumn: "span 2", background: "#eff6ff", border: "2px solid #0891b220" }}>
                        <div style={cardLabel}>Total N</div>
                        <div style={{ ...cardVal("#0891b2"), fontSize: 32 }}>{result.total?.toLocaleString()}</div>
                      </div>
                    </>
                  ) : (
                    <div style={{ ...cardStyle("#0891b2"), gridColumn: "span 2", background: "#eff6ff", border: "2px solid #0891b220" }}>
                      <div style={cardLabel}>Total N</div>
                      <div style={{ ...cardVal("#0891b2"), fontSize: 32 }}>{result.total?.toLocaleString()}</div>
                    </div>
                  )}
                  {result.perGroup && result.n2 === null && (
                    <div style={cardStyle("#4f46e5")}>
                      <div style={cardLabel}>Per group</div>
                      <div style={cardVal("#4f46e5")}>{result.perGroup?.toLocaleString()}</div>
                    </div>
                  )}
                  {result.events && (
                    <div style={cardStyle("#f59e0b")}>
                      <div style={cardLabel}>Events needed</div>
                      <div style={cardVal("#f59e0b")}>{result.events}</div>
                    </div>
                  )}
                  {result.effectSize != null && result.effectLabel && (
                    <div style={cardStyle("#10b981")}>
                      <div style={cardLabel}>{result.effectLabel}</div>
                      <div style={{ ...cardVal("#10b981"), fontSize: 18 }}>{typeof result.effectSize === "number" ? result.effectSize.toFixed(3) : result.effectSize}</div>
                    </div>
                  )}
                </div>

                {/* Note */}
                {result.note && (
                  <div style={{ padding: "10px 12px", background: "#f0f9ff", borderRadius: 10, border: "1px solid #bae6fd", marginBottom: 12, fontSize: 12, color: "#0c4a6e", lineHeight: 1.5 }}>
                    ℹ {result.note}
                  </div>
                )}

                {/* Output tabs */}
                <div style={{ display: "flex", gap: 0, background: "#f1f5f9", borderRadius: 12, padding: 4, marginBottom: 14 }}>
                  {[
                    { k: "results", icon: "📊", label: "Interpret" },
                    { k: "summary", icon: "📝", label: "Summary" },
                    { k: "code",    icon: "💻", label: "Code" },
                    { k: "citations", icon: "📚", label: "Refs" },
                  ].map(t => (
                    <button key={t.k} onClick={() => setActiveTab(t.k)} style={{
                      flex: 1, padding: "8px 4px", borderRadius: 9, border: "none", cursor: "pointer",
                      background: activeTab === t.k ? "#4f46e5" : "transparent",
                      color: activeTab === t.k ? "#fff" : "#64748b",
                      fontSize: 10.5, fontWeight: 700, display: "flex", flexDirection: "column", alignItems: "center", gap: 2,
                    }}>
                      <span style={{ fontSize: 14 }}>{t.icon}</span>
                      <span>{t.label}</span>
                    </button>
                  ))}
                </div>

                {/* Tab content */}
                {activeTab === "results" && (
                  <div style={{ background: "#fff", borderRadius: 12, padding: 16, border: "1px solid #e2e8f0" }}>
                    <div style={{ fontSize: 13, fontWeight: 600, color: "#1e293b", marginBottom: 10 }}>Interpretation</div>
                    <InterpretationText formulaKey={selectedKey} params={params} result={result} />
                  </div>
                )}

                {activeTab === "summary" && (
                  <SummaryTab manuscript={manuscript} formula={formula} params={params} result={result} />
                )}

                {activeTab === "code" && code && (
                  <div style={{ display: "flex", flexDirection: "column", gap: 10 }}>
                    <div style={{ display: "flex", gap: 6 }}>
                      {["r", "stata", "sas"].map(lang => {
                        const col = { r: "#2563eb", stata: "#0891b2", sas: "#7c3aed" };
                        return (
                          <button key={lang} onClick={() => setActiveCode(lang)} style={{
                            flex: 1, padding: "8px 4px", borderRadius: 8, border: "1.5px solid",
                            borderColor: activeCode === lang ? col[lang] : "#e2e8f0",
                            background: activeCode === lang ? col[lang] + "15" : "#fff",
                            color: activeCode === lang ? col[lang] : "#64748b",
                            fontSize: 12, fontWeight: 700, cursor: "pointer", textTransform: "uppercase",
                          }}>{lang}</button>
                        );
                      })}
                    </div>
                    <CodeBlock code={code[activeCode]} lang={activeCode} />
                    <div style={{ fontSize: 11, color: "#94a3b8" }}>
                      Run this code to independently verify the calculation.
                    </div>
                  </div>
                )}

                {activeTab === "citations" && (
                  <div style={{ background: "#fff", borderRadius: 12, padding: 14, border: "1px solid #e2e8f0" }}>
                    <div style={{ fontSize: 13, fontWeight: 600, color: "#1e293b", marginBottom: 10 }}>
                      References — {formula.label}
                    </div>
                    <div style={{ display: "flex", flexDirection: "column", gap: 3 }}>
                      {citations.map((c, i) => <CitationCard key={i} cite={c} idx={i} />)}
                    </div>
                    {citations.length > 0 && (
                      <div style={{ marginTop: 10, padding: "9px 11px", background: "#f0fdf4", borderRadius: 8, border: "1px solid #bbf7d0", fontSize: 11, color: "#166534" }}>
                        💡 Citation [1] is the primary reference for the manuscript. Use [2–{citations.length}] as backup.
                      </div>
                    )}
                  </div>
                )}
              </>
            )}
          </div>
        )}
      </div>

      {/* ── BOTTOM TAB BAR ── */}
      <div style={{
        display: "flex", background: "#fff",
        borderTop: "1px solid #e2e8f0",
        paddingBottom: "max(12px, env(safe-area-inset-bottom, 12px))",
        flexShrink: 0,
        boxShadow: "0 -4px 20px rgba(0,0,0,0.07)",
      }}>
        {[
          { screen: "design",  icon: "⊞", label: "Designs" },
          { screen: "params",  icon: "⚙", label: "Parameters" },
          { screen: "results", icon: "Σ", label: "Results" },
        ].map(({ screen, icon, label }) => {
          const active = mobileScreen === screen;
          const hasResult = screen === "results" && result && !result.error;
          return (
            <button key={screen} onClick={() => setMobileScreen(screen)} style={{
              flex: 1, display: "flex", flexDirection: "column", alignItems: "center",
              gap: 3, padding: "10px 4px 4px", background: "none", border: "none", cursor: "pointer",
              position: "relative",
            }}>
              {hasResult && !active && (
                <span style={{ position: "absolute", top: 8, right: "calc(50% - 14px)", width: 7, height: 7, borderRadius: "50%", background: "#10b981" }} />
              )}
              <span style={{ fontSize: 18, lineHeight: 1, color: active ? "#4f46e5" : "#94a3b8" }}>{icon}</span>
              <span style={{ fontSize: 10.5, fontWeight: 700, color: active ? "#4f46e5" : "#94a3b8", letterSpacing: "0.01em" }}>{label}</span>
              {active && (
                <span style={{ position: "absolute", top: 0, left: "20%", right: "20%", height: 2, background: "#4f46e5", borderRadius: "0 0 2px 2px" }} />
              )}
            </button>
          );
        })}
      </div>
    </div>
  );
}

// Card style helpers
const cardStyle = (accent) => ({
  background: "#fff", borderRadius: 12, padding: "12px 14px",
  border: `1.5px solid ${accent}22`,
});
const cardLabel = { fontSize: 10, fontWeight: 700, color: "#94a3b8", letterSpacing: "0.06em", textTransform: "uppercase", marginBottom: 4 };
const cardVal = (accent) => ({ fontSize: 24, fontWeight: 700, color: accent, fontVariantNumeric: "tabular-nums" });



function ResultCard({ label, value, accent, large, isText }) {
  return (
    <div style={{
      background: "#fff", borderRadius: 12, padding: "12px 16px",
      border: `1.5px solid ${accent}30`,
      boxShadow: large ? `0 0 0 3px ${accent}15` : "none",
      minWidth: 110,
    }}>
      <div style={{ fontSize: 10, fontWeight: 600, color: "#94a3b8", letterSpacing: "0.05em", marginBottom: 4, textTransform: "uppercase" }}>{label}</div>
      <div style={{ fontSize: large ? 28 : 20, fontWeight: 700, color: accent, fontVariantNumeric: "tabular-nums" }}>
        {isText ? value : Number(value).toLocaleString()}
      </div>
    </div>
  );
}

function InterpretationText({ formulaKey, params, result }) {
  const effectBenchmarks = {
    "Cohen's d": [[0.2, "small"], [0.5, "medium"], [0.8, "large"]],
    "Cohen's d (paired)": [[0.2, "small"], [0.5, "medium"], [0.8, "large"]],
    "Cohen's h": [[0.2, "small"], [0.5, "medium"], [0.8, "large"]],
    "Cohen's f": [[0.1, "small"], [0.25, "medium"], [0.4, "large"]],
    "Cohen's f²": [[0.02, "small"], [0.15, "medium"], [0.35, "large"]],
    "ρ (Pearson)": [[0.1, "small"], [0.3, "medium"], [0.5, "large"]],
  };
  const benchmarks = result.effectLabel ? effectBenchmarks[result.effectLabel] : null;
  let effectCat = "";
  if (benchmarks && result.effectSize) {
    for (const [threshold, label] of benchmarks) {
      if (result.effectSize >= threshold) effectCat = label;
    }
  }

  return (
    <div style={{ fontSize: 12.5, lineHeight: 1.75, color: "#475569" }}>
      {result.n2 !== null && result.n2 !== undefined && (
        <p style={{ margin: "0 0 8px" }}>
          You need <strong style={{ color: "#4f46e5" }}>{result.n1}</strong> participants in group 1 and{" "}
          <strong style={{ color: "#4f46e5" }}>{result.n2}</strong> in group 2,
          for a <strong style={{ color: "#0891b2" }}>total of {result.total}</strong>.
          {params.attrition > 0 && ` This already accounts for ${params.attrition}% attrition.`}
        </p>
      )}
      {(result.n2 === null || result.n2 === undefined) && (
        <p style={{ margin: "0 0 8px" }}>
          You need a <strong style={{ color: "#0891b2" }}>total of {result.total}</strong> participants.
          {params.attrition > 0 && ` This already accounts for ${params.attrition}% attrition.`}
        </p>
      )}
      {result.effectSize !== null && result.effectSize !== undefined && benchmarks && (
        <p style={{ margin: "0 0 8px" }}>
          The effect size ({result.effectLabel} = {result.effectSize.toFixed(3)}) is considered{" "}
          <strong style={{ color: "#10b981" }}>{effectCat || "very small"}</strong> by Cohen's (1988) conventions.
          {effectCat === "small" && " Consider whether this is the minimum clinically important difference."}
          {effectCat === "medium" && " This is a realistic and detectable effect for most clinical studies."}
          {effectCat === "large" && " A large effect — verify this reflects a realistic clinical difference."}
        </p>
      )}
      {result.events && (
        <p style={{ margin: "0 0 8px" }}>
          The study requires <strong style={{ color: "#f59e0b" }}>{result.events} events</strong>.
          With an event probability of {(params.p_event * 100).toFixed(0)}%, you need {result.total} participants enrolled.
        </p>
      )}
      <p style={{ margin: 0, color: "#94a3b8", fontSize: 11.5 }}>
        Alpha = {params.alpha} · Power = {params.power * 100}%
        {params.sides ? ` · ${params.sides === 2 ? "Two-sided" : "One-sided"} test` : ""}
      </p>
    </div>
  );
}
