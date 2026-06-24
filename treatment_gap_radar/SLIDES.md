---
marp: true
theme: default
paginate: true
title: Treatment Gap Radar
---

<!-- Render: `marp SLIDES.md --pdf` (or open in VS Code with the Marp extension). -->

# 🛰️ Treatment Gap Radar
### Does global antimicrobial R&D go where resistance is actually getting worse?

**Vivli AMR Surveillance Open Data Re-use Data Challenge** · AMR ID 00013367

Integrating 13 surveillance datasets with the Global AMR R&D Hub to find where
resistance need and R&D attention are misaligned.

---

## The problem

- Resistance keeps rising; antibiotic R&D is slow, expensive, and **finite**.
- Funders and companies concentrate on a few "famous" pathogens.
- **No one routinely checks, with data, whether R&D goes where resistance is worst.**

→ We built a radar that does exactly that.

---

## The idea — two scores per pathogen

| Index | Question | Source |
|---|---|---|
| **Resistance Need Index (RNI)** | How bad is this bug getting? | 6 surveillance indicators |
| **R&D Attention Index (RAI)** | How much is the world investing? | Global AMR R&D Hub |
| **Gap = RNI − RAI** | Is attention aligned with need? | the overlay |

Six indicators: prevalence · time-trend · multidrug resistance · geographic spread ·
therapeutic scarcity · pediatric involvement.

---

## The data

- **13 Vivli AMR datasets** — ATLAS, GEARS, SIDERO-WT, KEYSTONE, SOAR×3, Innoviva,
  DREAM (TB), Venus (PLEA/GASAR), SPIDAAR (Africa RWE).
- **~1.0 M isolates → 17.9 M isolate–antibiotic results (2004–2025).**
- **Global AMR R&D Hub** — 19,023 projects, **$18.9 B**.

One harmonized table; CLSI/TB breakpoints derive S/I/R for the 12 MIC-only datasets.

---

## It's validated

- **Cross-dataset:** *A. baumannii* meropenem 60.9% (ATLAS) vs 62.3% (Innoviva);
  *P. aeruginosa* meropenem 20–21% across ATLAS & GEARS.
- **vs WHO GLASS 2022:** MRSA **36.6% vs 35%** (near-exact). 3GC estimates sit below GLASS
  medians — consistent with ATLAS's known high-income bias; our Africa subset (82–90%) brackets
  the high end.

---

## The result — the gap

🔴 **Priority gaps (high need, low attention):** **Enterococcus faecium (VRE)**, *S. epidermidis*.
🟢 **Well-served:** *A. baumannii*, *K. pneumoniae*, **TB** (high need + high R&D).
🔵 Broad-spectrum Gram-negative R&D already covers the Enterobacterales.

*Robust:* identical priority gaps under both defensible R&D-attribution schemes (Spearman ≈ 1).

---

## It's rigorous, not arbitrary

- **PCA:** the six indicators load on one axis (PC1 = 72% variance) → equal weighting is *data-justified*.
- **Sensitivity:** rankings stable across weighting schemes (Spearman 0.98–0.99).
- **Trends:** logistic resistance~year, e.g. *A. baumannii* meropenem **OR 1.07/yr, p<0.001**.
- **Bootstrap 95% CIs** on every headline rate.

---

## Early warning — forecasting the future

**Klebsiella pneumoniae carbapenem resistance is projected to cross 50% by ≈2036**
(95% CI 2035–2038), rising from ~21% today.

*A. baumannii* already >50%, heading to ~86% by 2035. Declining drugs flagged as no-crossing.

→ A forward-looking trigger for stewardship and procurement.

---

## Predicting the blind spots (ML)

Gradient-boosted model predicts resistance from **generalizable** features (pathogen, drug, Gram,
**continent, income tier**, year) — never country identity.

- **R² = 0.73 on entirely unseen countries** (weighted MAE 8.1% vs 18.4% baseline).
- → estimates resistance where there is **no surveillance** (esp. LMICs).

---

## Findings & implications

1. **Targeted R&D for Gram-positive gaps** — VRE / *E. faecium*.
2. **Sustain investment** in critical Gram-negatives — resistance still climbing.
3. **Expand LMIC surveillance** — the biggest blind spot; Africa shows >80% 3GC resistance
   yet is barely monitored.

---

## Deliverables & openness

- 📓 Reproducible notebook · 🛰️ interactive 9-tab dashboard · 📄 full report
- One command: `python -m src.pipeline` · open GitHub repo
- Raw data kept local (Vivli DUA); only aggregate results published (n ≥ 30)

> Data from GSK, Innoviva, Johnson & Johnson, Paratek, Pfizer, Shionogi, Venatorx,
> Venus Remedies, obtained through https://amr.vivli.org

---

# Thank you

**Treatment Gap Radar** — aligning antimicrobial R&D with where resistance is worsening.

github.com/niyeldeii/AMR-VIVLI--treatment-gap-radar
