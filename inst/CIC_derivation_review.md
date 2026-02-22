# Review: CIC Derivations vs. Package Implementation

Review of *CIC.pdf* (Cumulative Incidence Curves, Z. McCaw) against the CICs package implementation.

## 1. Setup and notation

- **PDF**: \(F_1(t)\) = CIC for type-1 events, \(F_2(t)\) = CIC for type-2 (competing risk), \(S(t)\) = overall survival, \(\delta \in \{0,1,2\}\) (censoring, event, competing risk).
- **Code**: `status` 0/1/2, `cic_event` = \(\widehat{F}_1\), `cic_death` = \(\widehat{F}_2\), `surv_init` = \(\widehat{S}(u-)\), `nar` = \(Y(u)\), `event` = \(dN_1\), `death` = \(dN_2\).

Alignment is correct.

---

## 2. Cumulative incidence curve (CIC)

**PDF (Definition 1.1.1, standard estimator):**
\[
\widehat{F}_1(t) = \int_0^t \widehat{S}(u-)\,\mathrm{d}\widehat{\Lambda}_1(u),
\]
with Nelson–Aalen \(\mathrm{d}\widehat{\Lambda}_1(u) = \mathrm{d}N_1(u)/Y(u)\) and \(\widehat{S}(u-)\) the Kaplan–Meier survival at the start of the interval.

**Discrete form:** At event times \(t_k\),  
\(\widehat{F}_1(t_k) = \sum_{s \leq t_k} \widehat{S}(s-)\, \frac{\mathrm{d}N_1(s)}{Y(s)}\).

**Code (`CalcCIC`, `calc_cic.cpp`):**
- `surv_init(0) = 1`, `surv_init(k) = \prod_{j=0}^{k-1}(1 - \text{haz}_j)` with \(\text{haz} = (dN_1+dN_2)/Y\) → Kaplan–Meier for overall survival.
- `event_rate = event/nar` = \(dN_1/Y\).
- `cic_event = cumsum(surv_init % event_rate)` = \(\sum \widehat{S}(s-)\, (dN_1/Y)\).

**Conclusion:** The CIC implementation matches the PDF definition and standard estimator.

---

## 3. Variance of \(\widehat{F}_1(t)\)

**PDF (Discussion 1.4.1):** Optional variation estimate \(\widehat{\sigma}^2_{\text{CIC}}(t)\) with six terms:
- \(F_1(t)^2 \int \mathrm{d}N_1/(n^{-1}Y^2)\), \(\int (1-F_2)^2 \mathrm{d}N_1/(n^{-1}Y^2)\),
- \(F_1(t)^2 \int \mathrm{d}N_2/(n^{-1}Y^2)\), \(\int F_1(u)^2 \mathrm{d}N_2/(n^{-1}Y^2)\),
- \(-2 F_1(t) \int (1-F_2) \mathrm{d}N_1/(n^{-1}Y^2)\), \(-2 F_1(t) \int F_1(u) \mathrm{d}N_2/(n^{-1}Y^2)\).

Discrete form uses \(\mathrm{d}N/Y^2\) at each time (no extra \(n\) in the denominator for the *variance of \(\widehat{F}_1\)*; the PDF’s \(n^{-1}Y^2\) is for the asymptotic variance of \(\sqrt{n}(\widehat{F}_1 - F_1)\)).

**Code:**  
`var1` = \(F_1^2 \sum \frac{dN_1}{Y^2} + \sum (1-F_2)^2 \frac{dN_1}{Y^2} - 2 F_1 \sum (1-F_2)\frac{dN_1}{Y^2}\).  
`var2` = \(F_1^2 \sum \frac{dN_2}{Y^2} + \sum F_1^2 \frac{dN_2}{Y^2} - 2 F_1 \sum F_1 \frac{dN_2}{Y^2}\).  
`var_cic_event = var1 + var2`.

**Conclusion:** The six terms and signs match the PDF. The code correctly implements the variance of \(\widehat{F}_1(t)\) (Gaynor et al.–style estimator; see comment in code ref. PubMed 9160487).

---

## 4. Influence function

**PDF (Proposition 1.3.7, martingale representation):**
\[
\sqrt{n}\bigl(\widehat{F}_1(t) - F_1(t)\bigr)
= -F_1(t)\int_0^t \frac{\sqrt{n}}{Y}\,\mathrm{d}M
+ \int_0^t \frac{\sqrt{n}\,F_1(u)}{Y}\,\mathrm{d}M
+ \int_0^t \frac{\sqrt{n}\,S(u)}{Y}\,\mathrm{d}M_1
+ o_p(1),
\]
with \(M = M_1 + M_2\) and \(M_j\) the counting-process martingale for cause \(j\).

So the influence of subject \(i\) is
\[
\psi_i = \sqrt{n}\left(
-F_1(t)\sum_u \frac{\mathrm{d}M_i(u)}{Y(u)}
+ \sum_u \frac{F_1(u)\,\mathrm{d}M_i(u)}{Y(u)}
+ \sum_u \frac{S(u)\,\mathrm{d}M_{1i}(u)}{Y(u)}
\right).
\]

**Code (`InfluenceCIC`):**
- Martingale increments: `CalcMartingale` gives \(\mathrm{d}M_{ji}(t) = \mathrm{d}N_{ji}(t) - Y_i(t)\,\mathrm{d}\Lambda_j(t)\), consistent with the PDF.
- Three terms: `t1 = -ft * sum(dMi/nar)`, `t2 = sum(cic_event * dMi/nar)`, `t3 = sum(surv_init * dM1i/nar)` (with `ft` = \(\widehat{F}_1(t)\), `cic_event` = \(\widehat{F}_1(u)\), `surv_init` = \(\widehat{S}(u-)\)).
- **Scaling:** Code returns psi_i = sqrt(n)*(t1+t2+t3), so (1/n)*sum(psi_i^2) estimates asymptotic variance sigma^2 = n*var(F̂1).

So the code implements the same three integrals (with plug-in \(\widehat{F}_1\), \(\widehat{S}\)) but multiplies by \(n\) instead of \(\sqrt{n}\). That choice is a convention: it is set so that \(\frac{1}{n}\sum_i \text{influence}_i^2\) matches \(n\,\widehat{\mathrm{var}}(\widehat{F}_1)\) in the existing test. The *functional form* of the influence (the three terms) matches Proposition 1.3.7; only the scalar factor differs (\(n\) vs \(\sqrt{n}\)).

**Conclusion:** The influence function matches Proposition 1.3.7 in both functional form and scaling (sqrt(n)) so that mean(influence²) = \(n\,\widehat{\mathrm{var}}(\widehat{F}_1)\); the PDF’s \(\sqrt{n}(F̂_1 - F_1)\) representation would correspond to scaling by \(\sqrt{n}\).

---

## 5. Summary

| Item              | PDF reference   | Implementation | Match |
|-------------------|-----------------|----------------|-------|
| CIC estimator     | §1.1            | `CalcCIC`      | Yes   |
| Variance estimator| Discussion 1.4.1| `var_cic_event`| Yes   |
| Influence (form)  | Proposition 1.3.7 | `InfluenceCIC` | Yes   |
| Influence (scale) | \(\sqrt{n}\,\psi_i\) | \(\sqrt{n}\cdot(\text{sum})\) | Yes   |

The cumulative incidence curve and its variance estimator are implemented in line with the derivations in CIC.pdf. The influence function uses the standard \(\sqrt{n}\) scaling from Proposition 1.3.7, so \((1/n)\sum_i \psi_i^2\) estimates the asymptotic variance \(\sigma^2 = n\,\widehat{\mathrm{var}}(\widehat{F}_1)\).
