import pandas as pd
import numpy as np
from scipy import optimize, stats
import statsmodels.api as sm

df = pd.read_csv('meta_analysis_jegh_preprint/data/raw/part_a_african_cohorts.csv')
df_eaec = df[(df['eaec_pos'] > 0) & (df['eaec_total'] > 0)].copy()

# 1. Modality: WGS vs PCR vs Culture
def get_modality(m):
    m_lower = str(m).lower()
    if 'wgs' in m_lower:
        return 'WGS'
    elif 'colony' in m_lower or 'culture' in m_lower:
        return 'Culture'
    else:
        return 'PCR'

df_eaec['modality'] = df_eaec['modality'].apply(get_modality)

# 2. Era: Pre-2015 vs 2015-2024
import re
def get_era(study_id):
    match = re.search(r'20\d\d|19\d\d', study_id)
    if match:
        year = int(match.group(0))
        return 'Pre-2015' if year < 2015 else '2015-2024'
    return '2015-2024'

df_eaec['era_cat'] = df_eaec['study_id'].apply(get_era)

# 3. Age: Pediatric vs All_Ages
df_eaec['age_cat'] = df_eaec['age_strata'].apply(lambda x: 'All_Ages' if ('mixed' in str(x).lower() or 'all' in str(x).lower()) else 'Pediatric')

# Calculate effect size: logit prevalence
# p = (r + 0.5) / (n + 1)
# y = log(p / (1 - p))
# v = 1 / (n * p * (1 - p))
r = df_eaec['eaec_pos'].values
n = df_eaec['eaec_total'].values

# Freeman-Tukey arcsine transformed proportion
# t = arcsin(sqrt(r/(n+1))) + arcsin(sqrt((r+1)/(n+1)))
# v = 1 / (n + 0.5)
# Or logit:
p_corr = (r + 0.5) / (n + 1)
y_logit = np.log(p_corr / (1 - p_corr))
v_logit = 1.0 / (r + 0.5) + 1.0 / (n - r + 0.5)

# FT:
t_ft = np.arcsin(np.sqrt(r / (n + 1.0))) + np.arcsin(np.sqrt((r + 1.0) / (n + 1.0)))
v_ft = 1.0 / (n + 0.5)

def run_meta_reg(y, v, X_df, name=""):
    X = sm.add_constant(X_df).values
    k, p = X.shape
    
    # REML function
    def reml_objective(tau2):
        if tau2 < 0:
            return 1e12
        w = 1.0 / (v + tau2)
        W = np.diag(w)
        XtWX = X.T @ W @ X
        sign, logdet = np.linalg.slogdet(XtWX)
        if sign <= 0:
            return 1e12
        # beta hat
        beta = np.linalg.solve(XtWX, X.T @ (w * y))
        res = y - X @ beta
        # REML log likelihood
        loglik = -0.5 * (np.sum(np.log(v + tau2)) + logdet + np.sum(w * (res**2)))
        return -loglik # minimize negative loglik
    
    # Optimize tau2
    res_opt = optimize.minimize_scalar(reml_objective, bounds=(0, 2.0), method='bounded')
    tau2 = max(0.0, float(res_opt.x))
    
    # Final model
    w = 1.0 / (v + tau2)
    W = np.diag(w)
    XtWX = X.T @ W @ X
    inv_XtWX = np.linalg.inv(XtWX)
    beta = inv_XtWX @ (X.T @ (w * y))
    se = np.sqrt(np.diag(inv_XtWX))
    z = beta / se
    pval = 2.0 * (1.0 - stats.norm.cdf(np.abs(z)))
    
    # Omnibus test Q_M for moderators (excluding intercept)
    if p > 1:
        beta_mods = beta[1:]
        cov_mods = inv_XtWX[1:, 1:]
        QM = float(beta_mods.T @ np.linalg.inv(cov_mods) @ beta_mods)
        df_QM = p - 1
        p_QM = 1.0 - stats.chi2.cdf(QM, df_QM)
    else:
        QM, df_QM, p_QM = 0.0, 0, 1.0
        
    # Residual heterogeneity QE
    res = y - X @ beta
    QE = float(np.sum(w * (res**2)))
    df_QE = k - p
    p_QE = 1.0 - stats.chi2.cdf(QE, df_QE) if df_QE > 0 else 1.0
    
    # I2 residual
    # I2 = (QE - df_QE) / QE
    I2 = max(0.0, (QE - df_QE) / QE) * 100.0 if QE > df_QE else 0.0
    
    return {
        "name": name,
        "tau2": tau2,
        "I2": I2,
        "QE": QE,
        "p_QE": p_QE,
        "QM": QM,
        "df_QM": df_QM,
        "p_QM": p_QM,
        "beta": beta,
        "se": se,
        "pval": pval,
        "coef_names": ['Intercept'] + list(X_df.columns)
    }

print("=== LOGIT TRANSFORMED PREVALENCE ===")
# Total variance (intercept only)
null_mod = run_meta_reg(y_logit, v_logit, pd.DataFrame(index=df_eaec.index), "Null")
tau2_0 = null_mod['tau2']
print(f"Null Model: tau2_0 = {tau2_0:.5f}, I2 = {null_mod['I2']:.2f}%")

# 1. Modality
X_mod = pd.get_dummies(df_eaec['modality'], drop_first=True, dtype=float)
res_mod = run_meta_reg(y_logit, v_logit, X_mod, "Modality (Culture, WGS vs PCR)")
r2_mod = max(0.0, (tau2_0 - res_mod['tau2']) / tau2_0) * 100.0
print(f"\nModality alone: tau2 = {res_mod['tau2']:.5f}, R2 = {r2_mod:.2f}%, QM = {res_mod['QM']:.3f} (p={res_mod['p_QM']:.4f})")
for name, b, s, p in zip(res_mod['coef_names'], res_mod['beta'], res_mod['se'], res_mod['pval']):
    print(f"  {name}: beta={b:.4f}, se={s:.4f}, p={p:.4f}")

# 2. Era
X_era = pd.get_dummies(df_eaec['era_cat'], drop_first=True, dtype=float)
res_era = run_meta_reg(y_logit, v_logit, X_era, "Era (Pre-2015 vs 2015-2024)")
r2_era = max(0.0, (tau2_0 - res_era['tau2']) / tau2_0) * 100.0
print(f"\nEra alone: tau2 = {res_era['tau2']:.5f}, R2 = {r2_era:.2f}%, QM = {res_era['QM']:.3f} (p={res_era['p_QM']:.4f})")
for name, b, s, p in zip(res_era['coef_names'], res_era['beta'], res_era['se'], res_era['pval']):
    print(f"  {name}: beta={b:.4f}, se={s:.4f}, p={p:.4f}")

# 3. Age
X_age = pd.get_dummies(df_eaec['age_cat'], drop_first=True, dtype=float)
res_age = run_meta_reg(y_logit, v_logit, X_age, "Age (Pediatric vs All Ages)")
r2_age = max(0.0, (tau2_0 - res_age['tau2']) / tau2_0) * 100.0
print(f"\nAge alone: tau2 = {res_age['tau2']:.5f}, R2 = {r2_age:.2f}%, QM = {res_age['QM']:.3f} (p={res_age['p_QM']:.4f})")
for name, b, s, p in zip(res_age['coef_names'], res_age['beta'], res_age['se'], res_age['pval']):
    print(f"  {name}: beta={b:.4f}, se={s:.4f}, p={p:.4f}")

# 4. Multivariate model: Modality + Era + Age
X_full = pd.concat([X_mod, X_era, X_age], axis=1)
res_full = run_meta_reg(y_logit, v_logit, X_full, "Full Multivariate Model")
r2_full = max(0.0, (tau2_0 - res_full['tau2']) / tau2_0) * 100.0
print(f"\nMultivariate Model: tau2 = {res_full['tau2']:.5f}, R2 = {r2_full:.2f}%, QM = {res_full['QM']:.3f} (p={res_full['p_QM']:.4f}), residual I2 = {res_full['I2']:.2f}%")
for name, b, s, p in zip(res_full['coef_names'], res_full['beta'], res_full['se'], res_full['pval']):
    print(f"  {name}: beta={b:.4f}, se={s:.4f}, p={p:.4f}")


print("\n=== FREEMAN-TUKEY ARCSINE TRANSFORMED PREVALENCE ===")
null_mod_ft = run_meta_reg(t_ft, v_ft, pd.DataFrame(index=df_eaec.index), "Null FT")
tau2_0_ft = null_mod_ft['tau2']
print(f"Null Model FT: tau2_0 = {tau2_0_ft:.5f}")

res_mod_ft = run_meta_reg(t_ft, v_ft, X_mod, "Modality FT")
r2_mod_ft = max(0.0, (tau2_0_ft - res_mod_ft['tau2']) / tau2_0_ft) * 100.0
print(f"Modality alone FT: tau2 = {res_mod_ft['tau2']:.5f}, R2 = {r2_mod_ft:.2f}%, QM = {res_mod_ft['QM']:.3f} (p={res_mod_ft['p_QM']:.4f})")

res_era_ft = run_meta_reg(t_ft, v_ft, X_era, "Era FT")
r2_era_ft = max(0.0, (tau2_0_ft - res_era_ft['tau2']) / tau2_0_ft) * 100.0
print(f"Era alone FT: tau2 = {res_era_ft['tau2']:.5f}, R2 = {r2_era_ft:.2f}%, QM = {res_era_ft['QM']:.3f} (p={res_era_ft['p_QM']:.4f})")

res_age_ft = run_meta_reg(t_ft, v_ft, X_age, "Age FT")
r2_age_ft = max(0.0, (tau2_0_ft - res_age_ft['tau2']) / tau2_0_ft) * 100.0
print(f"Age alone FT: tau2 = {res_age_ft['tau2']:.5f}, R2 = {r2_age_ft:.2f}%, QM = {res_age_ft['QM']:.3f} (p={res_age_ft['p_QM']:.4f})")

res_full_ft = run_meta_reg(t_ft, v_ft, X_full, "Full Multivariate FT")
r2_full_ft = max(0.0, (tau2_0_ft - res_full_ft['tau2']) / tau2_0_ft) * 100.0
print(f"Multivariate Model FT: tau2 = {res_full_ft['tau2']:.5f}, R2 = {r2_full_ft:.2f}%, QM = {res_full_ft['QM']:.3f} (p={res_full_ft['p_QM']:.4f})")
for name, b, s, p in zip(res_full_ft['coef_names'], res_full_ft['beta'], res_full_ft['se'], res_full_ft['pval']):
    print(f"  {name}: beta={b:.4f}, se={s:.4f}, p={p:.4f}")


print("\n=== MODALITY 4 CATEGORIES (WGS vs qPCR/TAC vs PCR vs Culture) ===")
df_eaec['modality4'] = df_eaec['modality_cat4']
X_mod4 = pd.get_dummies(df_eaec['modality4'], drop_first=True, dtype=float)
res_mod4 = run_meta_reg(y_logit, v_logit, X_mod4, "Modality 4 (Culture ref)")
r2_mod4 = max(0.0, (tau2_0 - res_mod4['tau2']) / tau2_0) * 100.0
print(f"Modality 4 alone (Logit): tau2 = {res_mod4['tau2']:.5f}, R2 = {r2_mod4:.2f}%, QM = {res_mod4['QM']:.3f} (p={res_mod4['p_QM']:.4f})")

res_mod4_ft = run_meta_reg(t_ft, v_ft, X_mod4, "Modality 4 FT")
r2_mod4_ft = max(0.0, (tau2_0_ft - res_mod4_ft['tau2']) / tau2_0_ft) * 100.0
print(f"Modality 4 alone (FT): tau2 = {res_mod4_ft['tau2']:.5f}, R2 = {r2_mod4_ft:.2f}%, QM = {res_mod4_ft['QM']:.3f} (p={res_mod4_ft['p_QM']:.4f})")

X_full4 = pd.concat([X_mod4, X_era, X_age], axis=1)
res_full4 = run_meta_reg(y_logit, v_logit, X_full4, "Full Multivariate 4 (Logit)")
r2_full4 = max(0.0, (tau2_0 - res_full4['tau2']) / tau2_0) * 100.0
print(f"Full Multivariate 4 (Logit): tau2 = {res_full4['tau2']:.5f}, R2 = {r2_full4:.2f}%, QM = {res_full4['QM']:.3f} (p={res_full4['p_QM']:.4f})")

res_full4_ft = run_meta_reg(t_ft, v_ft, X_full4, "Full Multivariate 4 (FT)")
r2_full4_ft = max(0.0, (tau2_0_ft - res_full4_ft['tau2']) / tau2_0_ft) * 100.0
print(f"Full Multivariate 4 (FT): tau2 = {res_full4_ft['tau2']:.5f}, R2 = {r2_full4_ft:.2f}%, QM = {res_full4_ft['QM']:.3f} (p={res_full4_ft['p_QM']:.4f})")
for name, b, s, p in zip(res_full4_ft['coef_names'], res_full4_ft['beta'], res_full4_ft['se'], res_full4_ft['pval']):
    print(f"  {name}: beta={b:.4f}, se={s:.4f}, p={p:.4f}")

