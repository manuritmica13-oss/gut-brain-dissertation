#!/usr/bin/env python3
"""Train classifiers on Wallen MetaPhlAn relative abundances.

Steps:
- load metaphlan_rel_ab.tsv (features x samples) and transpose
- load subject_metadata.tsv and filter to Case/Control
- align tables, encode Case=1, Control=0
- train/test split 70/30 stratified
- apply SMOTE to training set
- train RandomForest and XGBoost
- evaluate AUPRC (primary), AUROC, sensitivity, specificity, confusion matrix
- 200-iteration permutation null for AUPRC
- SHAP for top 20 features, save plot
- write metrics to results/ml_wallen_results.txt

Requires: scikit-learn, xgboost, shap, imblearn, pandas, numpy, matplotlib
"""

from __future__ import annotations

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except Exception:
    XGBOOST_AVAILABLE = False
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve, auc, roc_auc_score, confusion_matrix, recall_score
from sklearn.preprocessing import LabelEncoder

import shap

from imblearn.over_sampling import SMOTE

ROOT = Path(__file__).resolve().parents[2]
feat_fp = ROOT / 'data' / 'processed' / 'wallen' / 'metaphlan_rel_ab.tsv'
meta_fp = ROOT / 'data' / 'processed' / 'wallen' / 'subject_metadata.tsv'
OUT_DIR = ROOT / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)
METRICS_F = OUT_DIR / 'ml_wallen_results.txt'
SHAP_PNG = OUT_DIR / 'shap_wallen.png'

def load_data():
    if not feat_fp.exists():
        raise FileNotFoundError(feat_fp)
    if not meta_fp.exists():
        raise FileNotFoundError(meta_fp)
    feat = pd.read_csv(feat_fp, sep='\t', index_col=0).T
    meta = pd.read_csv(meta_fp, sep='\t')
    return feat, meta


def clean_taxon_name(name):
    if not isinstance(name, str):
        return str(name)
    last = name.split('|')[-1]
    for prefix in ('s__', 'g__'):
        if last.startswith(prefix):
            return last[len(prefix):]
    return last


def prepare(feat: pd.DataFrame, meta: pd.DataFrame):
    # filter metadata to known Case_status
    if 'Case_status' not in meta.columns:
        raise KeyError('Case_status column not found in metadata')
    if 'sample_name' not in meta.columns:
        raise KeyError('sample_name column not found in metadata')
    sel = meta['Case_status'].isin(['PD', 'Control'])
    meta_f = meta.loc[sel, ['sample_name', 'Case_status']].copy()
    meta_f = meta_f.drop_duplicates(subset='sample_name')

    # align on sample IDs: feat.index are sample IDs (e.g., DC001)
    common_ids = feat.index.intersection(meta_f['sample_name'])
    X = feat.loc[common_ids].copy()
    y = meta_f.set_index('sample_name').loc[common_ids, 'Case_status']
    y = (y == 'PD').astype(int)
    return X, y

def evaluate_model(clf, X_test, y_test):
    y_prob = clf.predict_proba(X_test)[:,1]
    y_pred = clf.predict(X_test)
    precision, recall, _ = precision_recall_curve(y_test, y_prob)
    auprc = auc(recall, precision)
    auroc = roc_auc_score(y_test, y_prob)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    return {'auprc': auprc, 'auroc': auroc, 'sensitivity': sensitivity, 'specificity': specificity, 'confusion': (tn, fp, fn, tp)}

def permutation_null(clf_factory, X_train, y_train, X_test, y_test, niter=200, random_state=0):
    rng = np.random.RandomState(random_state)
    scores = []
    for i in range(niter):
        y_perm = rng.permutation(y_train)
        clf = clf_factory()
        clf.fit(X_train, y_perm)
        y_prob = clf.predict_proba(X_test)[:,1]
        precision, recall, _ = precision_recall_curve(y_test, y_prob)
        scores.append(auc(recall, precision))
    return np.array(scores)

def main():
    feat, meta = load_data()
    X, y = prepare(feat, meta)

    # print class counts in aligned dataset before splitting/SMOTE
    n_pd = int((y == 1).sum())
    n_ctrl = int((y == 0).sum())
    print(f'Aligned samples: {len(y)} total ({n_pd} PD, {n_ctrl} Control)')

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)

    # SMOTE on training set
    sm = SMOTE(random_state=42)
    X_train_sm, y_train_sm = sm.fit_resample(X_train, y_train)

    # classifiers
    rf = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)
    rf.fit(X_train_sm, y_train_sm)
    rf_metrics = evaluate_model(rf, X_test, y_test)

    xgb_metrics = None
    if XGBOOST_AVAILABLE:
        xgb = XGBClassifier(n_estimators=500, use_label_encoder=False, eval_metric='logloss', random_state=42)
        xgb.fit(X_train_sm, y_train_sm)
        xgb_metrics = evaluate_model(xgb, X_test, y_test)
    else:
        print('XGBoost not available; skipping XGBoost model (install libomp/libxgboost).')

    # permutation null for RF
    def rf_factory():
        return RandomForestClassifier(n_estimators=200, random_state=0, n_jobs=1)
    null_scores = permutation_null(rf_factory, X_train_sm, y_train_sm, X_test, y_test, niter=200, random_state=42)

    # SHAP for top features (use TreeExplainer)
    # prefer XGBoost explainer if available, otherwise use RF
    model_for_shap = None
    if XGBOOST_AVAILABLE:
        model_for_shap = xgb
    else:
        model_for_shap = rf
    explainer = shap.TreeExplainer(model_for_shap)
    shap_vals = explainer.shap_values(X_test)
    # shap_vals shape depends on binary/multi; normalize to 2D array for positive class
    if isinstance(shap_vals, list):
        shap_arr = shap_vals[1]
    else:
        shap_arr = shap_vals
    if isinstance(shap_arr, np.ndarray) and shap_arr.ndim == 3:
        shap_arr = shap_arr[:, :, 1]
    shap_arr = np.asarray(shap_arr)
    mean_abs_shap = np.mean(np.abs(shap_arr), axis=0)
    feat_names = X_test.columns.tolist()
    clean_names = [clean_taxon_name(f) for f in feat_names]
    df_shap = pd.DataFrame({
        'feature_orig': feat_names,
        'feature_label': clean_names,
        'mean_abs_shap': mean_abs_shap,
    })
    df_shap = df_shap.sort_values('mean_abs_shap', ascending=False).head(20)

    # plot SHAP summary for top 20 with readable species/genus labels
    top_features = df_shap['feature_orig'].tolist()
    top_idx = [X_test.columns.get_loc(f) for f in top_features]
    shap_vals_top = shap_arr[:, top_idx]
    X_top = X_test[top_features]
    plt.figure(figsize=(18, 10))
    shap.summary_plot(
        shap_vals_top,
        X_top,
        feature_names=df_shap['feature_label'].tolist(),
        show=False,
    )
    plt.gcf().subplots_adjust(left=0.38)
    plt.savefig(SHAP_PNG, dpi=300)

    # save metrics
    with open(METRICS_F, 'w') as fh:
        fh.write('RandomForest metrics:\n')
        for k,v in rf_metrics.items():
            fh.write(f"{k}\t{v}\n")
        if xgb_metrics is not None:
            fh.write('\nXGBoost metrics:\n')
            for k,v in xgb_metrics.items():
                fh.write(f"{k}\t{v}\n")
        else:
            fh.write('\nXGBoost metrics:\n')
            fh.write('not available in this environment\n')
        fh.write('\nPermutation null (AUPRC) summary:\n')
        fh.write(f"mean\t{null_scores.mean()}\nstd\t{null_scores.std()}\n")
    print('Saved metrics to', METRICS_F)
    print('Saved SHAP plot to', SHAP_PNG)

if __name__ == '__main__':
    main()
