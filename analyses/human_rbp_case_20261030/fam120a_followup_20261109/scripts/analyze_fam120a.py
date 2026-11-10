"""Audit FAM120A regional differences against strand, interval and RNA expression.

All counts use persisted peaks, the original hg19 reference and public baseline
RNA-seq. Genes, rather than individual peaks, are the unit of the paired summary.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT.parent
DATA = ROOT / 'data'
REGIONS = ['Promoter', 'UTR5', 'CDS', 'UTR3', 'Intron', 'Intergenic']
CELLS = ['K562', 'HepG2']
ACCESSIONS = {'K562': 'ENCFF368MXJ', 'HepG2': 'ENCFF741FRL'}
EXPRESSION = {'K562': ['ENCFF286GLL', 'ENCFF986DBN'],
              'HepG2': ['ENCFF533XPJ', 'ENCFF321JIT']}
BATCH_CONTROLS = {'K562': ['ENCFF912ATD', 'ENCFF539FKS'],
                  'HepG2': ['ENCFF813RCD', 'ENCFF065EZU']}


def expression_table(ids, folder, prefix):
    series = []
    for cell, accs in ids.items():
        for a in accs:
            d = pd.read_csv(folder / (a + '.TPM.tsv.gz'), sep='\t')
            series.append(d.set_index('gene_id').TPM.rename(prefix + '_' + cell + '_' + a))
    expr = pd.concat(series, axis=1).fillna(0)
    for cell in CELLS:
        expr[prefix + '_' + cell] = expr[[prefix + '_' + cell + '_' + a for a in ids[cell]]].mean(axis=1)
    return expr


def summarize(peaks, mode, cohort):
    rows = []
    for cell, p in peaks.groupby('cell', sort=False):
        counts = p.region.value_counts().reindex(REGIONS, fill_value=0)
        for reg, n in counts.items():
            rows.append(dict(mode=mode, cohort=cohort, cell=cell, region=reg,
                             n_peaks=len(p), n_genes=p.gene_id.nunique(), count=int(n),
                             fraction=n / len(p) if len(p) else np.nan))
    return rows


def gene_table(peaks):
    tab = peaks.groupby(['gene_id', 'cell', 'region']).size().unstack('region', fill_value=0)
    tab = tab.reindex(columns=REGIONS, fill_value=0)
    tab['n_peaks'] = tab.sum(axis=1)
    for r in REGIONS:
        tab[r + '_fraction'] = tab[r] / tab.n_peaks
    return tab.reset_index()


def main():
    expr = expression_table(EXPRESSION, DATA, 'baseline')
    expr = expr.join(expression_table(BATCH_CONTROLS, CASE / 'data/expression', 'batch_control'), how='outer').fillna(0)
    expr['baseline_fold_ratio'] = np.maximum(expr.baseline_K562, expr.baseline_HepG2) / np.minimum(expr.baseline_K562, expr.baseline_HepG2).replace(0, np.nan)
    expr['batch_control_fold_ratio'] = np.maximum(expr.batch_control_K562, expr.batch_control_HepG2) / np.minimum(expr.batch_control_K562, expr.batch_control_HepG2).replace(0, np.nan)
    expr.index.name = 'gene_id'
    expr.to_csv(DATA / 'gene_expression.tsv.gz', sep='\t')
    symbols = pd.concat([pd.read_csv(CASE / 'data/events' / (a + '.SE.tsv.gz'), sep='\t', usecols=['gene', 'gene_symbol'])
                         for a in ['ENCFF259KEJ', 'ENCFF630EAW']]).drop_duplicates('gene').set_index('gene').gene_symbol
    modes = ['original_center_unstranded', 'center_sense', 'interval_sense', 'interval_unstranded']
    all_peaks = []
    for mode in modes:
        for cell, acc in ACCESSIONS.items():
            p = pd.read_csv(DATA / f'{acc}_{mode}.tsv.gz', sep='\t')
            assert len(p) == {'K562': 5683, 'HepG2': 8081}[cell]
            p = p.join(expr, on='gene_id')
            p['gene_symbol'] = p.gene_id.map(symbols)
            all_peaks.append(p)
    peaks = pd.concat(all_peaks, ignore_index=True)
    peaks.to_csv(DATA / 'annotated_peaks_with_expression.tsv.gz', sep='\t', index=False)
    rows = []
    cohorts = {}
    for min_tpm in [0.1, 1, 5]:
        common = (expr.baseline_K562 >= min_tpm) & (expr.baseline_HepG2 >= min_tpm)
        for fc in [1.2, 2, 5, np.inf]:
            name = f'baseline_tpm{min_tpm:g}_fc{fc:g}'
            cohorts[name] = set(expr.index[common & (expr.baseline_fold_ratio <= fc)])
    common = (expr.batch_control_K562 >= 1) & (expr.batch_control_HepG2 >= 1)
    cohorts['batch_control_tpm1_fcinf'] = set(expr.index[common])
    cohorts['batch_control_tpm1_fc2'] = set(expr.index[common & (expr.batch_control_fold_ratio <= 2)])
    for mode, group in peaks.groupby('mode', sort=False):
        rows += summarize(group, mode, 'all_peaks')
        for name, genes in cohorts.items():
            sub = group[group.gene_id.isin(genes)]
            rows += summarize(sub, mode, name)
            if name in ['baseline_tpm1_fcinf', 'baseline_tpm1_fc2']:
                rows += summarize(sub[sub.n_overlap_genes.eq(1)], mode, name + '_unique_gene')
    summaries = pd.DataFrame(rows)
    summaries.to_csv(DATA / 'regional_expression_sensitivity.tsv', sep='\t', index=False)
    sense = peaks[peaks['mode'].eq('center_sense')].copy()
    stable_genes = cohorts['baseline_tpm1_fc2']
    per_gene = gene_table(sense[sense.gene_id.isin(stable_genes)])
    per_gene = per_gene.join(expr[['baseline_K562', 'baseline_HepG2', 'baseline_fold_ratio']], on='gene_id')
    per_gene['gene_symbol'] = per_gene.gene_id.map(symbols)
    per_gene.to_csv(DATA / 'per_gene_regions.tsv', sep='\t', index=False)
    gene_counts = per_gene.pivot(index='gene_id', columns='cell', values='n_peaks').fillna(0)
    rng = np.random.default_rng(20261109)
    paired_results = []
    for min_peaks in [1, 3, 5, 10]:
        ids = gene_counts.index[(gene_counts.K562 >= min_peaks) & (gene_counts.HepG2 >= min_peaks)]
        sub = per_gene[per_gene.gene_id.isin(ids)]
        subset_peaks = sense[sense.gene_id.isin(ids)]
        rows_shared = summarize(subset_peaks, 'center_sense', f'stable_common_bound_min{min_peaks}')
        summaries = pd.concat([summaries, pd.DataFrame(rows_shared)], ignore_index=True)
        for region in ['UTR3', 'Intron']:
            matrix = sub.pivot(index='gene_id', columns='cell', values=region + '_fraction')
            diff = (matrix.HepG2 - matrix.K562).to_numpy()
            if len(diff):
                boot = np.array([rng.choice(diff, len(diff), replace=True).mean() for _ in range(10000)])
                result = dict(min_peaks_per_cell=min_peaks, region=region, n_genes=len(diff),
                              mean_K562=float(matrix.K562.mean()), mean_HepG2=float(matrix.HepG2.mean()),
                              mean_difference_HepG2_minus_K562=float(diff.mean()),
                              bootstrap95_low=float(np.quantile(boot, .025)),
                              bootstrap95_high=float(np.quantile(boot, .975)),
                              wilcoxon_p=float(wilcoxon(diff).pvalue) if (diff != 0).any() else 1.)
                paired_results.append(result)
    summaries.to_csv(DATA / 'regional_expression_sensitivity.tsv', sep='\t', index=False)
    pd.DataFrame(paired_results).to_csv(DATA / 'paired_gene_comparison.tsv', sep='\t', index=False)
    # Predefine locus selection by expression comparability, gene identity and
    # support in both cell lines. Rank by the smaller peak count, not a P value.
    pivot = per_gene.pivot(index='gene_id', columns='cell', values=['n_peaks', 'UTR3_fraction', 'Intron_fraction'])
    pivot.columns = ['_'.join(c) for c in pivot.columns]
    pivot = pivot.dropna()
    pivot = pivot.join(expr[['baseline_K562', 'baseline_HepG2', 'baseline_fold_ratio']])
    pivot['gene_symbol'] = pivot.index.map(symbols)
    pivot['min_peak_support'] = pivot[['n_peaks_K562', 'n_peaks_HepG2']].min(axis=1)
    pivot['regional_difference'] = pivot.Intron_fraction_HepG2 - pivot.Intron_fraction_K562
    candidates = pivot[(pivot.min_peak_support >= 3) & (pivot.UTR3_fraction_K562 >= .5) &
                       (pivot.Intron_fraction_HepG2 >= .5)].sort_values(['min_peak_support', 'regional_difference'], ascending=False)
    candidates.to_csv(DATA / 'candidate_loci.tsv', sep='\t')
    report = dict(input_peak_counts={'K562': 5683, 'HepG2': 8081},
                  original_count_reproduction=True,
                  common_expressed_genes=len(cohorts['baseline_tpm1_fcinf']),
                  similarly_expressed_genes=len(stable_genes),
                  paired_gene_comparisons=paired_results,
                  candidate_loci=len(candidates),
                  top_candidates=candidates.head(12).reset_index().to_dict('records'))
    (DATA / 'analysis_report.json').write_text(json.dumps(report, indent=2))
    print(summaries[summaries['mode'].eq('center_sense') & summaries.cohort.isin(
        ['all_peaks', 'baseline_tpm1_fcinf', 'baseline_tpm1_fc2', 'stable_common_bound_min1'])].pivot_table(
        index=['cohort', 'cell', 'n_peaks'], columns='region', values='fraction').round(4).to_string())
    print('PAIRED', json.dumps(paired_results, indent=2))
    print('TOP_CANDIDATES', candidates.head(10).to_string())


if __name__ == '__main__':
    main()
