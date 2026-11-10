"""Post hoc threshold sensitivity and exact regional peak subsampling.

This extension reads the frozen outputs; it never reruns annotation or changes
the original results. See peak_count_sensitivity_plan_20261106.md.
"""
from pathlib import Path
import gzip
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import false_discovery_control, mannwhitneyu

from analyze_revision import REGIONS, PROFILE_KEYS, profile_mass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / 'data' / 'peak_count_sensitivity_20261106'
SEED = 20261106
CUTOFFS = [0, 100, 250, 500, 1000, 2000, 5000]
REPEATS = 200
PERMUTATIONS = 10000
BOOTSTRAPS = 10000


def paired_values(matrix):
    n = len(matrix)
    if n < 2 or matrix.shape != (n, n) or not np.isfinite(matrix).all():
        raise ValueError('A finite square matrix with at least two RBPs is required')
    off = matrix.copy()
    np.fill_diagonal(off, np.nan)
    same = np.diag(matrix)
    other = (np.nanmedian(off, axis=0) + np.nanmedian(off, axis=1)) / 2
    return same, other, other - same


def permutation_indices(n, rng):
    # Each row is an independent uniform permutation of protein identities.
    return np.stack([rng.permutation(n) for _ in range(PERMUTATIONS)])


def describe_matrix(matrix, permutations):
    same, other, delta = paired_values(matrix)
    null = matrix[np.arange(len(matrix))[None, :], permutations].mean(axis=1)
    return dict(n_rbps=len(matrix), same_median=float(np.median(same)),
                other_median=float(np.median(other)),
                median_difference=float(np.median(delta)),
                n_closer=int((delta > 0).sum()),
                fraction_closer=float((delta > 0).mean()),
                permutation_p=float((1 + (null <= same.mean()).sum()) /
                                    (1 + len(null))))


def threshold_result(matrix, ids, cutoff, view, rng):
    result = describe_matrix(matrix, permutation_indices(len(ids), rng))
    same, other, delta = paired_values(matrix)
    boot = np.median(delta[rng.integers(len(ids), size=(BOOTSTRAPS, len(ids)))], axis=1)
    result.update(cutoff=cutoff, view=view,
                  difference_ci_low=float(np.quantile(boot, .025)),
                  difference_ci_high=float(np.quantile(boot, .975)))
    pairs = pd.DataFrame(dict(cutoff=cutoff, view=view, protein=ids,
                              same_distance=same, other_median_distance=other,
                              difference=delta))
    return result, pairs


def function_tests(frame, cutoff):
    cols = [c for c in frame if c.startswith('function__')]
    means = frame.groupby('protein')[REGIONS].mean()
    labels = frame.groupby('protein')[cols].max()
    rows = []
    for col in cols:
        flag = labels[col].astype(bool)
        if flag.sum() < 5 or (~flag).sum() < 5:
            continue
        for region in REGIONS:
            a = means.loc[flag, region]
            b = means.loc[~flag, region]
            test = mannwhitneyu(a, b, alternative='two-sided')
            rows.append(dict(cutoff=cutoff, function=col.removeprefix('function__'),
                             region=region, n_positive=len(a), n_other=len(b),
                             positive_median=float(a.median()), other_median=float(b.median()),
                             cliffs_delta=float(2 * test.statistic / (len(a) * len(b)) - 1),
                             p=float(test.pvalue)))
    table = pd.DataFrame(rows)
    if len(table):
        table['q_within_cutoff'] = false_discovery_control(table.p)
    return table


def exact_samples(counts, depth, repeats, rng):
    if depth < 1 or np.any(counts < 0) or np.any(counts.sum(axis=1) < depth):
        raise ValueError('Subsampling depth must be positive and no larger than each input')
    draws = np.stack([rng.multivariate_hypergeometric(row, depth, size=repeats)
                      for row in counts], axis=1)
    assert np.all(draws.sum(axis=2) == depth)
    assert np.all(draws <= counts[None, :, :])
    return draws


def validate_sampler():
    # Edge cases of exact finite-population sampling, including zero categories.
    counts = np.array([[0, 2, 3, 0, 1, 0], [1, 0, 0, 5, 0, 0]])
    rng = np.random.default_rng(SEED)
    assert np.array_equal(exact_samples(counts, 6, 3, rng),
                          np.broadcast_to(counts, (3, *counts.shape)))
    sample = exact_samples(counts, 1, 12, rng)
    assert np.all((sample == 0) | (sample == 1))
    try:
        exact_samples(counts, 7, 1, rng)
    except ValueError:
        pass
    else:
        raise AssertionError('Invalid depth was accepted')
    # Check expectation and finite-population variance against exact moments.
    row = np.array([[2, 3, 5, 0, 0, 0]])
    draws = exact_samples(row, 4, 30000, rng)[:, 0, :]
    expected = row[0] / 10 * 4
    variance = 4 * row[0] / 10 * (1 - row[0] / 10) * (10 - 4) / (10 - 1)
    assert np.allclose(draws.mean(axis=0), expected, atol=.025)
    assert np.allclose(draws.var(axis=0), variance, atol=.025)


def main():
    OUT.mkdir(exist_ok=True)
    validate_sampler()
    data = pd.read_csv(ROOT / 'data/landscape/datasets.tsv', sep='\t')
    with gzip.open(ROOT / 'summary/annotation.json.gz', 'rt', encoding='utf-8') as stream:
        annotation = json.load(stream)
    count_map = {r['accession']:np.array([r['regions'].get(k, 0) for k in REGIONS], dtype=np.int64)
                 for r in annotation}
    for row in data.itertuples():
        counts = count_map[row.accession]
        assert counts.sum() == row.n_filtered
        assert np.allclose(counts / counts.sum(), data.set_index('accession').loc[row.accession, REGIONS].to_numpy(float))
    a = data[data.cell.eq('K562')].set_index('protein')
    b = data[data.cell.eq('HepG2')].set_index('protein')
    ids = sorted(a.index.intersection(b.index))
    peak_min = np.array([min(a.loc[p, 'n_filtered'], b.loc[p, 'n_filtered']) for p in ids])
    full_regions = cdist(a.loc[ids, REGIONS], b.loc[ids, REGIONS], metric='jensenshannon')
    saved = pd.read_csv(HERE / 'data/regional_distances.tsv', sep='\t', index_col=0)
    assert saved.index.tolist() == ids and np.allclose(saved, full_regions)
    profiles = np.load(ROOT / 'summary/profiles.npz')
    lookup = {k:{str(acc):profile_mass(profiles[k + '__x'], y, k)
                 for acc, y in zip(profiles[k + '__ids'], profiles[k])} for k in PROFILE_KEYS}
    complete = [p for p in ids if all(str(cell.loc[p, 'accession']) in lookup[k]
                                     for k in PROFILE_KEYS for cell in [a, b])]
    vectors = [np.stack([np.concatenate([np.sqrt(frame.loc[p, REGIONS].to_numpy(float)),
                *[np.sqrt(lookup[k][str(frame.loc[p, 'accession'])]) for k in PROFILE_KEYS]]) / 2
                for p in complete]) for frame in [a, b]]
    full_multi = cdist(*vectors) / np.sqrt(2)
    saved_multi = pd.read_csv(HERE / 'data/multiview_distances.tsv', sep='\t', index_col=0)
    assert saved_multi.index.tolist() == complete and np.allclose(saved_multi, full_multi)
    results, pairs, function_tables, cohorts, function_counts = [], [], [], [], []
    rng = np.random.default_rng(SEED)
    for cutoff in CUTOFFS:
        for view, names, matrix in [('Regional', ids, full_regions), ('Four views', complete, full_multi)]:
            keep = [p for p in names if min(a.loc[p, 'n_filtered'], b.loc[p, 'n_filtered']) >= cutoff]
            ix = [names.index(p) for p in keep]
            record, table = threshold_result(matrix[np.ix_(ix, ix)], keep, cutoff, view, rng)
            results.append(record)
            pairs.append(table)
            cohorts.extend(dict(view=view, cutoff=cutoff, protein=p,
                                K562_peaks=int(a.loc[p, 'n_filtered']),
                                HepG2_peaks=int(b.loc[p, 'n_filtered'])) for p in keep)
        selected = data[data.n_filtered >= cutoff]
        tests = function_tests(selected, cutoff)
        function_tables.append(tests)
        function_counts.append(dict(cutoff=cutoff, n_datasets=len(selected),
                                    n_rbps=selected.protein.nunique(), n_tests=len(tests),
                                    n_q_below_05=int((tests.q_within_cutoff < .05).sum())))
    thresholds = pd.DataFrame(results)
    thresholds['permutation_q_grid'] = false_discovery_control(thresholds.permutation_p)
    thresholds.to_csv(OUT / 'threshold_summary.tsv', sep='\t', index=False)
    pd.concat(pairs, ignore_index=True).to_csv(OUT / 'threshold_pairs.tsv', sep='\t', index=False)
    pd.DataFrame(cohorts).to_csv(OUT / 'threshold_cohorts.tsv', sep='\t', index=False)
    functions = pd.concat(function_tables, ignore_index=True)
    functions.to_csv(OUT / 'function_tests_all_thresholds.tsv', sep='\t', index=False)
    pd.DataFrame(function_counts).to_csv(OUT / 'function_cohorts.tsv', sep='\t', index=False)
    primary_functions = functions[((functions.function == 'Splicing regulation') & (functions.region == 'Intron')) |
                                 ((functions.function == 'Translation regulation') & (functions.region == 'UTR3'))]
    primary_functions.to_csv(OUT / 'highlighted_function_tests.tsv', sep='\t', index=False)
    depths = [int(peak_min.min()), *CUTOFFS[1:]]
    repeat_rows, subsample_summaries, protein_rows = [], [], []
    sampled_counts = {}
    for depth in depths:
        keep = [p for p, m in zip(ids, peak_min) if m >= depth]
        ix = [ids.index(p) for p in keep]
        permutations = permutation_indices(len(keep), rng)
        baseline = describe_matrix(full_regions[np.ix_(ix, ix)], permutations)
        sampled = []
        for frame in [a, b]:
            counts = np.stack([count_map[frame.loc[p, 'accession']] for p in keep])
            sampled.append(exact_samples(counts, depth, REPEATS, rng))
        sampled_counts[f'depth_{depth}_proteins'] = np.array(keep)
        sampled_counts[f'depth_{depth}_K562'] = sampled[0]
        sampled_counts[f'depth_{depth}_HepG2'] = sampled[1]
        metrics, differences = [], []
        for repeat in range(REPEATS):
            matrix = cdist(sampled[0][repeat] / depth, sampled[1][repeat] / depth,
                           metric='jensenshannon')
            record = describe_matrix(matrix, permutations)
            record.update(depth=depth, repeat=repeat + 1)
            repeat_rows.append(record)
            metrics.append(record)
            differences.append(paired_values(matrix)[2])
        summary = dict(depth=depth, n_rbps=len(keep), repeats=REPEATS)
        for key in ['same_median', 'other_median', 'median_difference', 'fraction_closer', 'permutation_p']:
            values = np.array([r[key] for r in metrics])
            summary[key] = float(np.median(values))
            summary[key + '_p025'] = float(np.quantile(values, .025))
            summary[key + '_p975'] = float(np.quantile(values, .975))
            summary['full_cohort_' + key] = baseline[key]
        subsample_summaries.append(summary)
        differences = np.stack(differences)
        full_differences = paired_values(full_regions[np.ix_(ix, ix)])[2]
        for i, protein in enumerate(keep):
            values = differences[:, i]
            protein_rows.append(dict(depth=depth, protein=protein, full_difference=full_differences[i],
                                     median_difference=float(np.median(values)),
                                     difference_p025=float(np.quantile(values, .025)),
                                     difference_p975=float(np.quantile(values, .975)),
                                     fraction_repeats_closer=float((values > 0).mean())))
        print(json.dumps(summary), flush=True)
    pd.DataFrame(repeat_rows).to_csv(OUT / 'subsampling_repetitions.tsv', sep='\t', index=False)
    pd.DataFrame(subsample_summaries).to_csv(OUT / 'subsampling_summary.tsv', sep='\t', index=False)
    pd.DataFrame(protein_rows).to_csv(OUT / 'subsampling_per_protein.tsv', sep='\t', index=False)
    np.savez_compressed(OUT / 'subsampled_region_counts.npz', **sampled_counts)
    paths = [ROOT / 'data/landscape/datasets.tsv', ROOT / 'summary/annotation.json.gz',
             ROOT / 'summary/profiles.npz', HERE / 'data/regional_distances.tsv',
             HERE / 'data/multiview_distances.tsv', Path(__file__).resolve(),
             HERE / 'peak_count_sensitivity_plan_20261106.md']
    report = dict(seed=SEED, cutoffs=CUTOFFS, depths=depths, repeats=REPEATS,
                  permutations=PERMUTATIONS, bootstraps=BOOTSTRAPS,
                  sampling='Exact multivariate hypergeometric, without replacement; regional composition only',
                  interval='RBP bootstrap CI for threshold effects; subsampling percentile range for rarefaction',
                  input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  threshold_results=thresholds.to_dict('records'),
                  highlighted_functions=primary_functions.to_dict('records'),
                  subsampling_results=subsample_summaries)
    (OUT / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('COMPLETE', OUT, flush=True)


if __name__ == '__main__':
    main()
