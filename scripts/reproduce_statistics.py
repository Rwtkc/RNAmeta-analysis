"""Recompute manuscript RBP statistics from the included annotation outputs."""
from pathlib import Path
import argparse, importlib.util, json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'analyses/human_rbp_case_20261030'

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/statistics')
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    spec=importlib.util.spec_from_file_location('original_analysis',CASE/'revision_20261030/analyze_revision.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.OUT=out
    frame=pd.read_csv(CASE/'data/landscape/datasets.tsv',sep='\t')
    a=frame[frame.cell.eq('K562')].set_index('protein');b=frame[frame.cell.eq('HepG2')].set_index('protein')
    ids=sorted(a.index.intersection(b.index))
    regional=cdist(a.loc[ids,module.REGIONS],b.loc[ids,module.REGIONS],metric='jensenshannon')
    report={'regional_cross_cell':module.summary_distances(regional,ids,'regional')}
    profiles=np.load(CASE/'summary/profiles.npz')
    lookup={k:{str(acc):module.profile_mass(profiles[k+'__x'],y,k) for acc,y in zip(profiles[k+'__ids'],profiles[k])} for k in module.PROFILE_KEYS}
    complete=[p for p in ids if all(str(cell.loc[p,'accession']) in lookup[k] for k in module.PROFILE_KEYS for cell in [a,b])]
    vectors=[np.stack([np.concatenate([np.sqrt(cell.loc[p,module.REGIONS].to_numpy(float)),
               *[np.sqrt(lookup[k][str(cell.loc[p,'accession'])]) for k in module.PROFILE_KEYS]])/2 for p in complete]) for cell in [a,b]]
    report['multiview_cross_cell']=module.summary_distances(cdist(*vectors)/np.sqrt(2),complete,'multiview')
    tests=module.functional_tests(frame,'all')
    assert len(tests)==114
    reference=json.loads((CASE/'revision_20261030/data/report.json').read_text())
    for view in ['regional_cross_cell','multiview_cross_cell']:
        for key in ['n','same_median','other_median','permutation_p','closer_than_other_median']:
            assert np.isclose(report[view][key],reference[view][key]),(view,key)
    expected=pd.read_csv(CASE/'revision_20261030/data/all_function_tests.tsv',sep='\t')
    np.testing.assert_allclose(tests[['p','q','cliffs_delta']],expected[['p','q','cliffs_delta']],atol=1e-12)
    genes=pd.read_csv(CASE/'fam120a_followup_20261109/data/per_gene_regions.tsv',sep='\t')
    paired=pd.read_csv(CASE/'fam120a_followup_20261109/data/paired_gene_comparison.tsv',sep='\t')
    counts=genes.pivot(index='gene_id',columns='cell',values='n_peaks').fillna(0)
    ids=counts.index[(counts.K562>=1)&(counts.HepG2>=1)]
    assert len(ids)==354
    rng=np.random.default_rng(20261109)
    report['fam120a']={}
    for region in ['UTR3','Intron']:
        matrix=genes[genes.gene_id.isin(ids)].pivot(index='gene_id',columns='cell',values=region+'_fraction')
        diff=(matrix.HepG2-matrix.K562).to_numpy()
        boot=np.array([rng.choice(diff,len(diff),replace=True).mean() for _ in range(10000)])
        ci=np.quantile(boot,[.025,.975]);saved=paired[(paired.min_peaks_per_cell==1)&(paired.region==region)].iloc[0]
        np.testing.assert_allclose([diff.mean(),*ci],[saved.mean_difference_HepG2_minus_K562,saved.bootstrap95_low,saved.bootstrap95_high],atol=1e-12)
        report['fam120a'][region]={'n_genes':len(diff),'mean_difference':float(diff.mean()),'bootstrap95':ci.tolist()}
    report['verification']='Matches the archived manuscript results'
    (out/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
