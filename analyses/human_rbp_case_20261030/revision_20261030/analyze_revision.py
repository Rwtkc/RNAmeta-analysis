"""Fixed, label-free positional comparisons and complete functional tests."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon, false_discovery_control, rankdata

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=HERE/'data'
OUT.mkdir(exist_ok=True)
REGIONS=['Promoter','UTR5','CDS','UTR3','Intron','Intergenic']
PROFILE_KEYS=['metagene_mrna','splicesite__5PSS','splicesite__3PSS']
SEED=20261030


def summary_distances(matrix,ids,label):
    n=len(ids)
    same=np.diag(matrix)
    median_row=np.array([np.median(np.delete(matrix[i],i)) for i in range(n)])
    median_col=np.array([np.median(np.delete(matrix[:,i],i)) for i in range(n)])
    other=(median_row+median_col)/2
    # A protein, rather than each off-diagonal comparison, is the paired unit.
    delta=other-same
    rng=np.random.default_rng(SEED)
    null=np.array([matrix[np.arange(n),rng.permutation(n)].mean() for _ in range(10000)])
    boot=np.median(delta[rng.integers(0,n,size=(10000,n))],axis=1)
    rank_forward=np.array([(matrix[i]<=matrix[i,i]).sum() for i in range(n)])
    rank_reverse=np.array([(matrix[:,i]<=matrix[i,i]).sum() for i in range(n)])
    table=pd.DataFrame(dict(protein=ids,same_distance=same,other_median_distance=other,difference=delta,
        K562_to_HepG2_rank=rank_forward,HepG2_to_K562_rank=rank_reverse))
    table.to_csv(OUT/(label+'_pairs.tsv'),sep='\t',index=False)
    pd.DataFrame(matrix,index=ids,columns=ids).to_csv(OUT/(label+'_distances.tsv'),sep='\t')
    return dict(n=n,same_median=float(np.median(same)),other_median=float(np.median(other)),
        median_difference=float(np.median(delta)),median_difference_ci=np.quantile(boot,[.025,.975]).tolist(),
        paired_p=float(wilcoxon(delta,alternative='two-sided').pvalue),
        permutation_p=float((1+np.sum(null<=same.mean()))/(1+len(null))),
        mean_true_distance=float(same.mean()),mean_permuted_distance=float(null.mean()),
        closer_than_other_median=int((delta>0).sum()),
        forward_top1=int((rank_forward==1).sum()),reverse_top1=int((rank_reverse==1).sum()),
        forward_top5=int((rank_forward<=5).sum()),reverse_top5=int((rank_reverse<=5).sum()),
        forward_median_rank=float(np.median(rank_forward)),reverse_median_rank=float(np.median(rank_reverse)))


def profile_mass(x,y,key):
    if key.startswith('metagene'):
        weights=np.ones(len(x))
    else:
        weights=np.empty(len(x));weights[0]=(x[1]-x[0])/2;weights[-1]=(x[-1]-x[-2])/2
        weights[1:-1]=(x[2:]-x[:-2])/2
    a=y*weights
    return a/a.sum()


def functional_tests(d,label):
    fcols=[c for c in d if c.startswith('function__')]
    means=d.groupby('protein')[REGIONS].mean()
    labels=d.groupby('protein')[fcols].max()
    rows=[]
    for c in fcols:
        flag=labels[c].astype(bool)
        if flag.sum()<5 or (~flag).sum()<5:continue
        for region in REGIONS:
            a=means.loc[flag,region];b=means.loc[~flag,region]
            u=mannwhitneyu(a,b,alternative='two-sided')
            rows.append(dict(function=c.removeprefix('function__'),region=region,n_positive=len(a),n_other=len(b),
                positive_median=float(a.median()),other_median=float(b.median()),
                cliffs_delta=float(2*u.statistic/(len(a)*len(b))-1),p=float(u.pvalue)))
    t=pd.DataFrame(rows);t['q']=false_discovery_control(t.p)
    t.to_csv(OUT/(label+'_function_tests.tsv'),sep='\t',index=False)
    return t


def response_analysis(d):
    rows=[];cache={}
    for e in json.loads((ROOT/'experiment_design.json').read_text()):
        p=ROOT/'data/events'/(e['accession']+'.SE.tsv.gz')
        if not p.exists():continue
        f=pd.read_csv(p,sep='\t');b=d[(d.protein==e['protein'])&(d.cell==e['cell'])].iloc[0]
        for a in e['control_rsem'].split(','):
            if a not in cache:cache[a]=pd.read_csv(ROOT/'data/expression'/(a+'.TPM.tsv.gz'),sep='\t').set_index('gene_id').TPM
        expr=np.mean([cache[a].get(b.gene,np.nan) for a in e['control_rsem'].split(',')])
        high=f[f.min_reads>=30]
        rows.append(dict(protein=e['protein'],cell=e['cell'],accession=e['accession'],intron=b.Intron,
            n_covered=len(f),n_peaks=b.n_filtered,response=f.regulated.mean(),n_regulated=int(f.regulated.sum()),
            response_reads30=high.regulated.mean(),n_covered_reads30=len(high),control_rbp_tpm=expr))
    t=pd.DataFrame(rows);t.to_csv(OUT/'response_experiments.tsv',sep='\t',index=False)
    m=t.groupby('protein').mean(numeric_only=True);m.to_csv(OUT/'response_proteins.tsv',sep='\t')
    tests=[]
    for key in ['response','response_reads30']:
        s=spearmanr(m.intron,m[key],nan_policy='omit')
        tests.append(dict(outcome=key,n=len(m),rho=float(s.statistic),p=float(s.pvalue)))
    good=m.dropna(subset=['intron','response','control_rbp_tpm'])
    nuisance=np.column_stack([np.ones(len(good)),*[rankdata(np.log1p(good[c])) for c in ['n_peaks','n_covered','control_rbp_tpm']]])
    def residual(v):return v-nuisance@np.linalg.lstsq(nuisance,v,rcond=None)[0]
    x=residual(rankdata(good.intron));y=residual(rankdata(good.response))
    rho=float(np.corrcoef(x,y)[0,1]);rng=np.random.default_rng(SEED)
    null=np.array([np.corrcoef(x,rng.permutation(y))[0,1] for _ in range(10000)])
    tests.append(dict(outcome='partial_rank_adjusted',n=len(good),rho=rho,p=float((1+(np.abs(null)>=abs(rho)).sum())/10001)))
    return tests


def main():
    d=pd.read_csv(ROOT/'data/landscape/datasets.tsv',sep='\t')
    a=d[d.cell.eq('K562')].set_index('protein');b=d[d.cell.eq('HepG2')].set_index('protein')
    ids=sorted(a.index.intersection(b.index))
    region_matrix=cdist(a.loc[ids,REGIONS],b.loc[ids,REGIONS],metric='jensenshannon')
    report={'regional_cross_cell':summary_distances(region_matrix,ids,'regional')}
    for cutoff in [500]:
        keep=[p for p in ids if min(a.loc[p,'n_filtered'],b.loc[p,'n_filtered'])>=cutoff]
        ix=[ids.index(p) for p in keep]
        report['regional_min500']=summary_distances(region_matrix[np.ix_(ix,ix)],keep,'regional_min500')
    profiles=np.load(ROOT/'summary/profiles.npz')
    profile_lookup={key:{str(acc):profile_mass(profiles[key+'__x'],y,key) for acc,y in zip(profiles[key+'__ids'],profiles[key])} for key in PROFILE_KEYS}
    complete=[p for p in ids if all(str(cell.loc[p,'accession']) in profile_lookup[k] for k in PROFILE_KEYS for cell in [a,b])]
    vectors=[]
    for frame in [a,b]:
        vectors.append(np.stack([np.concatenate([np.sqrt(frame.loc[p,REGIONS].to_numpy(float)),*[np.sqrt(profile_lookup[k][str(frame.loc[p,'accession'])]) for k in PROFILE_KEYS]])/2 for p in complete]))
    multiview=cdist(*vectors)/np.sqrt(2)
    regions_hellinger=cdist(np.sqrt(a.loc[complete,REGIONS].to_numpy(float)),np.sqrt(b.loc[complete,REGIONS].to_numpy(float)))/np.sqrt(2)
    report['multiview_cross_cell']=summary_distances(multiview,complete,'multiview')
    report['region_same_cohort']=summary_distances(regions_hellinger,complete,'region_same_cohort')
    keep=[p for p in complete if min(a.loc[p,'n_filtered'],b.loc[p,'n_filtered'])>=500]
    ix=[complete.index(p) for p in keep]
    report['multiview_min500']=summary_distances(multiview[np.ix_(ix,ix)],keep,'multiview_min500')
    forward=(OUT/'multiview_pairs.tsv');m=pd.read_csv(forward,sep='\t');r=pd.read_csv(OUT/'region_same_cohort_pairs.tsv',sep='\t')
    diff=((r.K562_to_HepG2_rank+r.HepG2_to_K562_rank)-(m.K562_to_HepG2_rank+m.HepG2_to_K562_rank))/2
    report['multiview_rank_comparison']=dict(n=len(diff),median_improvement=float(diff.median()),paired_p=float(wilcoxon(diff).pvalue))
    full=functional_tests(d,'all');filtered=functional_tests(d[d.n_filtered>=500],'min500')
    report['function_tests']=dict(total=len(full),fdr05=int((full.q<.05).sum()),min500_total=len(filtered),min500_fdr05=int((filtered.q<.05).sum()))
    report['response_tests']=response_analysis(d)
    p=[report[k]['paired_p'] for k in ['regional_cross_cell','multiview_cross_cell']]+[r['p'] for r in report['response_tests']]
    report['five_primary_secondary_tests_bh']=false_discovery_control(p).tolist()
    (OUT/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2));print(full[full.q<.05].sort_values('q').to_string(index=False))


if __name__=='__main__':main()
