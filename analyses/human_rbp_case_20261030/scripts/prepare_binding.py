"""Freeze official hg19 IDR inputs and their provenance for Online annotation."""
from pathlib import Path
from collections import Counter
import gzip
import hashlib
import json
import tarfile
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REMOTE = '/share/shiny/human_rbp_case_20261030'
ALGORITHM = 'e813511068942469272dabd321adb15c4e2f1cdef658f46ec3dbdb5b97e18948'


def main():
    manifest = json.loads((ROOT/'download_manifest.json').read_text())
    inputs = ROOT/'data/binding'; inputs.mkdir(parents=True, exist_ok=True)
    records = []; audit = []
    for r in manifest:
        if r['kind'] != 'binding': continue
        raw = ROOT/r['path']
        assert hashlib.md5(raw.read_bytes()).hexdigest() == r['md5sum']
        data = pd.read_csv(raw, sep='\t', header=None)
        assert data.shape[1] == 10 and (data[2] > data[1]).all()
        # ENCODE eCLIP IDR narrowPeak stores log2 fold enrichment and -log10 P.
        keep = data[(data[6] >= 3) & (data[7] >= 3)]
        assert len(keep) > 0
        path = inputs/(r['accession']+'.bed')
        keep.iloc[:, :6].to_csv(path, sep='\t', header=False, index=False)
        content = path.read_bytes()
        record = dict(r, job_id=r['accession'], species='hg19',
                      bed_path=REMOTE+'/data/binding/'+path.name,
                      input_md5=hashlib.md5(content).hexdigest(),
                      input_sha256=hashlib.sha256(content).hexdigest(),
                      algorithm_id=ALGORITHM, seed=20261030, n_raw=len(data), n_filtered=len(keep))
        records.append(record)
        audit.append({k:record[k] for k in ['accession','experiment','protein','cell','n_raw','n_filtered','input_sha256']})
    (ROOT/'binding_manifest.json').write_text(json.dumps(records, indent=2))
    pd.DataFrame(audit).to_csv(ROOT/'data/binding_audit.tsv', sep='\t', index=False)
    report = dict(files=len(records), proteins=len({r['protein'] for r in records}),
                  cells=dict(Counter(r['cell'] for r in records)),
                  raw_peaks=sum(r['n_raw'] for r in records),
                  filtered_peaks=sum(r['n_filtered'] for r in records),
                  filters={'minimum_log2_fold':3,'minimum_negative_log10_p':3})
    (ROOT/'data/binding_report.json').write_text(json.dumps(report, indent=2))
    with tarfile.open(ROOT/'binding_payload.tar.gz','w:gz') as tar:
        tar.add(inputs, arcname='data/binding')
        tar.add(ROOT/'binding_manifest.json',arcname='binding_manifest.json')
    print(json.dumps(report))


if __name__ == '__main__': main()
