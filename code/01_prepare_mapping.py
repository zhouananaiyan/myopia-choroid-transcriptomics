"""Prepare versioned gene mappings; deduplicate independently of significance."""
import argparse
import hashlib
import json
import shutil
import time
import zlib
from pathlib import Path

import pandas as pd
import requests


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fetch(url, path):
    if not path.exists():
        print('Downloading', url, flush=True)
        tmp = path.with_suffix(path.suffix + '.part')
        for attempt in range(6):
            start = tmp.stat().st_size if tmp.exists() else 0
            try:
                with requests.get(url, headers={'Range':f'bytes={start}-'}, stream=True,
                                  timeout=(30,180)) as response:
                    response.raise_for_status()
                    resumed = start > 0 and response.status_code == 206
                    total = int(response.headers.get('Content-Length',0)) + (start if resumed else 0)
                    with tmp.open('ab' if resumed else 'wb') as handle:
                        for chunk in response.iter_content(1024*1024):
                            handle.write(chunk)
                    if total and tmp.stat().st_size != total:
                        raise IOError('Incomplete download')
                tmp.replace(path)
                break
            except requests.RequestException as exc:
                print('Download retry',attempt+1,str(exc)[:130],flush=True)
                time.sleep(2)
        if not path.exists():
            raise RuntimeError(f'Download incomplete: {url}')
    return {'url': url, 'file': path.name, 'bytes': path.stat().st_size,
            'sha256': sha256(path)}


def features(path, species):
    x = pd.read_csv(path, sep='\t', header=None, names=['gene_id', 'symbol', 'type'])
    # Seurat make.unique uses .1, .2, ... for repeated feature names.
    seen = {}
    names = []
    used = set(x.symbol.astype(str))
    for value in x.symbol.astype(str):
        i = seen.get(value, 0)
        name = value if i == 0 else f'{value}.{i}'
        while i > 0 and name in used:
            i += 1
            name = f'{value}.{i}'
        seen[value] = i + 1
        names.append(name)
        used.add(name)
    x['analysis_feature'] = names
    x['species'] = species
    x['gene_id'] = x.gene_id.str.split('.').str[0]
    return x


def copy_valid_gzip(source, destination):
    """Copy the first complete gzip member and reject extra gzip members."""
    data = source.read_bytes()
    decoder = zlib.decompressobj(wbits=31)
    raw = decoder.decompress(data) + decoder.flush()
    if not raw:
        raise ValueError(f'Empty gzip payload: {source}')
    used = len(data) - len(decoder.unused_data)
    if decoder.unused_data.startswith(b'\x1f\x8b'):
        raise ValueError(f'Multiple gzip members require explicit review: {source}')
    destination.write_bytes(data[:used])
    return {
        'file': str(destination),
        'source_bytes': len(data),
        'copied_bytes': used,
        'trailing_bytes_removed': len(data) - used,
        'sha256': sha256(destination),
    }


def curate_gse256200_metadata(path):
    """Resolve the two swapped tissue labels using title and sample description."""
    df = pd.read_csv(path)
    expected = {
        'GSM8087498': {'title_prefix': 'retina,', 'tissue': 'Retina', 'count_column': 'A995_RR'},
        'GSM8087499': {'title_prefix': 'choroid,', 'tissue': 'Choroid', 'count_column': 'A995LC'},
    }
    changed = []
    for gsm, rule in expected.items():
        rows = df.index[df['gsm'].eq(gsm)].tolist()
        if len(rows) != 1:
            raise ValueError(f'Expected one GSE256200 metadata row for {gsm}')
        i = rows[0]
        title = str(df.at[i, 'title']).lower()
        if not title.startswith(rule['title_prefix']) or str(df.at[i, 'count_column']) != rule['count_column']:
            raise ValueError(f'Unexpected title/count mapping for {gsm}: {df.at[i, "title"]}, {df.at[i, "count_column"]}')
        before = {'source': df.at[i, 'source'], 'tissue': df.at[i, 'tissue']}
        after = {'source': rule['tissue'], 'tissue': rule['tissue']}
        if before != after:
            df.at[i, 'source'] = after['source']
            df.at[i, 'tissue'] = after['tissue']
            changed.append({'gsm': gsm, 'before': before, 'after': after})
    df.to_csv(path, index=False)
    return {
        'file': str(path),
        'basis': 'GEO sample title plus count-column description; source/tissue fields were internally inconsistent for these two samples',
        'changed': changed,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('config.json'))
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text(encoding='utf-8'))
    source = Path(cfg['source_root'])
    root = Path(cfg['output_root'])
    dest = root / 'analysis' / 'inputs'
    cache = source / 'revision_20261014_reference_cache'
    dest.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    url = 'https://ftp.ensembl.org/pub/release-105/tsv/ensembl-compara/homologies/mus_musculus/Compara.105.protein_default.homologies.tsv.gz'
    mouse_path = cache / 'mouse_Compara105_protein.tsv.gz'
    provenance = [fetch(url, mouse_path)]
    columns = ['gene_stable_id', 'species', 'homology_type', 'homology_gene_stable_id',
               'homology_species', 'identity', 'homology_identity', 'is_high_confidence']
    mouse_parts = []
    for chunk in pd.read_csv(mouse_path, sep='\t', usecols=columns, chunksize=300000):
        sub = chunk[(chunk.species == 'mus_musculus') & (chunk.homology_species == 'gallus_gallus')].copy()
        if len(sub):
            mouse_parts.append(sub)
    if not mouse_parts:
        raise ValueError('No direct mouse-chicken homology records in the source')
    mouse_map = pd.concat(mouse_parts).rename(columns={
        'gene_stable_id': 'mouse_gene_id', 'homology_gene_stable_id': 'chick_gene_id'})
    mouse_map['ensembl_release'] = 105
    mouse_map.to_csv(dest / 'chick_mouse_all_homologies.csv', index=False)
    strict = mouse_map[mouse_map.homology_type.eq('ortholog_one2one')].copy()
    strict = strict.drop_duplicates(['chick_gene_id', 'mouse_gene_id'])
    assert not strict.chick_gene_id.duplicated().any()
    assert not strict.mouse_gene_id.duplicated().any()
    strict.to_csv(dest / 'chick_mouse_one2one.csv', index=False)
    chick_path = source / 'external_validation/GSE227724/Compara.105.protein_default.homologies.tsv.gz'
    gp = pd.read_csv(chick_path, sep='\t', usecols=columns)
    gp = gp[(gp.species == 'gallus_gallus') & (gp.homology_species == 'cavia_porcellus')].copy()
    gp = gp.rename(columns={'gene_stable_id': 'chick_gene_id', 'homology_gene_stable_id': 'gp_gene_id'})
    gp['ensembl_release'] = 105
    gp.to_csv(dest / 'chick_gp_all_homologies.csv', index=False)
    gp = gp[gp.homology_type.eq('ortholog_one2one')].drop_duplicates(['chick_gene_id', 'gp_gene_id'])
    assert not gp.chick_gene_id.duplicated().any()
    assert not gp.gp_gene_id.duplicated().any()
    gp.to_csv(dest / 'chick_gp_one2one.csv', index=False)
    features(cfg['chick_features'], 'gallus_gallus').to_csv(dest / 'chick_features.csv', index=False)
    features(cfg['mouse_features'], 'mus_musculus').to_csv(dest / 'mouse_features.csv', index=False)
    for gse, counts in [('GSE227724','GSE227724_rawCounts_dnrs_142samp.txt.gz'),
                        ('GSE261232','GSE261232_raw_counts_DNRS_144samp.txt.gz'),
                        ('GSE256200','GSE256200_count_table.txt.gz')]:
        folder = source / 'external_validation' / gse
        if not (folder / counts).exists():
            print('Counts filename requires inspection', gse, flush=True)
            continue
        if gse == 'GSE256200':
            provenance.append({
                'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256200/suppl/GSE256200_count_table.txt.gz',
                'source': str(folder / counts),
                'normalization': copy_valid_gzip(folder / counts, dest / counts),
            })
        else:
            shutil.copy2(folder / counts, dest / counts)
        shutil.copy2(folder / 'sample_metadata.csv', dest / f'{gse}_metadata.csv')
        if gse == 'GSE256200':
            provenance.append({'metadata_curation': curate_gse256200_metadata(dest / f'{gse}_metadata.csv')})
    geo_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235684/soft/GSE235684_family.soft.gz'
    provenance.append(fetch(geo_url, dest / 'GSE235684_family.soft.gz'))
    provenance.append({'url':'https://ftp.ensembl.org/pub/release-105/tsv/ensembl-compara/homologies/gallus_gallus/Compara.105.protein_default.homologies.tsv.gz',
                       'file':str(chick_path), 'sha256':sha256(chick_path)})
    (dest / 'mapping_provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print('One-to-one records:', len(strict), 'chick-mouse;', len(gp), 'chick-guinea-pig', flush=True)


if __name__ == '__main__':
    main()
