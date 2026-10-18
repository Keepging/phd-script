"""Fetch canonical sequences for all accessions in site_meta.tsv from UniProt REST.
Writes uniprot_sequences.tsv (accession, entry_version, sequence_version, seq_modified, sequence_length, sequence)
and uniprot_release.json (UniProt release + date from response headers)."""
import csv, json, sys, time, requests, pathlib

meta_path, out_dir = sys.argv[1], pathlib.Path(sys.argv[2])
accs = sorted({r['protein'] for r in csv.DictReader(open(meta_path), delimiter='\t')})
print(len(accs), 'accessions')

S = requests.Session()
S.headers['User-Agent'] = 'phd-script/kinase-scoring (python-requests)'
BASE = 'https://rest.uniprot.org/uniprotkb/search'
fields = 'accession,sequence,sequence_version,version,date_sequence_modified,date_modified,reviewed,id,gene_primary,organism_id'
release = None
rows = {}
BATCH = 100
for i in range(0, len(accs), BATCH):
    chunk = accs[i:i+BATCH]
    q = ' OR '.join(f'accession:{a}' for a in chunk)
    for attempt in range(5):
        try:
            r = S.get(BASE, params={'query': q, 'fields': fields, 'format': 'tsv', 'size': 500}, timeout=120)
            r.raise_for_status(); break
        except Exception as e:
            print('retry', i, attempt, e, file=sys.stderr); time.sleep(2**attempt)
    else:
        sys.exit('failed batch %d' % i)
    if release is None:
        release = {'x-uniprot-release': r.headers.get('x-uniprot-release'),
                   'x-uniprot-release-date': r.headers.get('x-uniprot-release-date'),
                   'x-api-deployment-date': r.headers.get('x-api-deployment-date'),
                   'fetched_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
        print(release)
    rd = csv.DictReader(r.text.splitlines(), delimiter='\t')
    for row in rd:
        # entries can be returned under a secondary accession; keep mapping to requested acc
        rows[row['Entry']] = row
    print(i, len(rows), flush=True)

# handle any requested accession not returned as primary (merged/secondary accessions)
missing = [a for a in accs if a not in rows]
print('missing as primary:', len(missing), missing[:20])
sec_map = {}
for a in missing:
    r = S.get(BASE, params={'query': f'accession:{a}', 'fields': fields, 'format': 'tsv'}, timeout=60)
    rd = list(csv.DictReader(r.text.splitlines(), delimiter='\t'))
    if rd:
        sec_map[a] = rd[0]['Entry']; rows.setdefault(rd[0]['Entry'], rd[0])
    else:
        print('NOT FOUND', a)

with open(out_dir/'uniprot_sequences.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['requested_accession','uniprot_entry','entry_name','gene_primary','reviewed','organism_id','entry_version','sequence_version','date_sequence_modified','date_modified','length','sequence'])
    for a in accs:
        e = sec_map.get(a, a)
        row = rows.get(e)
        if row is None:
            w.writerow([a, '', '', '', '', '', '', '', '', '', 0, '']); continue
        w.writerow([a, row['Entry'], row['Entry Name'], row.get('Gene Names (primary)',''), row['Reviewed'], row['Organism (ID)'],
                    row['Entry version'], row['Sequence version'], row['Date of last sequence modification'], row['Date of last modification'],
                    len(row['Sequence']), row['Sequence']])
release['secondary_accession_map'] = sec_map
json.dump(release, open(out_dir/'uniprot_release.json','w'), indent=1)
print('done')
