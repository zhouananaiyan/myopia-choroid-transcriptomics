import json
import sys
from pathlib import Path
import pandas as pd
import requests

root = Path(sys.argv[1])
out = root / 'analysis' / 'results'
df = pd.read_csv(out / 'cross_species_common_genes.csv')
payload = {
    'organism':'mmusculus',
    'query': {'opposite_candidates':df.loc[df.candidate,'mouse_gene_id'].tolist(),
              'recovery_up':df.loc[df.module.eq('Recovery_up'),'mouse_gene_id'].tolist(),
              'recovery_down':df.loc[df.module.eq('Recovery_down'),'mouse_gene_id'].tolist()},
    'sources':['GO:BP','GO:MF','GO:CC','REAC','KEGG','WP'],
    'background':df.mouse_gene_id.tolist(),
    'domain_scope':'custom',
    'significance_threshold_method':'fdr',
    'user_threshold':0.05,
    'all_results':True,
    'no_evidences':False
}
(out/'enrichment_request.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
response=requests.post('https://biit.cs.ut.ee/gprofiler/api/gost/profile/',json=payload,timeout=180)
response.raise_for_status()
result=response.json()
(out/'enrichment_response.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
data=pd.DataFrame(result.get('result',[]))
if data.empty:
    data=pd.DataFrame(columns=['query','source','native','name','p_value','significant'])
for col in data.columns:
    if data[col].map(lambda x:isinstance(x,(list,dict))).any():
        data[col]=data[col].map(lambda x:json.dumps(x) if isinstance(x,(list,dict)) else x)
data.to_csv(out/'enrichment_complete.csv',index=False)
print(data.groupby('query').size().to_string() if len(data) else 'No terms returned')
if len(data):
    print(data.loc[data.p_value<.05,['query','source','native','name','p_value']].head(12).to_string(index=False))
