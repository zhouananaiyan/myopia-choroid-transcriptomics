"""Run the R analysis and save stdout and stderr to analysis/run_log.txt."""
import argparse
import subprocess
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('root',type=Path)
p.add_argument('--rscript',required=True)
p.add_argument('--r-library',default='')
a=p.parse_args()
script=Path(__file__).with_name('03_reanalyse.R')
command=[a.rscript,str(script),str(a.root)]
if a.r_library:command.append(a.r_library)
log=a.root/'analysis'/'run_log.txt'
with log.open('w',encoding='utf-8') as handle:
    result=subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT)
print('R exit status:',result.returncode,'Log:',log)
raise SystemExit(result.returncode)
