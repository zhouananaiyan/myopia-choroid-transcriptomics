import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';

const root=process.argv[2], number=process.argv[3];
const specs=JSON.parse(await fs.readFile(path.join(root,'analysis/workbook_data',`Additional_file_${number}.json`),'utf8'));
const wb=Workbook.create();
const previews=path.join(root,'qa','sheets');
await fs.mkdir(previews,{recursive:true});
const safe=x=>typeof x==='string'&&x.startsWith('=')?"'"+x:x;
for(const spec of specs){
  const sheet=wb.worksheets.add(spec.name);
  const nr=spec.rows.length+1,nc=spec.columns.length;
  sheet.getRangeByIndexes(0,0,1,nc).values=[spec.columns];
  for(let i=0;i<spec.rows.length;i+=4000){
    const block=spec.rows.slice(i,i+4000).map(row=>row.map(safe));
    sheet.getRangeByIndexes(i+1,0,block.length,nc).values=block;
  }
  const used=sheet.getRangeByIndexes(0,0,nr,nc);
  used.format.font.name='Arial';used.format.font.size=10;
  used.format.columnWidth=20;used.format.rowHeight=17;
  const header=sheet.getRangeByIndexes(0,0,1,nc);
  header.format.fill='#E6ECEF';header.format.font.bold=true;
  header.format.rowHeight=44;header.format.wrapText=true;
  sheet.freezePanes.freezeRows(1);
  for(let j=0;j<nc;j++){
    const name=spec.columns[j];
    const column=sheet.getRangeByIndexes(1,j,Math.max(1,nr-1),1);
    if(/PValue|p_value|FDR|MWU/.test(name))column.setNumberFormat('0.000E+00');
    else if(/logFC|logCPM|effect|correlation|Correlation|jaccard|agreement|precision|recall/.test(name))column.setNumberFormat('0.0000');
    else if(/gene_id|stable_id/.test(name))sheet.getRangeByIndexes(0,j,nr,1).format.columnWidth=25;
  }
  if(spec.name==='README'){
    sheet.getRangeByIndexes(0,0,nr,1).format.columnWidth=35;
    sheet.getRangeByIndexes(0,1,nr,1).format.columnWidth=108;
    used.format.wrapText=true;used.format.rowHeight=42;
  }
  const lastCol=Math.min(nc,7);
  const range=sheet.getRangeByIndexes(0,0,Math.min(nr,7),lastCol);
  // Render the header and up to six data rows for layout inspection.
  try{
    const end=String.fromCharCode(64+lastCol)+Math.min(nr,7);
    const pic=await wb.render({sheetName:spec.name,range:`A1:${end}`,scale:1,format:'png'});
    await fs.writeFile(path.join(previews,`${number}_${spec.name}.png`),new Uint8Array(await pic.arrayBuffer()));
  }catch(e){console.log('PREVIEW_FAILED',spec.name,String(e).slice(0,160));}
  console.log('Added',spec.name,nr-1,'rows');
}
wb.recalculate();
const first=wb.worksheets.getItem(specs[0].name);
const check={sheet:specs[0].name,rows:specs[0].rows.length+1,columns:specs[0].columns.length,
  firstRow:first.getRangeByIndexes(0,0,1,specs[0].columns.length).values[0]};
await fs.writeFile(path.join(previews,`book${number}_inspect.json`),JSON.stringify(check,null,2));
const output=await SpreadsheetFile.exportXlsx(wb);
const dest=path.join(root,'submission',`Additional_file_${number}.xlsx`);
await fs.mkdir(path.dirname(dest),{recursive:true});
await output.save(dest);
console.log('SAVED',dest);
