import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root=process.cwd();
const out=path.join(root,'outputs/gtex_gene_cohort_audit_2026-10-26');
const d=JSON.parse(await fs.readFile(path.join(root,'tmp/gtex_cohort_workbook_data.json'),'utf8'));
const wb=Workbook.create();
const summary=wb.worksheets.add('Summary');
const missing=wb.worksheets.add('Missing symbols');
const extras=wb.worksheets.add('Extra symbols');
const dup=wb.worksheets.add('Duplicate symbols');
const ink='#213547', light='#E9EFF5';
function base(s,range) {
  s.showGridLines=false;
  s.getRange(range).format.font={name:'Arial',size:11,color:ink};
  s.getRange(range).format.rowHeight=23;
  s.getRange(range).format.verticalAlignment='center';
}
function title(s,text,end) {
  s.getRange('A2').values=[[text]];
  s.getRange('A2').format.font={name:'Arial',size:16,bold:true,color:ink};
  s.getRange(`A2:${end}2`).format.borders={bottom:{style:'thin',color:'#A8BAC9'}};
  s.getRange('A2').format.rowHeight=31;
}
function header(s,range) {
  s.getRange(range).format={fill:ink,font:{name:'Arial',size:11,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,rowHeight:34};
}
function dataSheet(s,name,titleText,headers,rows,widths) {
  const end=String.fromCharCode(64+headers.length), bottom=5+rows.length;
  base(s,`A1:${end}${bottom}`); title(s,titleText,end);
  s.getRange('A3').values=[['Reference: GENCODE v26. Comparison: Ensembl 103 and the frozen EM gene list.']];
  s.getRange('A3').format.font={italic:true,color:'#5F6F7B'};
  s.getRange(`A5:${end}5`).values=[headers];
  s.getRange(`A6:${end}${bottom}`).values=rows;
  const table=s.tables.add(`A5:${end}${bottom}`,true,name);
  table.showFilterButton=true;
  header(s,`A5:${end}5`);
  widths.forEach((width,i)=>s.getRange(`${String.fromCharCode(65+i)}1:${String.fromCharCode(65+i)}${bottom}`).format.columnWidth=width);
  s.freezePanes.freezeRows(5);
  s.freezePanes.freezeColumns(1);
  return bottom;
}
const mr=d.missing.map(r=>[r.gtex_v26_symbol,r.reason,r.in_local_expression_symbols,r.ensembl103_symbols,r.ensembl103_biotypes,r.chromosomes,r.gencode_v26_ids,r.reference_gene_ids,r.note]);
const mend=dataSheet(missing,'MissingSymbols','Genes absent from the current scan',
  ['GTEx v26 symbol','Reason','In expression list','Ensembl 103 symbol','Ensembl 103 biotype','Chr','GENCODE v26 IDs','ID count','Notes'],mr,[23,44,19,25,40,10,51,11,100]);
missing.getRange(`H6:H${mend}`).setNumberFormat('0');
for(let i=0;i<mr.length;i++) if(mr[i][8]) {
  missing.getRange(`I${i+6}`).format.wrapText=true;
  missing.getRange(`A${i+6}:I${i+6}`).format.rowHeight=54;
}
const er=d.extras.map(r=>[r.current_symbol,r.gencode_v26_biotypes,r.ensembl103_biotypes,r.gencode_v26_ids,r.ensembl103_ids]);
const eend=dataSheet(extras,'ExtraSymbols','Selected genes outside the v26 protein-coding set',
  ['Selected symbol','GENCODE v26 biotype','Ensembl 103 biotype','GENCODE v26 IDs','Ensembl 103 IDs'],er,[24,40,27,49,49]);
const dr=d.duplicates.map(r=>[r.symbol,r.reference_gene_ids,r.ids_minus_symbols,r.in_current_scan,r.gencode_v26_ids]);
const dend=dataSheet(dup,'DuplicateSymbols','Reference symbols shared by multiple gene IDs',
  ['Symbol','Gene IDs','IDs minus symbols','In current scan','GENCODE v26 IDs'],dr,[25,16,22,22,80]);
dup.getRange(`B6:C${dend}`).setNumberFormat('0');

base(summary,'A1:C38'); title(summary,'GTEx EM gene cohort audit','C'); summary.tabColor=ink;
summary.getRange('A3').values=[['26 October 2026. Analysis of synced metadata; no changes to the EM run.']];
summary.getRange('A3').format.font={italic:true,color:'#5F6F7B'};
summary.getRange('A1:A38').format.columnWidth=56;
summary.getRange('B1:B38').format.columnWidth=18;
summary.getRange('C1:C38').format.columnWidth=93;
summary.getRange('A5:C5').values=[['Count reconciliation','Count','Definition']]; header(summary,'A5:C5');
summary.getRange('A6:C12').values=[
  ['Autosomal protein-coding reference',18901,'GENCODE v26 gene IDs on chromosomes 1–22.'],
  ['Repeated symbols',null,'Several reference IDs share a gene symbol.'],
  ['Distinct reference symbols',null,'The pipeline uses gene symbols as its keys.'],
  ['Reference symbols absent from scan',null,'The complete list is on Missing symbols.'],
  ['Selected symbols with another v26 biotype',null,'Ensembl 103 classifies these as protein-coding.'],
  ['Reconciled current scan',null,'18,901 − 17 − 1,334 + 67 = 17,617.'],
  ['Frozen scan after excluding X/Y/MT',17617,'18,468 original symbols − 851 chromosome exclusions.'],
];
summary.getRange('B7').formulas=[[`=-SUM('Duplicate symbols'!C6:C${dend})`]];
summary.getRange('B8').formulas=[['=SUM(B6:B7)']];
summary.getRange('B9').formulas=[[`=-COUNTA('Missing symbols'!A6:A${mend})`]];
summary.getRange('B10').formulas=[[`=COUNTA('Extra symbols'!A6:A${eend})`]];
summary.getRange('B11').formulas=[['=SUM(B8:B10)']];
summary.getRange('A11:C11').format.fill=light;
summary.getRange('A11:C11').format.font.bold=true;
summary.getRange('A14:C14').values=[['Why reference symbols are absent','Symbols','Interpretation']]; header(summary,'A14:C14');
const categories=[
  ['Renamed protein-coding gene','GTEx old names do not match Ensembl 103 names. CTAGE5/MIA2 needs ID review.'],
  ['Absent from GTEx expression list','Not a low-expression test. These names are absent from the expression rows.'],
  ['Reclassified as non-protein-coding','The same stable ID has a different biotype in Ensembl 103.'],
  ['Stable ID absent from Ensembl 103','Old IDs are not present in release 103. Merge/retirement histories unresolved.'],
];
summary.getRange('A15:C18').values=categories.map(([label,note])=>[label,null,note]);
summary.getRange('B15:B18').formulas=categories.map((_,i)=>[`=COUNTIF('Missing symbols'!B6:B${mend},A${15+i})`]);
summary.getRange('A20:C23').values=[
  ['Validation','Result','Evidence'],
  ['Ensembl 103 symbol filter','Exact match','Reproduces all 18,468 original manifest symbols with zero differences.'],
  ['Public GTEx expression symbol set','Exact match','All 54,592 symbols match the local expression list.'],
  ['Scope of this audit','Gene selection','No posterior fits were rerun. No scheduler or EM settings were changed.'],
]; header(summary,'A20:C20');
summary.getRange('A25:C25').values=[['Sources and limits','','']];
summary.getRange('A25:C25').format.fill=light;
summary.getRange('A25').format.font.bold=true;
summary.getRange('A26:C33').values=[
  ['GTEx reference API','v26',d.summary.public_sources['GTEx reference gene API']],
  ['Ensembl gene annotation','103',d.summary.public_sources['Ensembl 103 reference']],
  ['GTEx expression row metadata','v8',d.summary.public_sources['GTEx v8 expression metadata']],
  ['GTEx collapsed model rules','',d.summary.public_sources['GTEx collapsed gene model rules']],
  ['Local frozen gene lists','','genes_all.txt, genes_protein_coding.txt, iteration_000/manifest.csv, chromosome_exclusions.csv'],
  ['RCC custom annotation','','Not synced locally. The public reference reproduces its selected gene set exactly.'],
  ['Absent expression rows','','The precise model-construction reason was not reconstructed for each excluded ID.'],
  ['Counts describe selected entries','','A result file can hold a QC exclusion. Synced metadata flags ZBTB9 with no eligible tissues.'],
];
summary.getRange('C26:C30').format.wrapText=true;
summary.getRange('A26:C30').format.rowHeight=46;
summary.getRange('B6:B18').setNumberFormat('#,##0;[Red](#,##0);0');
summary.getRange('B6:B18').format.horizontalAlignment='right';
summary.getRange('B21:B23').format.horizontalAlignment='left';
wb.recalculate();
const inspect=await wb.inspect({kind:'table',range:'Summary!A6:C18',include:'values,formulas',tableMaxRows:15,tableMaxCols:3,maxChars:6500});
console.log(inspect.ndjson);
const values=summary.getRange('B6:B18').values;
if(values[5][0]!==17617 || values[9][0]!==679 || values[10][0]!==505 || values[11][0]!==113 || values[12][0]!==37) throw new Error(`Reconciliation failed: ${JSON.stringify(values)}`);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'Formula errors',maxChars:2000});
console.log(errors.ndjson);
for (const [sheetName,range,name] of [['Summary','A1:C23','summary'],['Missing symbols','A1:H14','missing'],['Extra symbols','A1:E13','extras'],['Duplicate symbols','A1:E14','duplicates']]) {
  const preview=await wb.render({sheetName,range,scale:1,format:'png'});
  await fs.writeFile(path.join(out,`${name}_preview.png`),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(out,'gtex_missing_genes_audit.xlsx'));
console.log('Exported gtex_missing_genes_audit.xlsx');
