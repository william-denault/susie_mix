"""Reuse the established vector review renderer for the expanded sections."""
from pathlib import Path
source = Path(__file__).with_name('render_manuscript_simulations.py').read_text(encoding='utf-8')
changes = {
    "'tmp/manuscript_simulation_update'": "'tmp/manuscript_lead_pip_update'",
    "'tmp/pdfs/manuscript_simulations'": "'tmp/pdfs/manuscript_lead_pip'",
    'KeepTogether, PageBreak': 'KeepTogether, PageBreak, Table, TableStyle',
    '10 September 2026.': '11 September 2026.',
    "fontSize=12,leading=16.5,spaceAfter=9": "fontSize=11.5,leading=15.5,spaceAfter=9",
    "r'\\times': ' x ', r'\\pm': '+/-', r'\\%': '%'": "r'\\times': ' x ', r'\\pm': '+/-', r'\\%': '%'",
    "==19": "==len(DATA['figures'])",
    "+19": "+len(DATA['figures'])",
    "pdf.pages[-19:]": "pdf.pages[-len(DATA['figures']):]",
    "assert f['title'] in extracted.replace('\\n',' '), f['title']": "assert html.unescape(markup(f['title'])) in extracted.replace('\\n',' '), f['title']",
    "checked 19 figure paths, labels, captions and PDF hashes.": "checked all figure paths, labels, captions and PDF hashes.",
}
for old,new in changes.items():
    assert old in source, old
    source = source.replace(old,new)
needle = "        elif b['kind']=='equation':"
insert = '''        elif b['kind']=='table':
            story.append(Paragraph('<b>Table '+b['number']+'. '+markup(b['title'])+'</b>',styles['SimHeading']))
            grid = Table([b['columns']]+b['rows'],colWidths=[105,75,95,110,95],hAlign='LEFT',repeatRows=1)
            grid.setStyle(TableStyle([
                ('FONTNAME',(0,0),(-1,0),'Times-Bold'),('FONTNAME',(0,1),(-1,-1),'Times-Roman'),
                ('FONTSIZE',(0,0),(-1,-1),10),('LEADING',(0,0),(-1,-1),13),
                ('ALIGN',(1,0),(-1,-1),'RIGHT'),('BOTTOMPADDING',(0,0),(-1,-1),6),
                ('TOPPADDING',(0,0),(-1,-1),6),('LINEABOVE',(0,0),(-1,0),.6,HexColor('#555555')),
                ('LINEBELOW',(0,0),(-1,0),.4,HexColor('#999999')),
                ('LINEBELOW',(0,-1),(-1,-1),.6,HexColor('#555555'))]))
            story.append(grid)
            story.append(Spacer(1,8))
            story.append(Paragraph(markup(b['text']),styles['SimNote']))
        elif b['kind']=='equation':'''
assert source.count(needle)==1
source = source.replace(needle,insert)
exec(compile(source,str(Path(__file__).resolve()),'exec'))
