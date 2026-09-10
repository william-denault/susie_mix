"""Reading PDF of the new sections, with original vector figure pages and captions."""
from pathlib import Path
from io import BytesIO
import html
import json
import re
import hashlib
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, KeepTogether, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.colors import HexColor
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / 'tmp/manuscript_simulation_update'
DATA = json.loads((STAGE / 'simulation_sections.json').read_text(encoding='utf-8'))
TMP = ROOT / 'tmp/pdfs/manuscript_simulations'
TMP.mkdir(parents=True, exist_ok=True)

def plain_math(text):
    text = text.replace(r'\widehat{\mathrm{Var}}', 'sample Var')
    text = text.replace(r'\mathrm{Var}_n', 'sample Var')
    for cmd in ('mathrm', 'mathbf', 'mathcal', 'overline'):
        text = re.sub(r'\\' + cmd + r'\{([^{}]*)\}', r'\1', text)
    symbols = {r'\ldots': '...', r'\geq': '>=', r'\in': ' in ',
               r'\alpha': 'alpha', r'\ell': 'l', r'\epsilon': 'error',
               r'\times': ' x ', r'\pm': '+/-', r'\%': '%',
               r'\{': '{', r'\}': '}', r'\_': '_'}
    for a,b in symbols.items(): text = text.replace(a,b)
    text = text.replace('{,}', ',')
    text = re.sub(r'([_^])\{([^{}]*)\}', r'\1\2', text)
    return text.replace('{','').replace('}','')

def markup(tex):
    tex = re.sub(r'\\ref\{([^}]+)\}', lambda m: DATA['refs'][m[1]], tex)
    tex = re.sub(r'\$([^$]+)\$', lambda m: plain_math(m[1]), tex)
    tex = re.sub(r'\\texttt\{([^{}]*)\}', r'\1', tex)
    tex = re.sub(r'\\textbf\{([^{}]*)\}', r'\1', tex)
    tex = tex.replace(r'\%', '%').replace(r'\_', '_').replace('~',' ')
    tex = tex.replace('--','-').replace('``','"').replace("''",'"')
    tex = re.sub(r'\s+', ' ', tex).strip()
    if '\\' in tex: raise ValueError('Unrendered LaTeX: ' + tex)
    return html.escape(tex)

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='SimBody',fontName='Times-Roman',fontSize=12,leading=16.5,spaceAfter=9))
styles.add(ParagraphStyle(name='SimHeading',fontName='Times-Bold',fontSize=13,leading=17,spaceBefore=12,spaceAfter=7,keepWithNext=True))
styles.add(ParagraphStyle(name='SimTitle',fontName='Times-Bold',fontSize=21,leading=25,spaceAfter=14))
styles.add(ParagraphStyle(name='SimSection',fontName='Times-Bold',fontSize=17,leading=21,spaceBefore=10,spaceAfter=12,keepWithNext=True))
styles.add(ParagraphStyle(name='SimEquation',fontName='Courier',fontSize=10.3,leading=15,spaceBefore=4,spaceAfter=12,leftIndent=8))
styles.add(ParagraphStyle(name='SimNote',fontName='Times-Italic',fontSize=10.5,leading=14,spaceAfter=14,textColor=HexColor('#444444')))
styles.add(ParagraphStyle(name='SimCaption',fontName='Times-Roman',fontSize=12,leading=15.5))

story = [Paragraph('SuSiE-mix: simulation sections',styles['SimTitle']),
         Paragraph('Reading copy of the new manuscript text and captioned figures. '
                   '10 September 2026. Figure numbers match the updated LaTeX manuscript. '
                   'This file contains the simulation sections only; it is not a compilation '
                   'of the full manuscript.',styles['SimNote'])]
for part in ('results','methods','discussion'):
    if part != 'results': story.append(PageBreak())
    story.append(Paragraph(part.capitalize(),styles['SimSection']))
    for b in DATA['blocks'][part]:
        if b['kind']=='heading': story.append(Paragraph(markup(b['text']),styles['SimHeading']))
        elif b['kind']=='equation':
            story.append(Paragraph(html.escape(b['plain']).replace('\n','<br/>'),styles['SimEquation']))
        else: story.append(Paragraph(markup(b['text']),styles['SimBody']))

def footer(c,doc):
    c.saveState()
    c.setStrokeColor(HexColor('#bbbbbb')); c.line(54,751,558,751)
    c.setFont('Helvetica',8.5); c.setFillColor(HexColor('#555555'))
    c.drawString(54,762,'SuSiE-mix | Simulation sections')
    c.drawRightString(558,30,str(doc.page))
    c.restoreState()

text_path = TMP / 'text_pages.pdf'
SimpleDocTemplate(str(text_path),pagesize=(612,792),leftMargin=54,rightMargin=54,
                  topMargin=54,bottomMargin=50,title='SuSiE-mix simulation sections',
                  author='William R.P. Denault and coauthors').build(story,onFirstPage=footer,onLaterPages=footer)
writer = PdfWriter()
writer.append(PdfReader(text_path))
for f in DATA['figures']:
    source = STAGE / 'simulation_figures' / (f['stem'] + '.pdf')
    reader = PdfReader(source)
    assert len(reader.pages)==1
    figure_page = reader.pages[0]
    w,h = float(figure_page.mediabox.width),float(figure_page.mediabox.height)
    # Preserve the full original vector plot at its original physical size.
    cap = Paragraph('<b>Figure ' + f['number'] + '. ' + markup(f['title']) + '</b> ' + markup(f['caption']),styles['SimCaption'])
    _, cap_h = cap.wrap(w-80,1000)
    caption_region = cap_h+68
    page_h = h+caption_region+26
    layer = BytesIO()
    c = canvas.Canvas(layer,pagesize=(w,page_h))
    c.setFont('Helvetica',9); c.setFillColor(HexColor('#555555'))
    c.drawString(40,page_h-18,'SuSiE-mix | ' + ('Supplementary simulation figure' if f['supplementary'] else 'Main simulation figure'))
    cap.drawOn(c,40,caption_region-22-cap_h)
    c.setStrokeColor(HexColor('#bbbbbb')); c.line(40,caption_region-12,w-40,caption_region-12)
    c.setFont('Helvetica',8.5); c.drawRightString(w-40,22,str(len(writer.pages)+1))
    c.save(); layer.seek(0)
    page = PdfReader(layer).pages[0]
    page.merge_transformed_page(figure_page,Transformation().translate(tx=0,ty=caption_region))
    writer.add_page(page)
writer.add_metadata({'/Title':'SuSiE-mix: simulation sections and captioned figures',
                     '/Author':'William R.P. Denault and coauthors',
                     '/Subject':'Reading copy of the simulation update; not the full manuscript'})
output = STAGE / 'simulation_sections_review.pdf'
writer.write(output)

# Cross-reference, environment, and file-integrity checks on the edited source.
tex = (STAGE/'susie_mix_manuscript.tex').read_text(encoding='utf-8')
labels = re.findall(r'\\label\{([^}]+)\}',tex)
refs = re.findall(r'\\ref\{([^}]+)\}',tex)
assert len(labels)==len(set(labels)), 'Duplicate labels'
assert set(refs)<=set(labels), 'Unresolved figure references: '+str(set(refs)-set(labels))
assert len(re.findall(r'\\begin\{figure\}',tex))==19
stack=[]
for match in re.finditer(r'\\(begin|end)\{([^}]+)\}',re.sub(r'(?<!\\)%[^\n]*','',tex)):
    if match[1]=='begin': stack.append(match[2])
    else:
        assert stack and stack.pop()==match[2], 'Mismatched environment'
assert not stack
for path in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',tex):
    assert (STAGE/path).is_file(), 'Missing PDF '+path
for m in json.loads((STAGE/'simulation_tables/figure_manifest.json').read_text()):
    assert hashlib.sha256((STAGE/m['file']).read_bytes()).hexdigest()==m['sha256']
for obsolete in ('n=200','90 scenarios','20260908+1000003','standardization, EM estimation','Example supplement heading'):
    assert obsolete not in tex, 'Obsolete benchmark text: '+obsolete
pdf = PdfReader(output)
assert len(pdf.pages)==len(PdfReader(text_path).pages)+19
for f,page in zip(DATA['figures'],pdf.pages[-19:]):
    extracted=page.extract_text()
    assert 'Figure '+f['number']+'.' in extracted
    assert f['title'] in extracted.replace('\n',' '), f['title']
print(f'Created {len(pdf.pages)}-page reading copy; checked 19 figure paths, labels, captions and PDF hashes.')
