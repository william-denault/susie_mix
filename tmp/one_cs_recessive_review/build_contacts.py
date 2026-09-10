from pathlib import Path
import csv
from PIL import Image, ImageDraw, ImageFont
root=Path(r"C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix")
src=root/"plot/one_cs_recessive"
out=root/"tmp/one_cs_recessive_review"
out.mkdir(parents=True,exist_ok=True)
rows=sorted(csv.DictReader((src/"fit_mix_one_cs_recessive_plot_summary.csv").open(encoding="utf-8-sig")),key=lambda r:(r["gene"],r["tissue"]))
font=ImageFont.truetype("C:/Windows/Fonts/arial.ttf",22)
for start in range(0,len(rows),6):
    canvas=Image.new("RGB",(2100,1300),"white")
    draw=ImageDraw.Draw(canvas)
    for k,r in enumerate(rows[start:start+6]):
        x,y=(k%3)*700,(k//3)*650
        img=Image.open(src/Path(r["output_file"]).name).convert("RGB")
        img.thumbnail((700,615))
        canvas.paste(img,(x+(700-img.width)//2,y+33))
        draw.text((x+5,y+3),f'{start+k+1} {r["gene"]} / {r["tissue"]}',fill="black",font=font)
    canvas.save(out/f"contact_{start//6+1:02d}.jpg",quality=92)
print(f"Created {(len(rows)+5)//6} contact sheets for {len(rows)} plots")

