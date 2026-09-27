from pathlib import Path

path = Path(r"C:\Document\Serieux\Travail\Package\git\susieR\vignettes\susie_slide_prior.tex")
text = path.read_text(encoding="utf-8")
changes = {
    """dataset across methods and settings. Slider probabilities remain
$1/17$ in every fit. The slider settings are""":
    """dataset across methods and settings. Finite-prior fits keep every slider
probability at $1/17$; the original slider uses continuous plug-in estimates.
The slider settings are""",
    r""" b_0&=\bar y-\sum_j\bar x_j b_j^x-\sum_j\bar h_j b_j^h,&
 \widehat y_i&=b_0+\sum_jx_{ij}b_j^x+\sum_j\ind\{x_{ij}=1\}b_j^h.""":
    r""" b_0&=\bar y-\sum_j\bar x_j b_j^x-\sum_j\bar h_j b_j^h,\nonumber\\
 \widehat y_i&=b_0+\sum_jx_{ij}b_j^x+\sum_j\ind\{x_{ij}=1\}b_j^h.""",
    """Reversing the counted allele maps $(\\beta,\\delta)$ to $(-\\beta,-\\delta)$,
with an intercept adjustment when needed. An asymmetric prior requires
a consistent allele convention. To enforce invariance to allele flips,""":
    """With centered coding or a fitted intercept, reversing the counted allele
maps $(\\beta,\\delta)$ to $(-\\beta,-\\delta)$, with an intercept adjustment
on the original scale. Without an intercept, raw-coding allele flips need
not preserve the model. An asymmetric prior requires a consistent allele
convention. In the intercept-adjusted setting, to enforce invariance to allele flips,""",
    """The conditional coefficient distribution given a SNP is a finite
Gaussian mixture. A single Gaussian with matched moments generally""":
    """Conditional moments for zero-prior pairs are only storage placeholders;
inference uses pairs with positive posterior mass.
The conditional coefficient distribution given a SNP is a finite
Gaussian mixture. A single Gaussian with matched moments generally""",
}
for old, new in changes.items():
    assert old in text, old[:100]
    text = text.replace(old, new, 1)
path.write_text(text, encoding="utf-8")
print("Clarified benchmark priors, conditional moments and allele orientation; split the prediction display.")
