# Updated slider manuscript

`slider_model.tex` is the complete updated document. Section 5, "From slider SER to variational empirical Bayes", is embedded directly in it; no `\input` of another source is required.

The new section is also provided separately as `slider_veb_section.tex` for convenient copying into another version of the manuscript. It uses the existing preamble macros and the bibliography key `wang2018susie`, included in the full document.

`additive_delta_histograms.png` is the original embedded histogram image extracted from page 10 of the supplied `slider_model-4.pdf`. Keep it beside `slider_model.tex`. The optional `output/slider_genotype_boxplots.png` retains the original conditional inclusion and is not required.

Compile from this directory in your existing LaTeX installation, or upload the ZIP contents to Overleaf and select `slider_model.tex` as the main document:

```text
pdflatex slider_model.tex
pdflatex slider_model.tex
```

The bibliography is included in the document, so BibTeX is not required. The original package requirements are retained, including `bbold`.

The update adds the model definition, two propositions and proofs, a computable conditional ELBO, the residual SER updates, variance updates, and the distinction between component-specific and shared sliders. The introduction and one-SNP scope paragraph have been reconciled with the extension. Minor wording errors and the endpoints typo have also been corrected. Existing simulation numbers are retained; those simulations were not rerun for this writing task.

Validation: all 63 labels are unique; references, citations, braces, and environments pass source checks. Independent numerical checks on 30 designs (360 component updates) agree with dense Gaussian integration and the block-ELBO identity within 1.5e-14, and the ELBO stays below the exactly enumerated full-model log evidence. These checks support the formulas; they do not establish calibration or a global optimum. A local LaTeX compiler was unavailable, so the document has not been compiled in this session.
