from pathlib import Path
root = Path(r"C:\Document\Serieux\Travail\Package\git\susieR")
source = (root / "inst/examples/benchmark_slider_prior.R").read_text()
old_cases = '''cases <- data.frame(n=c(500L,500L,1000L),p=c(1000L,1000L,3000L),
                    variance=c("fixed","estimated","fixed"))'''
new_cases = '''cases <- data.frame(n=c(500L,500L),p=c(1000L,1000L),
                    variance=c("fixed","residual_only"))'''
assert old_cases in source
source = source.replace(old_cases,new_cases)
source = source.replace('estimated <- cfg$variance=="estimated"',
                        'estimated <- cfg$variance=="residual_only"')
source = source.replace('estimate_prior_variance=estimated,estimate_residual_variance=estimated',
                        'estimate_prior_variance=FALSE,estimate_residual_variance=estimated')
source = source.replace('"Fixed: neither Gaussian variance estimated. Estimated: both estimated, prior method optim.",',
    '"Fixed: neither variance estimated. Residual_only: sigma2 updated after each sweep; all effect-prior V values fixed.",')
source = source.replace('"The new slider probabilities remain fixed at 1/17 in both variance settings.",',
    '"The new slider probabilities remain fixed at 1/17 in both variance settings. No prior-variance optimization is performed.",')
target = Path(r"C:\Document\Serieux\Travail\Data_analysis_and_papers\susie_mix\tmp\slider_prior_impl\benchmark_residual_variance.R")
target.write_text(source)
print(target)
