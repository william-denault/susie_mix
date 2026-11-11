# Validation — 2026-11-10

The user reported iteration 024, batch 25, 117/708 committed stages, 257/299
validated genes in the current batch, and no workflow error. The installer
checks the complete installed code manifest against the released 2026-11-09
extension before modifying any configuration.

Executed checks:

- 29 Python upgrade tests using the full 18,263-gene, 708-stage baseline fixture:
  exact 215-stage prefix preservation, unchanged batch membership and ordering,
  unchanged runtime estimates, all-gene coverage, the requested rho boundaries,
  new full 029/030 stages, total 400 stages, and controller termination at 030.
- Current or completed 025 is accepted because its definition is unchanged.
  Any prepared/submitted 026 work is rejected, including an empty owned folder
  or a task that has not yet produced results.
- Live-controller hold/drain/release, idle continuation, no cancellation,
  confirmation of the recovery job, waiting for the previous upgrade watchdog,
  interrupted code/configuration/plan installation, unchanged scientific state,
  duplicate launch protection, and resume using the installed frozen schedule.
- A normal in-flight submission can finish while the gate drains workers. An
  unresolved submission or unrelated error after draining is not cleared.
- 64 actual susieRSlidePrior fits on four synthetic genes through iterations
  015–030 with reduced batch counts, including the 0.8/0.5/0.3 sequence and both
  final full EM steps. Exact sufficient-statistic updates, previous-fit/prior
  warm-start validation, T/F gene names, history fields, and actual current/previous
  result cleanup were checked.
- Summary metadata checks for phase, rho, full-EM batch zero, input/update roles,
  historical versus completed/planned/running labels and submission ordering.
- All executable runtime files, scheduler code, R engine and scientific helpers
  are byte-identical to the 2026-11-09 runtime. Only its README is updated.
- ZIP integrity, nested SHA256 manifests, LF line endings, Python 3.6 parsing,
  and syntax checks for all shell launchers and job scripts.

Slurm interactions in tests are mocked. The model tests use local synthetic
data and do not measure GTEx convergence or RCC performance. No live RCC jobs
were submitted or cancelled during local validation. Completing the requested
schedule does not itself certify convergence.
