# Slider 015–025 transition audit — 2026-11-06

Audited runtime: output/slider_online_queue_patch_2026-11-06/runtime.
RCC confirmed identical minibatch.py and engine.R SHA256 digests, with seven
committed stages and batch 8 still fitting. The scientific engine was unchanged.

## Batch sizes

Every online pass partitions all 18,263 genes exactly once: batches 1–2 contain
2,030 genes each; batches 3–9 contain 2,029 each. There are 54 online updates
across 015–020, then five full-cohort updates in 021–025.

## Verified warm-start and cleanup behavior

| New iteration | Warm starts | Results retired |
|---|---|---|
| 015 | Certified/pinned 014 | 013 and earlier |
| 016 | Each gene's 015 fit and its own batch input prior | Source 014 links and 014 results, after all nine 015 batches commit |
| 017–020 | Each gene's immediately previous iteration | Results older than the previous iteration |
| 021 | All 020 fits, with their individual batch priors | 019 and earlier |
| 022–025 | Immediately previous full-cohort step | Results older than the previous iteration |

The new fit uses the current stage's frozen prior. Validation of the warm start
uses the prior under which that particular gene was fitted, not the last prior
of its previous scan. State advances only through a durable commit. Recovered
commits are idempotent. Completed genes are adopted from saved outputs/receipts.

Deletion leaves all history, count summaries and state metadata intact. At 025,
024 and 025 result files remain, with no 026 stage. The final queued controller
exits on DONE.json without submitting further work. Result files in 025 were
fitted under that step's input prior; its after-update prior is also recorded.

## Tests performed

- Full 18,263-gene plan: membership, coverage, exact sizes, numbering through 025.
- Combined R engine and actual Python cleanup across all 59 stages, with 18
  fixture genes, two tissues, renamed IDs and literal T/F identifiers.
- 99 actual susieRSlidePrior fits across all 59 stages on nine tiny synthetic
  genes. Batch priors differ, and per-gene prior provenance is validated.
- Source links kept through the ninth update; old links deleted before 016
  workers; latest warm-start files and priors protected at every stage.
- Interrupted cleanup after partial deletion at 016 and 017; repeated cleanup
  resumes safely and changes no prior/history/count bytes.
- Interrupted durable commits at all online pass boundaries and at full steps.
- Online/full history columns, final state 59, and idempotent publication.
- Guard checks for active legacy jobs, unknown files, directories, and final
  controller completion without extra submissions.

Logs are tmp/transition_integrated_tests.log, tmp/transition_real_tests.log,
tmp/transition_orphan_integration.log, tmp/transition_cleanup_tests.log and
tmp/cleanup_upgrade_tests.log. The expected missing-receipt warning in the R
fixture comes from a negative test, not a failed audit.

## Concrete issue found and fixed in the candidate

Old cleanup rejected every non-.rds file, including .em_write_<alphanumeric>
files left by a worker killed during atomic saving. This can stop a later scan
even when every final result validates. The test reproduced that rejection.

RCC listed three such files in active iteration 015. Their presence during
fitting does not establish that they are abandoned; none should be removed from
the current/previous iteration. No such files were listed in 014, so this did
not reveal an immediate 015→016 cleanup obstruction.

The hardening recognizes the writer's reserved temporary format only while
retiring obsolete results. Other unexpected files/directories/symlinks remain
blocked. Tests inject real R tempfile names at every transition through 025;
cleanup succeeds while retaining current/previous outputs and all metadata.

The hardened candidate is packaged separately as slider_online_cleanup_patch;
no running RCC files were changed by this audit. Its installation waits for a
safe pause and includes a recovery job. All statistical R code is byte-identical
to the confirmed RCC runtime. Local tests cannot guarantee RCC uptime, disk
capacity, or convergence on all real data.
