param([switch]$Apply)
$ErrorActionPreference = 'Stop'
$taskStage = (Resolve-Path -LiteralPath 'C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/manuscript_lead_pip_update').Path
$taskTarget = (Resolve-Path -LiteralPath 'C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix').Path
$taskRootPrefix = $taskTarget.TrimEnd('\') + '\'
$taskOriginal = Get-Content -LiteralPath (Join-Path $taskStage 'original_hashes.json') -Raw | ConvertFrom-Json -AsHashtable
$taskFiles = @('susie_mix_manuscript.tex', 'simulation_sections.json', 'simulation_sections_review.pdf', 'SIMULATION_UPDATE.md',
  'simulation_figures/lead_pip_accuracy_pooled_K.pdf', 'simulation_figures/additive_only_lead_coding_compact.pdf',
  'simulation_figures/additive_only_pip_allocation_compact.pdf')
$taskFiles += Get-ChildItem -LiteralPath (Join-Path $taskStage 'simulation_tables') -File | ForEach-Object { 'simulation_tables/' + $_.Name }
$taskChanges = @()
foreach ($taskRelative in $taskFiles) {
  $taskSource = Join-Path $taskStage $taskRelative
  $taskDestination = [IO.Path]::GetFullPath((Join-Path $taskTarget $taskRelative))
  if (-not $taskDestination.StartsWith($taskRootPrefix, [StringComparison]::OrdinalIgnoreCase)) { throw "Path outside manuscript: $taskDestination" }
  $taskNewHash = (Get-FileHash -LiteralPath $taskSource -Algorithm SHA256).Hash.ToLowerInvariant()
  if (Test-Path -LiteralPath $taskDestination) {
    $taskOldHash = (Get-FileHash -LiteralPath $taskDestination -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($taskOldHash -eq $taskNewHash) { continue }
    if (-not $taskOriginal.ContainsKey($taskRelative) -or $taskOriginal[$taskRelative] -ne $taskOldHash) {
      throw "Manuscript changed since review, or unexpected existing file: $taskRelative"
    }
  } elseif ($taskOriginal.ContainsKey($taskRelative)) { throw "Expected manuscript file disappeared: $taskRelative" }
  $taskChanges += [pscustomobject]@{ File = $taskRelative; Source = $taskSource; Destination = $taskDestination; SHA256 = $taskNewHash }
}
if (-not $Apply) {
  $taskChanges | Select-Object File
  Write-Output "Ready to update $($taskChanges.Count) files; original contents are unchanged."
  return
}
foreach ($taskChange in $taskChanges) {
  if (Test-Path -LiteralPath $taskChange.Destination) {
    $taskBackup = Join-Path (Join-Path $taskStage 'before_apply') $taskChange.File
    New-Item -ItemType Directory -Path (Split-Path -Parent $taskBackup) -Force | Out-Null
    Copy-Item -LiteralPath $taskChange.Destination -Destination $taskBackup -Force
  }
  New-Item -ItemType Directory -Path (Split-Path -Parent $taskChange.Destination) -Force | Out-Null
  Copy-Item -LiteralPath $taskChange.Source -Destination $taskChange.Destination -Force
  if ((Get-FileHash -LiteralPath $taskChange.Destination -Algorithm SHA256).Hash.ToLowerInvariant() -ne $taskChange.SHA256) {
    throw "Copied file failed verification: $($taskChange.File)"
  }
}
$taskChanges | ConvertTo-Json -Depth 3 | Set-Content -LiteralPath (Join-Path $taskStage 'applied_manifest.json')
Write-Output "Updated and verified $($taskChanges.Count) files in $taskTarget."
