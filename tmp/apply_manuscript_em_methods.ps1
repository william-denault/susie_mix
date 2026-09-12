param([switch]$Apply)
$ErrorActionPreference = 'Stop'
$taskStage = (Resolve-Path -LiteralPath 'C:/Document/Serieux/Travail/Data_analysis_and_papers/susie_mix/tmp/manuscript_em_methods_update').Path
$taskTarget = (Resolve-Path -LiteralPath 'C:/Document/Serieux/Travail/Package/git/Science-SuSiE-mix').Path
$taskRootPrefix = $taskTarget.TrimEnd('\') + '\'
$taskOriginal = Get-Content -LiteralPath (Join-Path $taskStage 'original_hashes.json') -Raw | ConvertFrom-Json -AsHashtable
$taskFiles = @('susie_mix_manuscript.tex', 'library.bib')
$taskChanges = @()
foreach ($taskRelative in $taskFiles) {
  $taskSource = Join-Path $taskStage $taskRelative
  $taskDestination = [IO.Path]::GetFullPath((Join-Path $taskTarget $taskRelative))
  if (-not $taskDestination.StartsWith($taskRootPrefix, [StringComparison]::OrdinalIgnoreCase)) { throw "Path outside manuscript: $taskDestination" }
  $taskNewHash = (Get-FileHash -LiteralPath $taskSource -Algorithm SHA256).Hash.ToLowerInvariant()
  if (-not (Test-Path -LiteralPath $taskDestination -PathType Leaf)) { throw "Expected manuscript file disappeared: $taskRelative" }
  $taskOldHash = (Get-FileHash -LiteralPath $taskDestination -Algorithm SHA256).Hash.ToLowerInvariant()
  if ($taskOldHash -eq $taskNewHash) { continue }
  if (-not $taskOriginal.ContainsKey($taskRelative) -or $taskOriginal[$taskRelative] -ne $taskOldHash) {
    throw "Manuscript changed since review: $taskRelative"
  }
  $taskChanges += [pscustomobject]@{ File = $taskRelative; Source = $taskSource; Destination = $taskDestination; SHA256 = $taskNewHash }
}
if (-not $Apply) {
  $taskChanges | Select-Object File
  Write-Output "Ready to update $($taskChanges.Count) files; original contents are unchanged."
  return
}
foreach ($taskChange in $taskChanges) {
  $taskBackup = Join-Path (Join-Path $taskStage 'before_apply') $taskChange.File
  New-Item -ItemType Directory -Path (Split-Path -Parent $taskBackup) -Force | Out-Null
  Copy-Item -LiteralPath $taskChange.Destination -Destination $taskBackup -Force
  Copy-Item -LiteralPath $taskChange.Source -Destination $taskChange.Destination -Force
  if ((Get-FileHash -LiteralPath $taskChange.Destination -Algorithm SHA256).Hash.ToLowerInvariant() -ne $taskChange.SHA256) {
    throw "Copied file failed verification: $($taskChange.File)"
  }
}
$taskChanges | ConvertTo-Json -Depth 3 | Set-Content -LiteralPath (Join-Path $taskStage 'applied_manifest.json')
Write-Output "Updated and verified $($taskChanges.Count) files in $taskTarget."
