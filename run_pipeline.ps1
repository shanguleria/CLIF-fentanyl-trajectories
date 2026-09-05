# Run the full pipeline. Stops on first error.
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

if (-not (Test-Path "config/config.json")) {
  Write-Error "config/config.json not found. Run: Copy-Item config/config_template.json config/config.json, then edit it for your site."
}

$py = "python"
if (Test-Path ".venv/Scripts/python.exe") { $py = ".venv/Scripts/python.exe" }

# Preflight. Fails here, in a second, rather than an hour into a table read.
Write-Host "== Preflight =="
& $py code/check_config.py
if ($LASTEXITCODE -ne 0) { Write-Error "config is not usable. Nothing was run." }
& $py tests/test_covariates.py
if ($LASTEXITCODE -ne 0) { Write-Error "config/covariates.json failed its integrity checks. Nothing was run." }
& $py tests/test_fio2.py
if ($LASTEXITCODE -ne 0) { Write-Error "fio2 unit handling failed its tests. Nothing was run." }
Write-Host ""

Write-Host "== Phase 0: build cohort (Python) =="
& $py code/01_build_cohort.py
if ($LASTEXITCODE -ne 0) { Write-Error "Phase 0 failed" }

$steps = @("02_descriptive_trajectory","03_landmark_cohort","04_gbmt_classes",
           "05_lcmm_classes","06_two_indicator","07_outcomes")
foreach ($s in $steps) {
  Write-Host "== $s (R) =="
  & Rscript "code/$s.R"
  if ($LASTEXITCODE -ne 0) { Write-Error "$s failed" }
}

Write-Host "Pipeline complete. Review output/final_no_phi/ before sharing."
