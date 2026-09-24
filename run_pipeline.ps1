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
& $py tests/test_outliers.py
if ($LASTEXITCODE -ne 0) { Write-Error "outlier bounds failed their tests. Nothing was run." }
& $py tests/test_waterfall_cache.py
if ($LASTEXITCODE -ne 0) { Write-Error "the waterfall cache failed its tests. Nothing was run." }
& $py tests/test_doses.py
if ($LASTEXITCODE -ne 0) { Write-Error "dose unit conversion failed its tests. Nothing was run." }
& $py tests/test_pooling.py
if ($LASTEXITCODE -ne 0) { Write-Error "the federated-pooling exports are not poolable" }
& $py tests/test_paths.py
if ($LASTEXITCODE -ne 0) { Write-Error "output locations failed their tests. Nothing was run." }
& $py tests/test_build_cohort.py
if ($LASTEXITCODE -ne 0) { Write-Error "Phase 0 logic failed its tests. Nothing was run." }
Write-Host ""

Write-Host "== Phase 0: build cohort (Python) =="
& $py code/01_build_cohort.py
if ($LASTEXITCODE -ne 0) { Write-Error "Phase 0 failed" }

# Steps 3-6 are TABLED as of 2026-09-23 and deliberately excluded; 07_outcomes is
# an unconditional stop(). Run a tabled script by hand if you need it.
$steps = @("02_descriptive_cohort","03_delivery_states","04_landmark_cohort","05_exemplar")
foreach ($s in $steps) {
  Write-Host "== $s (R) =="
  & Rscript "code/$s.R"
  if ($LASTEXITCODE -ne 0) { Write-Error "$s failed" }
}

Write-Host "Pipeline complete. Review output/final_no_phi/ before sharing."
