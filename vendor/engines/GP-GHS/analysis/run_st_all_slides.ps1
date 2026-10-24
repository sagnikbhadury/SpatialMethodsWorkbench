# Fits GP-GHS to every slide of the Visium data with the gene panel of the
# index slide (S5_Rec_A121573_Rep2, already fitted in st_results/fit_m2).
# One slide at a time; st_fit.R resumes from its task files if interrupted.
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $root
$R = "C:\Users\bhadury\AppData\Local\Programs\R\R-4.4.1\bin\Rscript.exe"
$slides = "S5_Rec_A121573_Rep1","S1_Cec_A551763_Rep1","S1_Cec_A551763_Rep2","S2_ColR_A595688_Rep1","S2_ColR_A595688_Rep2",
          "S3_ColR_A416371_Rep1","S3_ColR_A416371_Rep2","S4_ColSig_A120838_Rep1","S4_ColSig_A120838_Rep2",
          "S6_Rec_A938797_Rep1","S6_Rec_A938797_Rep2","S7_RecSig_A798015_Rep1","S7_RecSig_A798015_Rep2"
$env:ST_N_CORES = "18"
foreach ($s in $slides) {
  $env:ST_SLIDE = $s
  & $R st_prepare_panel.R
  $dir = "st_results\slides\$s"
  if (Test-Path -LiteralPath "$dir\st_prepared.rds") {
    $env:ST_OUT_DIR = "st_results/slides/$s"
    & $R st_fit.R 2> "$dir\fit_run.err.log" | Tee-Object -FilePath "$dir\fit_run.log"
  }
  "=== finished $s at $(Get-Date -Format 'HH:mm:ss')"
}
"ALL SLIDES DONE"
