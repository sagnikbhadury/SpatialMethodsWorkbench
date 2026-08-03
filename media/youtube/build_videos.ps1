param(
  [ValidateSet('Public','Interpretation','All')]
  [string]$VideoSet = 'All',
  [switch]$NoVideo
)

$ErrorActionPreference = 'Stop'
$mediaRoot = $PSScriptRoot
$outputRoot = Join-Path $mediaRoot 'output'
$audioRoot = Join-Path $mediaRoot 'audio'
New-Item -ItemType Directory -Force -Path $outputRoot, $audioRoot | Out-Null

function Get-OfficeRgb([int]$r, [int]$g, [int]$b) {
  return $r + (256 * $g) + (65536 * $b)
}

function Get-WavDuration([string]$path) {
  $stream = [System.IO.File]::OpenRead($path)
  $reader = New-Object System.IO.BinaryReader($stream)
  try {
    $riff = -join ([char[]]$reader.ReadBytes(4))
    [void]$reader.ReadUInt32()
    $wave = -join ([char[]]$reader.ReadBytes(4))
    if ($riff -ne 'RIFF' -or $wave -ne 'WAVE') { throw "Not a PCM WAV file: $path" }
    $byteRate = 0
    $dataSize = 0
    while ($stream.Position -lt $stream.Length) {
      $chunkId = -join ([char[]]$reader.ReadBytes(4))
      if ($chunkId.Length -lt 4) { break }
      $chunkSize = $reader.ReadUInt32()
      $chunkStart = $stream.Position
      if ($chunkId -eq 'fmt ') {
        [void]$reader.ReadUInt16()
        [void]$reader.ReadUInt16()
        [void]$reader.ReadUInt32()
        $byteRate = $reader.ReadUInt32()
      } elseif ($chunkId -eq 'data') {
        $dataSize = $chunkSize
      }
      $stream.Position = $chunkStart + $chunkSize + ($chunkSize % 2)
    }
    if ($byteRate -le 0 -or $dataSize -le 0) { throw "WAV duration metadata missing: $path" }
    return [double]$dataSize / [double]$byteRate
  } finally {
    $reader.Close()
    $stream.Close()
  }
}

function New-Narration([string]$text, [string]$path) {
  if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
  $voice = New-Object -ComObject SAPI.SpVoice
  $stream = New-Object -ComObject SAPI.SpFileStream
  try {
    $david = @($voice.GetVoices()) | Where-Object { $_.GetDescription() -match 'David' } | Select-Object -First 1
    if ($null -ne $david) { $voice.Voice = $david }
    $voice.Rate = -1
    $voice.Volume = 100
    $stream.Open($path, 3, $false)
    $voice.AudioOutputStream = $stream
    [void]$voice.Speak($text)
  } finally {
    try { $stream.Close() } catch {}
    [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($stream)
    [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($voice)
  }
}

function Format-SrtTime([double]$seconds) {
  $span = [TimeSpan]::FromSeconds($seconds)
  return ('{0:00}:{1:00}:{2:00},{3:000}' -f [math]::Floor($span.TotalHours), $span.Minutes, $span.Seconds, $span.Milliseconds)
}

function Write-Srt($slides, [double[]]$durations, [string]$path) {
  $items = New-Object System.Collections.Generic.List[string]
  $counter = 1
  $timeline = 0.0
  for ($i = 0; $i -lt $slides.Count; $i++) {
    $sentences = @([regex]::Split($slides[$i].narration.Trim(), '(?<=[.!?])\s+') | Where-Object { $_.Trim().Length -gt 0 })
    $available = [math]::Max(1.0, $durations[$i] - 0.35)
    $weights = @($sentences | ForEach-Object { [math]::Max(10, $_.Length) })
    $weightTotal = ($weights | Measure-Object -Sum).Sum
    $local = 0.0
    for ($j = 0; $j -lt $sentences.Count; $j++) {
      $segment = $available * $weights[$j] / $weightTotal
      $start = $timeline + $local
      $end = [math]::Min($timeline + $available, $start + $segment)
      $items.Add([string]$counter)
      $items.Add("$(Format-SrtTime $start) --> $(Format-SrtTime $end)")
      $items.Add($sentences[$j].Trim())
      $items.Add('')
      $counter++
      $local += $segment
    }
    $timeline += $durations[$i]
  }
  [System.IO.File]::WriteAllLines($path, $items, (New-Object System.Text.UTF8Encoding($false)))
}

function Write-Chapters($slides, [double[]]$durations, [string]$path) {
  $items = New-Object System.Collections.Generic.List[string]
  $timeline = 0.0
  for ($i = 0; $i -lt $slides.Count; $i++) {
    $span = [TimeSpan]::FromSeconds([math]::Floor($timeline))
    $stamp = if ($span.TotalHours -ge 1) {
      '{0}:{1:00}:{2:00}' -f [math]::Floor($span.TotalHours), $span.Minutes, $span.Seconds
    } else {
      '{0}:{1:00}' -f $span.Minutes, $span.Seconds
    }
    $items.Add("$stamp $($slides[$i].title)")
    $timeline += $durations[$i]
  }
  [System.IO.File]::WriteAllLines($path, $items, (New-Object System.Text.UTF8Encoding($false)))
}

function Add-TextBox($slide, [string]$text, [single]$left, [single]$top, [single]$width, [single]$height, [single]$fontSize, [int]$color, [bool]$bold = $false) {
  $shape = $slide.Shapes.AddTextbox(1, $left, $top, $width, $height)
  $shape.TextFrame2.MarginLeft = 0
  $shape.TextFrame2.MarginRight = 0
  $shape.TextFrame2.MarginTop = 0
  $shape.TextFrame2.MarginBottom = 0
  $shape.TextFrame2.WordWrap = -1
  $shape.TextFrame2.TextRange.Text = $text
  $shape.TextFrame2.TextRange.Font.Name = 'Aptos'
  $shape.TextFrame2.TextRange.Font.Size = $fontSize
  $shape.TextFrame2.TextRange.Font.Fill.ForeColor.RGB = $color
  $shape.TextFrame2.TextRange.Font.Bold = $(if ($bold) { -1 } else { 0 })
  return $shape
}

function Build-VideoSet([string]$name, [string]$jsonFile, [string]$baseName) {
  Write-Host "Loading $name slide specification..."
  # Windows PowerShell 5.1 can preserve the JSON top-level array as one
  # pipeline object when conversion is nested directly inside @(...).
  # Parse first, then materialize its elements so each slide is indexed.
  $parsedSlides = Get-Content -LiteralPath (Join-Path $mediaRoot $jsonFile) -Raw -Encoding UTF8 | ConvertFrom-Json
  $slides = @($parsedSlides | ForEach-Object { $_ })
  $setAudio = Join-Path $audioRoot $baseName
  New-Item -ItemType Directory -Force -Path $setAudio | Out-Null
  $durations = New-Object double[] $slides.Count

  for ($i = 0; $i -lt $slides.Count; $i++) {
    $wav = Join-Path $setAudio ('slide-{0:00}.wav' -f ($i + 1))
    Write-Host ("Narrating slide {0}/{1}: {2}" -f ($i + 1), $slides.Count, $slides[$i].title)
    New-Narration $slides[$i].narration $wav
    $durations[$i] = (Get-WavDuration $wav) + 0.8
  }

  $srtPath = Join-Path $outputRoot ($baseName + '.srt')
  Write-Srt $slides $durations $srtPath
  $chaptersPath = Join-Path $outputRoot ($baseName + '-chapters.txt')
  Write-Chapters $slides $durations $chaptersPath

  $powerPoint = New-Object -ComObject PowerPoint.Application
  $presentation = $null
  try {
    $powerPoint.Visible = -1
    $presentation = $powerPoint.Presentations.Add()
    $presentation.PageSetup.SlideWidth = 960
    $presentation.PageSetup.SlideHeight = 540
    $presentation.SlideShowSettings.AdvanceMode = 2

    $ink = Get-OfficeRgb 28 40 48
    $soft = Get-OfficeRgb 78 91 98
    $cream = Get-OfficeRgb 255 253 248
    $teal = Get-OfficeRgb 23 107 104
    $coral = Get-OfficeRgb 212 98 69
    $white = Get-OfficeRgb 255 255 255

    for ($i = 0; $i -lt $slides.Count; $i++) {
      $spec = $slides[$i]
      $slide = $presentation.Slides.Add($i + 1, 12)
      $slide.FollowMasterBackground = 0
      $slide.Background.Fill.Solid()
      $slide.Background.Fill.ForeColor.RGB = $cream

      $bar = $slide.Shapes.AddShape(1, 0, 0, 14, 540)
      $bar.Fill.ForeColor.RGB = $(if (($i % 2) -eq 0) { $teal } else { $coral })
      $bar.Line.Visible = 0

      [void](Add-TextBox $slide $spec.kicker 48 26 860 20 11 $coral $true)
      $titleSize = $(if ($spec.title.Length -gt 52) { 22 } elseif ($spec.title.Length -gt 40) { 27 } else { 34 })
      [void](Add-TextBox $slide $spec.title 48 52 860 62 $titleSize $ink $true)

      $hasImage = $null -ne $spec.PSObject.Properties['image'] -and
        -not [string]::IsNullOrWhiteSpace([string]$spec.image)
      if ($hasImage) {
        $imagePath = Join-Path $mediaRoot $spec.image
        $panel = $slide.Shapes.AddShape(5, 432, 126, 480, 318)
        $panel.Fill.ForeColor.RGB = $white
        $panel.Line.ForeColor.RGB = Get-OfficeRgb 220 220 216
        $panel.Line.Weight = 1
        [void]$slide.Shapes.AddPicture($imagePath, 0, -1, 440, 134, 464, 261)
        $bodyFont = $(if ($spec.body.Length -gt 340) { 16 } else { 18 })
        [void](Add-TextBox $slide $spec.body 48 134 345 315 $bodyFont $soft $false)
      } else {
        $bodyFont = $(if ($spec.body.Length -gt 430) { 19 } elseif ($spec.body.Length -gt 280) { 22 } else { 25 })
        [void](Add-TextBox $slide $spec.body 58 145 830 300 $bodyFont $soft $false)
      }

      $footer = "Spatial Methods Workbench v0.2.1  |  doi.org/10.5281/zenodo.21764196  |  $($i + 1)/$($slides.Count)"
      [void](Add-TextBox $slide $footer 48 505 860 16 9 $soft $false)

      $wav = Join-Path $setAudio ('slide-{0:00}.wav' -f ($i + 1))
      $media = $slide.Shapes.AddMediaObject2($wav, 0, -1, 930, 520, 1, 1)
      $media.AnimationSettings.PlaySettings.PlayOnEntry = -1
      $media.AnimationSettings.PlaySettings.HideWhileNotPlaying = -1
      $slide.SlideShowTransition.AdvanceOnClick = 0
      $slide.SlideShowTransition.AdvanceOnTime = -1
      $slide.SlideShowTransition.AdvanceTime = [single]$durations[$i]
    }

    $pptxPath = Join-Path $outputRoot ($baseName + '.pptx')
    if (Test-Path -LiteralPath $pptxPath) { Remove-Item -LiteralPath $pptxPath -Force }
    $presentation.SaveAs($pptxPath, 24)

    $thumbnailPath = Join-Path $outputRoot ($baseName + '-thumbnail.png')
    if (Test-Path -LiteralPath $thumbnailPath) { Remove-Item -LiteralPath $thumbnailPath -Force }
    $presentation.Slides.Item(1).Export($thumbnailPath, 'PNG', 1280, 720)

    if (-not $NoVideo) {
      $mp4Path = Join-Path $outputRoot ($baseName + '.mp4')
      if (Test-Path -LiteralPath $mp4Path) { Remove-Item -LiteralPath $mp4Path -Force }
      Write-Host "Exporting $name MP4. This may take several minutes..."
      $presentation.CreateVideo($mp4Path, $true, 5, 1080, 30, 85)
      $deadline = (Get-Date).AddMinutes(45)
      do {
        Start-Sleep -Seconds 5
        $status = [int]$presentation.CreateVideoStatus
        Write-Host "CreateVideoStatus=$status"
        if ((Get-Date) -gt $deadline) { throw "PowerPoint video export timed out for $name." }
      } while ($status -in 0,1,2)
      if ($status -ne 3 -or -not (Test-Path -LiteralPath $mp4Path)) { throw "PowerPoint video export failed for $name (status $status)." }
    }
  } finally {
    if ($null -ne $presentation) {
      try { $presentation.Close() } catch {}
      [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($presentation)
    }
    try { $powerPoint.Quit() } catch {}
    [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($powerPoint)
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
  }
}

if ($VideoSet -in @('Public','All')) {
  Build-VideoSet 'Public usage tutorial' 'public_usage_slides.json' 'SpatialMethodsWorkbench-Public-Usage-v0.2.1'
}
if ($VideoSet -in @('Interpretation','All')) {
  Build-VideoSet 'Interpretation training' 'interpretation_slides.json' 'SpatialMethodsWorkbench-Interpretation-Training-v0.2.1'
}

Get-ChildItem -LiteralPath $outputRoot | Select-Object Name,Length,LastWriteTime
