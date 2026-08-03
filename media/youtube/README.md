# Spatial Methods Workbench video package

This folder is the editable production package for two narrated courses.

## Public usage tutorial

Published video: <https://www.youtube.com/watch?v=Ft9groM2ilI>

- `output/SpatialMethodsWorkbench-Public-Usage-v0.2.1.mp4`: narrated 1080p video intended for YouTube
- `output/SpatialMethodsWorkbench-Public-Usage-v0.2.1.pptx`: editable 25-slide presentation with embedded narration
- `output/SpatialMethodsWorkbench-Public-Usage-v0.2.1.srt`: uploadable captions
- `output/SpatialMethodsWorkbench-Public-Usage-v0.2.1-chapters.txt`: YouTube chapters
- `output/SpatialMethodsWorkbench-Public-Usage-v0.2.1-thumbnail.png`: thumbnail source
- `PUBLIC_USAGE_VIDEO_GUIDE.md`: written companion and operational checklist
- `YOUTUBE_UPLOAD.md`: upload title, description, links, and tags

The public tutorial explains workflow selection, input requirements, controls, installation, and reproducibility. It does **not** interpret a user's scientific results.

## Local interpretation training

- `output/SpatialMethodsWorkbench-Interpretation-Training-v0.2.1.mp4`: narrated 1080p local course
- `output/SpatialMethodsWorkbench-Interpretation-Training-v0.2.1.pptx`: editable 29-slide presentation with embedded narration
- `output/SpatialMethodsWorkbench-Interpretation-Training-v0.2.1.srt`: captions
- `output/SpatialMethodsWorkbench-Interpretation-Training-v0.2.1-chapters.txt`: chapters
- `METHOD_INTERPRETATION_GUIDE.md`: detailed written interpretation guide
- `LITERATURE_CONTEXT.md`: publication-to-method evidence map

The interpretation course is deliberately separate from the public application and should remain local unless Sagnik Bhadury explicitly decides otherwise.

## Rebuild

Run from the repository root in Windows PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\media\youtube\build_videos.ps1 -VideoSet All
```

Edit `public_usage_slides.json` or `interpretation_slides.json` to change slide text or narration. The script regenerates narration WAV files, editable PowerPoints, captions, chapters, thumbnails, and MP4s. PowerPoint desktop is required for rendering and video export.

The screenshots and synthetic bundle contain demonstration data only. Do not add non-public or confidential material to this package.
