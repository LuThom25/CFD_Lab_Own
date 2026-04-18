#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════════
# make_videos.sh  –  erzeugt MP4-Videos aus den Animationsframes
# Aufruf: ./make_videos.sh   (aus dem LidDrivenCavity-Ordner)
# ══════════════════════════════════════════════════════════════════════════════

OUT_DIR="$(dirname "$0")/LidDrivenCavity_Output"
FRM_DIR="$OUT_DIR/frames"
VID_DIR="$OUT_DIR/videos"
mkdir -p "$VID_DIR"

# Frames pro Sekunde → bei 101 Frames und 10 fps = ~10s Video
FPS=10

declare -A TITLES=(
  [u]="u-Velocity (x-Komponente)"
  [v]="v-Velocity (y-Komponente)"
  [p]="Pressure (Druck)"
  [vel]="Velocity Magnitude (Geschwindigkeitsbetrag)"
)

for PREFIX in u v p vel; do
  TITLE="${TITLES[$PREFIX]}"
  OUT="$VID_DIR/video_${PREFIX}.mp4"
  echo "Erstelle: video_${PREFIX}.mp4  ($TITLE)"

  ffmpeg -y \
    -framerate $FPS \
    -i "$FRM_DIR/${PREFIX}_%04d.png" \
    -vf "scale=800:700:flags=lanczos,format=yuv420p" \
    -c:v libx264 \
    -preset slow \
    -crf 18 \
    "$OUT" 2>/dev/null

  if [ $? -eq 0 ]; then
    echo "  ✓  $OUT"
  else
    echo "  ✗  Fehler bei $PREFIX"
  fi
done

# ── Bonus: 2×2 Panel-Video (alle 4 Größen gleichzeitig) ──────────────────────
echo ""
echo "Erstelle: video_panel_2x2.mp4  (alle 4 Größen gleichzeitig)"

ffmpeg -y \
  -framerate $FPS \
  -i "$FRM_DIR/u_%04d.png" \
  -framerate $FPS \
  -i "$FRM_DIR/v_%04d.png" \
  -framerate $FPS \
  -i "$FRM_DIR/p_%04d.png" \
  -framerate $FPS \
  -i "$FRM_DIR/vel_%04d.png" \
  -filter_complex "
    [0:v]scale=600:500[tl];
    [1:v]scale=600:500[tr];
    [2:v]scale=600:500[bl];
    [3:v]scale=600:500[br];
    [tl][tr]hstack[top];
    [bl][br]hstack[bot];
    [top][bot]vstack[out]
  " \
  -map "[out]" \
  -vf "format=yuv420p" \
  -c:v libx264 \
  -preset slow \
  -crf 18 \
  "$VID_DIR/video_panel_2x2.mp4" 2>/dev/null

if [ $? -eq 0 ]; then
  echo "  ✓  $VID_DIR/video_panel_2x2.mp4"
else
  echo "  ✗  Fehler beim Panel-Video"
fi

echo ""
echo "━━━ Alle Videos in: $VID_DIR ━━━"
ls -lh "$VID_DIR"
echo ""
echo "Öffnen mit:"
echo "  open \"$VID_DIR\""
