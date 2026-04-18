#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# make_videos.sh  -  assemble MP4 videos from animation frames
# Run from the LidDrivenCavity directory: bash make_videos.sh
# ══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="$SCRIPT_DIR/LidDrivenCavity_Output"
FRM_DIR="$OUT_DIR/frames"
VID_DIR="$OUT_DIR/videos"

# Remove old videos so nothing from a previous run survives
# Use both rm -rf on the directory AND a glob delete to catch any macOS-synced
# duplicates with spaces in their names (e.g. "video_u 3.mp4").
if [ -d "$VID_DIR" ]; then
    rm -rf "$VID_DIR"
    echo "Removed old videos/"
fi
mkdir -p "$VID_DIR"
# Extra safety: remove any stray .mp4 files macOS may have re-synced back
find "$VID_DIR" -name "*.mp4" -delete 2>/dev/null

FPS=10   # frames per second -> 101 frames = ~10 s video

echo "=== Creating individual videos ==="

for PREFIX in u v p vel; do
    case "$PREFIX" in
        u)   TITLE="u-Velocity (x-component)" ;;
        v)   TITLE="v-Velocity (y-component)" ;;
        p)   TITLE="Pressure (relative, per-frame rescaled)" ;;
        vel) TITLE="Velocity Magnitude |u|" ;;
    esac

    OUT="$VID_DIR/video_${PREFIX}.mp4"
    echo "  Creating video_${PREFIX}.mp4  ($TITLE)"

    ffmpeg -y \
        -framerate $FPS \
        -i "$FRM_DIR/${PREFIX}_%04d.png" \
        -vf "scale=800:700:flags=lanczos,format=yuv420p" \
        -c:v libx264 -preset slow -crf 18 \
        "$OUT" 2>/dev/null

    [ $? -eq 0 ] && echo "    OK: $OUT" || echo "    FAILED: $PREFIX"
done

# ── 2x2 panel video (all four quantities simultaneously) ──────────────────────
echo ""
echo "=== Creating 2x2 panel video ==="

ffmpeg -y \
    -i "$VID_DIR/video_u.mp4" \
    -i "$VID_DIR/video_v.mp4" \
    -i "$VID_DIR/video_p.mp4" \
    -i "$VID_DIR/video_vel.mp4" \
    -filter_complex \
        "[0:v][1:v]hstack[top];
         [2:v][3:v]hstack[bot];
         [top][bot]vstack,format=yuv420p[out]" \
    -map "[out]" \
    -c:v libx264 -preset slow -crf 18 \
    "$VID_DIR/video_panel_2x2.mp4" 2>/dev/null

[ $? -eq 0 ] && echo "    OK: $VID_DIR/video_panel_2x2.mp4" || echo "    FAILED: panel video"

echo ""
echo "=== All videos in: $VID_DIR ==="
ls -lh "$VID_DIR"
echo ""
echo "Open folder with:"
echo "  open \"$VID_DIR\""
