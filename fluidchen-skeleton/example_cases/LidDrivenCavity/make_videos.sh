#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# make_videos.sh  -  assemble MP4 videos from animation frames
# Run from the LidDrivenCavity directory: bash make_videos.sh
# ══════════════════════════════════════════════════════════════════════════════

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="$SCRIPT_DIR/LidDrivenCavity_Output"
TASK4_DIR="$OUT_DIR/task4"
FRM_DIR="$TASK4_DIR/frames"
VID_DIR="$TASK4_DIR/videos"

# Remove old videos so nothing from a previous run survives
if [ -d "$VID_DIR" ]; then
    rm -rf "$VID_DIR"
    echo "Removed old task4/videos/"
fi
mkdir -p "$VID_DIR"
# Extra safety: remove any stray .mp4 files macOS may have re-synced back
find "$VID_DIR" -name "*.mp4" -delete 2>/dev/null

FPS=10   # frames per second -> 101 frames = ~10 s video

echo "=== Creating individual videos ==="

for PREFIX in u v p vel vec vec_clean; do
    case "$PREFIX" in
        u)         TITLE="u-Velocity (x-component)" ;;
        v)         TITLE="v-Velocity (y-component)" ;;
        p)         TITLE="Pressure (relative, per-frame rescaled)" ;;
        vel)       TITLE="Velocity Magnitude |u|" ;;
        vec)       TITLE="Velocity Vectors (Jet-coloured by magnitude)" ;;
        vec_clean) TITLE="Velocity Vectors (direction only, white arrows)" ;;
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

# ── 2x3 panel video (all 6 unique quantities) ────────────────────────────────
# Row 1: u-velocity | v-velocity | pressure
# Row 2: velocity magnitude | Jet-coloured vectors | direction-only vectors
echo ""
echo "=== Creating 2x3 panel video (u, v, p  |  vel, vec, vec_clean) ==="

ffmpeg -y \
    -i "$VID_DIR/video_u.mp4" \
    -i "$VID_DIR/video_v.mp4" \
    -i "$VID_DIR/video_p.mp4" \
    -i "$VID_DIR/video_vel.mp4" \
    -i "$VID_DIR/video_vec.mp4" \
    -i "$VID_DIR/video_vec_clean.mp4" \
    -filter_complex \
        "[0:v][1:v][2:v]hstack=inputs=3[top];
         [3:v][4:v][5:v]hstack=inputs=3[bot];
         [top][bot]vstack,format=yuv420p[out]" \
    -map "[out]" \
    -c:v libx264 -preset slow -crf 18 \
    "$VID_DIR/video_panel_2x3.mp4" 2>/dev/null

[ $? -eq 0 ] && echo "    OK: $VID_DIR/video_panel_2x3.mp4" || echo "    FAILED: 2x3 panel video"

echo ""
echo "=== All videos in: $VID_DIR ==="
ls -lh "$VID_DIR"
echo ""
echo "Open folder with:"
echo "  open \"$VID_DIR\""
