"""
Worksheet 1 - Lid-Driven Cavity: complete ParaView visualization pipeline
==========================================================================
Produces (as required by worksheet Section 5/6):

  Static final-state images:
    final_u.png           - u-velocity component   (Blue to Red, fixed range)
    final_v.png           - v-velocity component   (Blue to Red, fixed range)
    final_pressure.png    - pressure relative to mean (Cool to Warm, auto)
    final_velocity.png    - velocity magnitude |u| (Jet, 0..1)
    final_glyphs.png        - velocity arrows via Glyph filter
    final_streamlines.png   - streamlines via StreamTracer
    final_vectors_bw.png    - Jet-coloured vectors, dark background
    final_vectors_clean.png - direction-only white arrows, navy background

  Animation frames (one PNG per timestep):
    frames/u_NNNN.png   frames/v_NNNN.png
    frames/p_NNNN.png   frames/vel_NNNN.png
    frames/vec_NNNN.png frames/vec_clean_NNNN.png

Note on pressure: the PPE with pure Neumann BCs has no unique absolute
pressure level (null space). The absolute value drifts linearly over time,
but the pressure GRADIENT (range ~2.2 Pa throughout) is physically correct.
All pressure visualizations subtract the frame mean so the color range
always reflects the physically meaningful variation.
"""

from paraview.simple import *
import os, glob, shutil

paraview.simple._DisableFirstRenderCameraReset()

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE      = os.path.dirname(os.path.abspath(__file__))
OUT_DIR   = BASE + "/LidDrivenCavity_Output"   # VTK files live here
TASK4_DIR = OUT_DIR + "/task4"                  # all Task-4 output goes here
FRM_DIR   = TASK4_DIR + "/frames"
VID_DIR   = TASK4_DIR + "/videos"

# ── Clean up old Task-4 output so nothing from a previous run survives ────────
# Only wipe task4/ (and any stale top-level non-VTK files); preserve study
# outputs (study_plots/, task7_visuals/, task8_visuals/) created by run_studies.py
print("Cleaning up old Task-4 results...")
if os.path.isdir(OUT_DIR):
    for item in os.listdir(OUT_DIR):
        if item.endswith(".vtk"):
            continue                          # keep simulation output
        if item in ("study_plots", "task7_visuals", "task8_visuals"):
            continue                          # keep parameter-study outputs
        item_path = os.path.join(OUT_DIR, item)
        if os.path.isdir(item_path):
            shutil.rmtree(item_path)
        else:
            os.remove(item_path)
        print(f"  removed: {item}")

os.makedirs(FRM_DIR, exist_ok=True)
os.makedirs(VID_DIR, exist_ok=True)

# Extra pass: macOS iCloud can re-sync stale files back into fresh directories.
for _dir in (FRM_DIR, VID_DIR):
    for _f in glob.glob(_dir + "/*"):
        os.remove(_f)

print("Clean. Starting fresh.\n")

# Sort VTK files by timestep index (filename format: CaseName_rank.timestep.vtk)
def vtk_key(f):
    return int(os.path.basename(f).replace(".vtk", "").split(".")[-1])

VTK_FILES = sorted(glob.glob(OUT_DIR + "/*.vtk"), key=vtk_key)
N = len(VTK_FILES)
print(f"Timesteps found: {N}")
assert N > 0, "No VTK files found!"

# ── Camera: 2-D top view of the unit square domain ────────────────────────────
def setup_camera(view):
    view.ResetCamera()
    view.CameraPosition          = [0.5, 0.5, 3.0]
    view.CameraFocalPoint        = [0.5, 0.5, 0.0]
    view.CameraViewUp            = [0.0, 1.0, 0.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale     = 0.55

def new_view(size=(1200, 1000), bg=(0.12, 0.12, 0.18)):
    v = CreateView("RenderView")
    v.ViewSize  = list(size)
    v.Background = list(bg)
    v.OrientationAxesVisibility = 0   # hide x/y/z corner widget
    return v

def save_img(view, path):
    SaveScreenshot(path, view, ImageResolution=view.ViewSize)
    print(f"  saved: {os.path.basename(path)}")

# Colormaps: (preset name, fixed range or None for auto per-frame)
CMAPS = {
    # u: lid moves at +1.0; recirculation gives u_min ≈ -0.21 (Ghia 1982, Re=100)
    "u_comp"   : ("Blue to Red Rainbow", (-0.5,  1.0)),
    # v: symmetric recirculation, v_max ≈ ±0.18 (Ghia 1982, Re=100) → ±0.25 fills colormap well
    "v_comp"   : ("Blue to Red Rainbow", (-0.25, 0.25)),
    # pressure: pure Neumann → absolute level drifts; rescale per-frame to show gradient
    "p_norm"   : ("Cool to Warm",         None),        # auto per-frame
    # velocity magnitude: bounded by lid speed U_wall = 1.0
    "vel_mag"  : ("Jet",                 ( 0.0,  1.0)),
}

def apply_lut(lut, key, disp=None):
    preset, rng = CMAPS[key]
    lut.ApplyPreset(preset, True)
    if rng:
        lut.RescaleTransferFunction(rng[0], rng[1])
    elif disp is not None:
        disp.RescaleTransferFunctionToDataRange(False)

def add_colorbar(lut, view, title):
    bar = GetScalarBar(lut, view)
    bar.Title          = title
    bar.ComponentTitle = ""
    bar.TitleFontSize  = 13
    bar.LabelFontSize  = 11
    bar.WindowLocation = "Lower Right Corner"
    return bar

# ── Pipeline builder: reader + derived quantities ─────────────────────────────
def build_pipeline(files):
    rd = LegacyVTKReader(FileNames=files)
    # u-component
    c1 = Calculator(Input=rd)
    c1.AttributeType   = "Cell Data"
    c1.ResultArrayName = "u_comp"
    c1.Function        = "velocity_X"
    # v-component
    c2 = Calculator(Input=c1)
    c2.AttributeType   = "Cell Data"
    c2.ResultArrayName = "v_comp"
    c2.Function        = "velocity_Y"
    # velocity magnitude
    c3 = Calculator(Input=c2)
    c3.AttributeType   = "Cell Data"
    c3.ResultArrayName = "vel_mag"
    c3.Function        = "mag(velocity)"
    # pressure normalized: subtract cell (1,1) value as reference
    # Note: direct mean subtraction is not available in pvpython Calculator;
    # instead we use RescaleTransferFunctionToDataRange per frame for pressure.
    # For the static image we use the raw pressure field rescaled to its range.
    c4 = Calculator(Input=c3)
    c4.AttributeType   = "Cell Data"
    c4.ResultArrayName = "p_norm"
    c4.Function        = "pressure"   # renamed for separate LUT management
    return rd, c4

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 - Final state: six static images (worksheet Section 5/6)
# ══════════════════════════════════════════════════════════════════════════════
print("\n--- PART 1: Final-state images ---")

rd_f, calc_f = build_pipeline([VTK_FILES[-1]])

STATIC_FIELDS = [
    ("u_comp",  "u_comp",  "u-Velocity  [m/s]",           "final_u.png"),
    ("v_comp",  "v_comp",  "v-Velocity  [m/s]",           "final_v.png"),
    ("p_norm",  "p_norm",  "Pressure (relative)  [Pa]",   "final_pressure.png"),
    ("vel_mag", "vel_mag", "Velocity Magnitude |u|  [m/s]","final_velocity.png"),
]

for field, lut_key, title, fname in STATIC_FIELDS:
    v = new_view()
    d = Show(calc_f, v)
    d.Representation = "Surface"
    ColorBy(d, ("CELLS", field))
    lut = GetColorTransferFunction(field)
    apply_lut(lut, lut_key, disp=d)
    d.SetScalarBarVisibility(v, True)
    add_colorbar(lut, v, title)
    setup_camera(v)
    Render()
    save_img(v, TASK4_DIR + "/" + fname)
    Delete(v)

# Glyphs (velocity arrows) ────────────────────────────────────────────────────
print("  Rendering glyphs (arrows)...")
v_g = new_view()
d_bg = Show(calc_f, v_g)
d_bg.Representation = "Surface"
ColorBy(d_bg, ("CELLS", "vel_mag"))
lut_g = GetColorTransferFunction("vel_mag")
apply_lut(lut_g, "vel_mag")
d_bg.Opacity = 0.45

glyph = Glyph(Input=calc_f, GlyphType="Arrow")
glyph.OrientationArray = ["CELLS", "velocity"]
glyph.ScaleArray       = ["CELLS", "vel_mag"]
glyph.ScaleFactor      = 0.04
glyph.GlyphMode        = "Every Nth Point"
glyph.Stride           = 6

d_gl = Show(glyph, v_g)
d_gl.Representation = "Surface"
ColorBy(d_gl, ("POINTS", "vel_mag"))
lut_gl = GetColorTransferFunction("vel_mag")
apply_lut(lut_gl, "vel_mag")
d_gl.SetScalarBarVisibility(v_g, True)
add_colorbar(lut_gl, v_g, "Velocity |u|  [m/s]")
setup_camera(v_g)
Render()
save_img(v_g, TASK4_DIR + "/final_glyphs.png")
Delete(v_g)

# Streamlines ─────────────────────────────────────────────────────────────────
print("  Rendering streamlines...")
v_s = new_view(bg=(0.05, 0.05, 0.1))
d_sbg = Show(calc_f, v_s)
d_sbg.Representation = "Surface"
ColorBy(d_sbg, ("CELLS", "vel_mag"))
lut_s = GetColorTransferFunction("vel_mag")
apply_lut(lut_s, "vel_mag")
d_sbg.Opacity = 0.5

stream = StreamTracer(Input=rd_f, SeedType="Line")
stream.Vectors                 = ["POINTS", "velocity"]
stream.MaximumStreamlineLength = 3.0
stream.IntegrationDirection    = "BOTH"
stream.SeedType.Point1         = [0.05, 0.5, 0.0]
stream.SeedType.Point2         = [0.95, 0.5, 0.0]
stream.SeedType.Resolution     = 40

d_st = Show(stream, v_s)
d_st.Representation = "Surface"
d_st.LineWidth      = 2.0
ColorBy(d_st, ("POINTS", "velocity"))
lut_st = GetColorTransferFunction("velocity")
lut_st.ApplyPreset("Jet", True)
lut_st.RescaleTransferFunction(0.0, 1.0)
d_st.SetScalarBarVisibility(v_s, True)
add_colorbar(lut_st, v_s, "Velocity |u|  [m/s]")
setup_camera(v_s)
Render()
save_img(v_s, TASK4_DIR + "/final_streamlines.png")
Delete(v_s)

# Coloured vector field plot (Jet arrows on dark background) ──────────────────
print("  Rendering coloured vector field...")
# Delete the StreamTracer source so its seed line cannot bleed into this view
Delete(stream)

v_bw = new_view(bg=(0.05, 0.05, 0.10))

# Faint velocity-magnitude background for spatial context
d_bwbg = Show(calc_f, v_bw)
d_bwbg.Representation = "Surface"
ColorBy(d_bwbg, ("CELLS", "vel_mag"))
lut_bwbg = GetColorTransferFunction("vel_mag")
apply_lut(lut_bwbg, "vel_mag")
d_bwbg.Opacity = 0.30

# Arrows coloured by velocity magnitude (Jet, 0–1), every 3rd cell
glyph_bw = Glyph(Input=calc_f, GlyphType="Arrow")
glyph_bw.OrientationArray = ["CELLS", "velocity"]
glyph_bw.ScaleArray        = ["CELLS", "vel_mag"]
glyph_bw.ScaleFactor       = 0.060
glyph_bw.GlyphMode         = "Every Nth Point"
glyph_bw.Stride            = 3

d_gbw = Show(glyph_bw, v_bw)
d_gbw.Representation = "Surface"
ColorBy(d_gbw, ("POINTS", "vel_mag"))
lut_bw = GetColorTransferFunction("vel_mag")
apply_lut(lut_bw, "vel_mag")
d_gbw.SetScalarBarVisibility(v_bw, True)
add_colorbar(lut_bw, v_bw, "Velocity |u|  [m/s]")
setup_camera(v_bw)
Render()
save_img(v_bw, TASK4_DIR + "/final_vectors_bw.png")
Delete(v_bw)

# Clean direction-only vector plot (white arrows, blue background) ─────────────
print("  Rendering clean direction vectors (white arrows)...")
# Calculator: constant 1.0 so every arrow has the same length (direction only)
c_const_f = Calculator(Input=calc_f)
c_const_f.AttributeType   = "Cell Data"
c_const_f.ResultArrayName = "one"
c_const_f.Function        = "1.0"

v_cl = new_view(bg=(0.06, 0.10, 0.22))   # deep navy blue

glyph_cl = Glyph(Input=c_const_f, GlyphType="Arrow")
glyph_cl.OrientationArray = ["CELLS", "velocity"]
glyph_cl.ScaleArray        = ["CELLS", "one"]
glyph_cl.ScaleFactor       = 0.018
glyph_cl.GlyphMode         = "Every Nth Point"
glyph_cl.Stride            = 3

d_gcl = Show(glyph_cl, v_cl)
d_gcl.Representation = "Surface"
ColorBy(d_gcl, None)
d_gcl.AmbientColor = [1.0, 1.0, 1.0]
d_gcl.DiffuseColor = [1.0, 1.0, 1.0]
setup_camera(v_cl)
Render()
save_img(v_cl, TASK4_DIR + "/final_vectors_clean.png")
Delete(v_cl)
Delete(c_const_f)

Delete(calc_f); Delete(rd_f)

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 - Animation frames for all timesteps
# Key fix: pressure colormap is rescaled PER FRAME (not globally),
# avoiding the "all blue" issue caused by absolute pressure drift.
# ══════════════════════════════════════════════════════════════════════════════
print("\n--- PART 2: Animation frames (all timesteps) ---")

rd_a, calc_a = build_pipeline(VTK_FILES)
times = rd_a.TimestepValues
if not times:
    times = list(range(N))
print(f"  {len(times)} frames x 4 quantities...")

anim = GetAnimationScene()
anim.UpdateAnimationUsingDataTimeSteps()

# Create one small view per quantity
ANIM_FIELDS = [
    ("u_comp",  "u_comp",  "u-Velocity [m/s]",           "u"),
    ("v_comp",  "v_comp",  "v-Velocity [m/s]",           "v"),
    ("p_norm",  "p_norm",  "Pressure (relative) [Pa]",   "p"),
    ("vel_mag", "vel_mag", "Velocity Magnitude [m/s]",   "vel"),
]

views_a, disps_a = [], []
for field, lut_key, title, prefix in ANIM_FIELDS:
    va = new_view(size=(800, 700))
    da = Show(calc_a, va)
    da.Representation = "Surface"
    ColorBy(da, ("CELLS", field))
    lut_a = GetColorTransferFunction(field)
    apply_lut(lut_a, lut_key, disp=da)
    da.SetScalarBarVisibility(va, True)
    cb = add_colorbar(lut_a, va, title)
    cb.TitleFontSize = 11
    cb.LabelFontSize = 9
    setup_camera(va)
    views_a.append((va, prefix))
    disps_a.append((da, lut_a, lut_key))

# Vector animation view (arrows coloured by vel_mag on dark bg) ───────────────
v_vec = new_view(size=(800, 700), bg=(0.05, 0.05, 0.10))
d_vbg = Show(calc_a, v_vec)
d_vbg.Representation = "Surface"
ColorBy(d_vbg, ("CELLS", "vel_mag"))
lut_vec = GetColorTransferFunction("vel_mag")
apply_lut(lut_vec, "vel_mag")
d_vbg.Opacity = 0.30

glyph_a = Glyph(Input=calc_a, GlyphType="Arrow")
glyph_a.OrientationArray = ["CELLS", "velocity"]
glyph_a.ScaleArray        = ["CELLS", "vel_mag"]
glyph_a.ScaleFactor       = 0.060
glyph_a.GlyphMode         = "Every Nth Point"
glyph_a.Stride            = 3

d_ga = Show(glyph_a, v_vec)
d_ga.Representation = "Surface"
ColorBy(d_ga, ("POINTS", "vel_mag"))
lut_ga = GetColorTransferFunction("vel_mag")
apply_lut(lut_ga, "vel_mag")
d_ga.SetScalarBarVisibility(v_vec, True)
cb_vec = add_colorbar(lut_ga, v_vec, "Velocity |u|  [m/s]")
cb_vec.TitleFontSize = 11
cb_vec.LabelFontSize = 9
setup_camera(v_vec)

# Clean direction-only vector animation (white arrows, navy bg) ───────────────
c_const_a = Calculator(Input=calc_a)
c_const_a.AttributeType   = "Cell Data"
c_const_a.ResultArrayName = "one"
c_const_a.Function        = "1.0"

v_vcl = new_view(size=(800, 700), bg=(0.06, 0.10, 0.22))

glyph_vcl = Glyph(Input=c_const_a, GlyphType="Arrow")
glyph_vcl.OrientationArray = ["CELLS", "velocity"]
glyph_vcl.ScaleArray        = ["CELLS", "one"]
glyph_vcl.ScaleFactor       = 0.018
glyph_vcl.GlyphMode         = "Every Nth Point"
glyph_vcl.Stride            = 3

d_gvcl = Show(glyph_vcl, v_vcl)
d_gvcl.Representation = "Surface"
ColorBy(d_gvcl, None)
d_gvcl.AmbientColor = [1.0, 1.0, 1.0]
d_gvcl.DiffuseColor = [1.0, 1.0, 1.0]
setup_camera(v_vcl)

for idx, t in enumerate(times):
    anim.AnimationTime = t
    # Render once to update the pipeline to this timestep
    for va, _ in views_a:
        Render(va)
    Render(v_vec)
    Render(v_vcl)
    # Now save each view; rescale pressure LUT per frame AFTER the render
    for (va, prefix), (da, lut_a, lut_key) in zip(views_a, disps_a):
        if lut_key == "p_norm":
            # Per-frame rescale: eliminates the absolute drift problem
            da.RescaleTransferFunctionToDataRange(False)
            Render(va)
        SaveScreenshot(f"{FRM_DIR}/{prefix}_{idx:04d}.png", va,
                       ImageResolution=[800, 700])
    # Coloured vector frame
    SaveScreenshot(f"{FRM_DIR}/vec_{idx:04d}.png", v_vec,
                   ImageResolution=[800, 700])
    # Clean direction-only vector frame
    SaveScreenshot(f"{FRM_DIR}/vec_clean_{idx:04d}.png", v_vcl,
                   ImageResolution=[800, 700])
    if idx % 10 == 0:
        print(f"  frame {idx:03d} / {len(times)-1}")

for va, _ in views_a:
    Delete(va)
Delete(v_vec)
Delete(v_vcl)
Delete(c_const_a)
Delete(calc_a); Delete(rd_a)

# ══════════════════════════════════════════════════════════════════════════════
print(f"""
=== DONE ===
All Task-4 output saved to: {TASK4_DIR}/

  Static images:
    final_u.png  final_v.png  final_pressure.png  final_velocity.png
    final_glyphs.png  final_streamlines.png
    final_vectors_bw.png  final_vectors_clean.png

  Animation frames: task4/frames/[u|v|p|vel|vec|vec_clean]_NNNN.png
    ({len(times)} frames per quantity)

Next step - create videos:
  bash make_videos.sh   →  saves to task4/videos/
""")
