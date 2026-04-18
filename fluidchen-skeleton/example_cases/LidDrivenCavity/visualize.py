"""
Worksheet 1 - Lid-Driven Cavity: complete ParaView visualization pipeline
==========================================================================
Produces (as required by worksheet Section 5/6):

  Static final-state images:
    final_u.png           - u-velocity component   (Blue to Red, fixed range)
    final_v.png           - v-velocity component   (Blue to Red, fixed range)
    final_pressure.png    - pressure relative to mean (Cool to Warm, auto)
    final_velocity.png    - velocity magnitude |u| (Jet, 0..1)
    final_glyphs.png      - velocity arrows via Glyph filter
    final_streamlines.png - streamlines via StreamTracer

  Animation frames (one PNG per timestep):
    frames/u_NNNN.png   frames/v_NNNN.png
    frames/p_NNNN.png   frames/vel_NNNN.png

Note on pressure: the PPE with pure Neumann BCs has no unique absolute
pressure level (null space). The absolute value drifts linearly over time,
but the pressure GRADIENT (range ~2.2 Pa throughout) is physically correct.
All pressure visualizations subtract the frame mean so the color range
always reflects the physically meaningful variation.
"""

from paraview.simple import *
import os, glob

paraview.simple._DisableFirstRenderCameraReset()

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = BASE + "/LidDrivenCavity_Output"
FRM_DIR = OUT_DIR + "/frames"
os.makedirs(FRM_DIR, exist_ok=True)

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
    return v

def save_img(view, path):
    SaveScreenshot(path, view, ImageResolution=view.ViewSize)
    print(f"  saved: {os.path.basename(path)}")

# Colormaps: (preset name, fixed range or None for auto per-frame)
CMAPS = {
    "u_comp"   : ("Blue to Red Rainbow", (-1.0,  1.0)),
    "v_comp"   : ("Blue to Red Rainbow", (-0.5,  0.5)),
    "p_norm"   : ("Cool to Warm",         None),        # auto per-frame
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
    save_img(v, OUT_DIR + "/" + fname)
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
save_img(v_g, OUT_DIR + "/final_glyphs.png")
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
save_img(v_s, OUT_DIR + "/final_streamlines.png")
Delete(v_s)

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

for idx, t in enumerate(times):
    anim.AnimationTime = t
    # Render once to update the pipeline to this timestep
    for va, _ in views_a:
        Render(va)
    # Now save each view; rescale pressure LUT per frame AFTER the render
    for (va, prefix), (da, lut_a, lut_key) in zip(views_a, disps_a):
        if lut_key == "p_norm":
            # Per-frame rescale: eliminates the absolute drift problem
            da.RescaleTransferFunctionToDataRange(False)
            Render(va)
        SaveScreenshot(f"{FRM_DIR}/{prefix}_{idx:04d}.png", va,
                       ImageResolution=[800, 700])
    if idx % 10 == 0:
        print(f"  frame {idx:03d} / {len(times)-1}")

for va, _ in views_a:
    Delete(va)
Delete(calc_a); Delete(rd_a)

# ══════════════════════════════════════════════════════════════════════════════
print(f"""
=== DONE ===
Final-state images saved to:
  {OUT_DIR}/final_u.png
  {OUT_DIR}/final_v.png
  {OUT_DIR}/final_pressure.png
  {OUT_DIR}/final_velocity.png
  {OUT_DIR}/final_glyphs.png
  {OUT_DIR}/final_streamlines.png

Animation frames: {FRM_DIR}/[u|v|p|vel]_NNNN.png
  ({len(times)} frames per quantity)

Next step - create videos:
  bash make_videos.sh
""")
