"""
Worksheet 1 – Lid-Driven Cavity: vollständige Visualisierung mit pvpython
=========================================================================
Erzeugt laut Arbeitsblatt Section 5/6:
  ① Endzustand-Bilder (static):
      final_u.png           – u-Komponente (Surface)
      final_v.png           – v-Komponente (Surface)
      final_pressure.png    – Druck p (Surface)
      final_velocity.png    – Geschwindigkeitsbetrag |u| (Surface)
      final_glyphs.png      – Geschwindigkeitspfeile (Glyph-Filter)
      final_streamlines.png – Stromlinien (StreamTracer)
      final_panel.png       – 2×3 Übersichtsbild aller Größen

  ② Animationsframes (je Zeitschritt):
      frames/u_NNNN.png
      frames/v_NNNN.png
      frames/p_NNNN.png
      frames/vel_NNNN.png

  ③ Videos (via ffmpeg, separat aufrufbar):
      video_u.mp4
      video_v.mp4
      video_pressure.mp4
      video_velocity.mp4
"""

from paraview.simple import *
import os, glob

paraview.simple._DisableFirstRenderCameraReset()

# ── Pfade ──────────────────────────────────────────────────────────────────────
BASE      = os.path.dirname(os.path.abspath(__file__))
OUT_DIR   = BASE + "/LidDrivenCavity_Output"
FRM_DIR   = OUT_DIR + "/frames"
os.makedirs(FRM_DIR, exist_ok=True)

# VTK-Dateien sortiert nach Zeitschritt (Format: CaseName_rank.timestep.vtk)
def vtk_key(f):
    return int(os.path.basename(f).replace(".vtk","").split(".")[-1])

VTK_FILES = sorted(glob.glob(OUT_DIR + "/*.vtk"), key=vtk_key)
N = len(VTK_FILES)
print(f"Zeitschritte gefunden: {N}")
assert N > 0, "Keine VTK-Dateien gefunden!"

# ── Kamera-Preset (2-D Draufsicht, Einheitsdomäne 1×1) ───────────────────────
def setup_camera(view):
    view.ResetCamera()
    view.CameraPosition         = [0.5, 0.5, 3.0]
    view.CameraFocalPoint       = [0.5, 0.5, 0.0]
    view.CameraViewUp           = [0.0, 1.0, 0.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale    = 0.55

def new_view(size=(1200, 1000), bg=(0.12, 0.12, 0.18)):
    v = CreateView("RenderView")
    v.ViewSize  = list(size)
    v.Background = list(bg)
    return v

def colorbar(lut, view, title):
    disp_dummy = GetDisplayProperties(view=view)
    bar = GetScalarBar(lut, view)
    bar.Title          = title
    bar.ComponentTitle = ""
    bar.TitleFontSize  = 13
    bar.LabelFontSize  = 11
    bar.WindowLocation = "Lower Right Corner"
    return bar

def save(view, path):
    SaveScreenshot(path, view, ImageResolution=view.ViewSize)
    print(f"  ✓  {os.path.basename(path)}")

# Colormap-Presets pro Größe
CMAPS = {
    "u"        : ("Blue to Red Rainbow", (-1.0,  1.0)),  # u: blau negativ, rot positiv
    "v"        : ("Blue to Red Rainbow", (-0.5,  0.5)),  # v: blau negativ, rot positiv
    "pressure" : ("Cool to Warm",         None),          # p: auto-range
    "vel_mag"  : ("Jet",                 ( 0.0,  1.0)),  # |u|: 0..U_wall
}

def apply_cmap(lut, name, data_range=None):
    cname, fixed_range = CMAPS[name]
    lut.ApplyPreset(cname, True)
    if fixed_range:
        lut.RescaleTransferFunction(fixed_range[0], fixed_range[1])
    elif data_range:
        lut.RescaleTransferFunction(data_range[0], data_range[1])

# ══════════════════════════════════════════════════════════════════════════════
# Funktion: einen Zeitschritt laden + Calculator für alle Größen
# ══════════════════════════════════════════════════════════════════════════════
def load_step(files):
    rd = LegacyVTKReader(FileNames=files)
    # u-Komponente
    cu = Calculator(Input=rd)
    cu.AttributeType   = "Cell Data"
    cu.ResultArrayName = "u_comp"
    cu.Function        = "velocity_X"
    # v-Komponente
    cv = Calculator(Input=cu)
    cv.AttributeType   = "Cell Data"
    cv.ResultArrayName = "v_comp"
    cv.Function        = "velocity_Y"
    # Geschwindigkeitsbetrag
    cm = Calculator(Input=cv)
    cm.AttributeType   = "Cell Data"
    cm.ResultArrayName = "vel_mag"
    cm.Function        = "mag(velocity)"
    return rd, cm   # reader + Endfilter

# ══════════════════════════════════════════════════════════════════════════════
# TEIL 1 – Endzustand: einzelne Surface-Plots (worksheet-konform)
# ══════════════════════════════════════════════════════════════════════════════
print("\n━━━ TEIL 1: Endzustand-Bilder ━━━")

rd_final, calc_final = load_step([VTK_FILES[-1]])

for field, array, title, fname in [
    ("u_comp",   "CELLS", "u-Velocity  [m/s]",          "final_u.png"),
    ("v_comp",   "CELLS", "v-Velocity  [m/s]",          "final_v.png"),
    ("pressure", "CELLS", "Pressure  [Pa]",              "final_pressure.png"),
    ("vel_mag",  "CELLS", "Velocity Magnitude |u| [m/s]","final_velocity.png"),
]:
    v = new_view()
    d = Show(calc_final, v)
    d.Representation = "Surface"
    ColorBy(d, ("CELLS", field))
    lut = GetColorTransferFunction(field)
    key = "pressure" if field == "pressure" else \
          "u" if field == "u_comp" else \
          "v" if field == "v_comp" else "vel_mag"
    if key == "pressure":
        d.RescaleTransferFunctionToDataRange(True)
        Render()
    else:
        apply_cmap(lut, key)
    d.SetScalarBarVisibility(v, True)
    bar = GetScalarBar(lut, v)
    bar.Title = title; bar.ComponentTitle = ""
    bar.TitleFontSize = 13; bar.LabelFontSize = 11
    setup_camera(v); Render()
    save(v, OUT_DIR + "/" + fname)
    Delete(v)

# ── Glyphen (Pfeile) ─────────────────────────────────────────────────────────
print("  Rendere Glyphen (Pfeile)…")
v_g = new_view()
d_bg = Show(calc_final, v_g)
d_bg.Representation = "Surface"
ColorBy(d_bg, ("CELLS", "vel_mag"))
lut_g = GetColorTransferFunction("vel_mag")
apply_cmap(lut_g, "vel_mag")
d_bg.Opacity = 0.4

glyph = Glyph(Input=calc_final, GlyphType="Arrow")
glyph.OrientationArray  = ["CELLS", "velocity"]
glyph.ScaleArray        = ["CELLS", "vel_mag"]
glyph.ScaleFactor       = 0.04
glyph.MaximumNumberOfSamplePoints = 400
glyph.GlyphMode         = "Every Nth Point"
glyph.Stride            = 6

d_gl = Show(glyph, v_g)
d_gl.Representation = "Surface"
ColorBy(d_gl, ("POINTS", "vel_mag"))
lut_gl = GetColorTransferFunction("vel_mag")
apply_cmap(lut_gl, "vel_mag")
d_gl.SetScalarBarVisibility(v_g, True)
bar_gl = GetScalarBar(lut_gl, v_g)
bar_gl.Title = "Velocity |u| [m/s]"; bar_gl.ComponentTitle = ""

setup_camera(v_g); Render()
save(v_g, OUT_DIR + "/final_glyphs.png")
Delete(v_g)

# ── Stromlinien ───────────────────────────────────────────────────────────────
print("  Rendere Stromlinien…")
v_s = new_view(bg=(0.05, 0.05, 0.1))
d_sbg = Show(calc_final, v_s)
d_sbg.Representation = "Surface"
ColorBy(d_sbg, ("CELLS", "vel_mag"))
lut_s = GetColorTransferFunction("vel_mag")
apply_cmap(lut_s, "vel_mag")
d_sbg.Opacity = 0.5

stream = StreamTracer(Input=rd_final, SeedType="Line")
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
bar_st = GetScalarBar(lut_st, v_s)
bar_st.Title = "Velocity |u| [m/s]"; bar_st.ComponentTitle = ""

setup_camera(v_s); Render()
save(v_s, OUT_DIR + "/final_streamlines.png")
Delete(v_s)

# Speicher freigeben
Delete(calc_final); Delete(rd_final)

# ══════════════════════════════════════════════════════════════════════════════
# TEIL 2 – Animationsframes für alle 4 Größen (→ später Video)
# ══════════════════════════════════════════════════════════════════════════════
print("\n━━━ TEIL 2: Animationsframes (alle Zeitschritte) ━━━")

rd_anim, calc_anim = load_step(VTK_FILES)
times = rd_anim.TimestepValues
if not times:
    times = list(range(N))
print(f"  {len(times)} Frames × 4 Größen…")

anim = GetAnimationScene()
anim.UpdateAnimationUsingDataTimeSteps()

FIELDS_ANIM = [
    ("u_comp",   "u",        "u-Velocity [m/s]",           "u"),
    ("v_comp",   "v",        "v-Velocity [m/s]",           "v"),
    ("pressure", "pressure", "Pressure [Pa]",               "p"),
    ("vel_mag",  "vel_mag",  "Velocity Magnitude |u| [m/s]","vel"),
]

views_anim = []
disps_anim = []
for field, cmap_key, title, prefix in FIELDS_ANIM:
    va = new_view(size=(800, 700))
    da = Show(calc_anim, va)
    da.Representation = "Surface"
    ColorBy(da, ("CELLS", field))
    lut_a = GetColorTransferFunction(field)
    if cmap_key == "pressure":
        da.RescaleTransferFunctionToDataRange(True)
    else:
        apply_cmap(lut_a, cmap_key)
    da.SetScalarBarVisibility(va, True)
    bar_a = GetScalarBar(lut_a, va)
    bar_a.Title = title; bar_a.ComponentTitle = ""
    bar_a.TitleFontSize = 11; bar_a.LabelFontSize = 9
    setup_camera(va)
    views_anim.append((va, prefix))
    disps_anim.append((da, lut_a, cmap_key))

for idx, t in enumerate(times):
    anim.AnimationTime = t
    for (va, prefix), (da, lut_a, cmap_key) in zip(views_anim, disps_anim):
        if cmap_key == "pressure":
            da.RescaleTransferFunctionToDataRange(True)
        Render()
        SaveScreenshot(f"{FRM_DIR}/{prefix}_{idx:04d}.png", va,
                       ImageResolution=[800, 700])
    if idx % 10 == 0:
        print(f"    Frame {idx:03d}/{len(times)-1}")

for va, _ in views_anim:
    Delete(va)
Delete(calc_anim); Delete(rd_anim)

# ══════════════════════════════════════════════════════════════════════════════
# ABSCHLUSS
# ══════════════════════════════════════════════════════════════════════════════
print(f"""
━━━ FERTIG ━━━
Endzustand-Bilder:
  {OUT_DIR}/final_u.png
  {OUT_DIR}/final_v.png
  {OUT_DIR}/final_pressure.png
  {OUT_DIR}/final_velocity.png
  {OUT_DIR}/final_glyphs.png
  {OUT_DIR}/final_streamlines.png

Animationsframes: {FRM_DIR}/[u|v|p|vel]_NNNN.png

Jetzt Videos erzeugen mit:
  ./make_videos.sh
""")
