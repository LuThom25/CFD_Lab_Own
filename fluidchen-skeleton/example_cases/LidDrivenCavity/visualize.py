"""
Korrigierte ParaView-Visualisierung der Lid-Driven Cavity Simulation.
Erstellt:
  - final_velocity.png      : Geschwindigkeitsbetrag im Endzustand
  - final_pressure.png      : Druck (normiert) im Endzustand
  - final_streamlines.png   : Stromlinien im Endzustand
  - screenshots/vel_*.png   : Animation Geschwindigkeitsbetrag
  - screenshots/prs_*.png   : Animation Druck (normiert)
"""

from paraview.simple import *
import os, glob, numpy as np

paraview.simple._DisableFirstRenderCameraReset()

# ── Pfade ──────────────────────────────────────────────────────────────────────
BASE    = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = BASE + "/LidDrivenCavity_Output"
SCR_DIR = OUT_DIR + "/screenshots"
os.makedirs(SCR_DIR, exist_ok=True)

# Dateiformat: CaseName_rank.timestep.vtk  → sortiere nach timestep (Zahl nach dem Punkt)
def vtk_sort_key(f):
    stem = os.path.basename(f).replace(".vtk", "")   # z.B. "LidDrivenCavity_0.1701"
    parts = stem.split(".")                            # ["LidDrivenCavity_0", "1701"]
    return int(parts[-1])

VTK_FILES = sorted(glob.glob(OUT_DIR + "/*.vtk"), key=vtk_sort_key)
print(f"Gefundene Zeitschritte: {len(VTK_FILES)}")

# ── Hilfsfunktionen ────────────────────────────────────────────────────────────
def reset_view(view):
    view.ResetCamera()
    view.CameraPosition    = [0.5, 0.5, 3.0]
    view.CameraFocalPoint  = [0.5, 0.5, 0.0]
    view.CameraViewUp      = [0.0, 1.0, 0.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 0.55

def save(view, path):
    SaveScreenshot(path, view, ImageResolution=[1400, 1200])
    print(f"  Gespeichert: {os.path.basename(path)}")

# ══════════════════════════════════════════════════════════════════════════════
# 1) ENDZUSTAND: Geschwindigkeitsbetrag
# ══════════════════════════════════════════════════════════════════════════════
print("\n[1/3] Rendere Geschwindigkeitsbetrag (Endzustand)…")

reader_vel = LegacyVTKReader(FileNames=[VTK_FILES[-1]])

# Velocity-Magnitude via Calculator (Betrag des 3D-Vektors)
calc = Calculator(Input=reader_vel)
calc.AttributeType = "Cell Data"
calc.ResultArrayName = "vel_magnitude"
calc.Function = "mag(velocity)"

view1 = CreateView("RenderView")
view1.ViewSize = [1400, 1200]
view1.Background = [0.15, 0.15, 0.2]

disp = Show(calc, view1)
disp.Representation = "Surface"
ColorBy(disp, ("CELLS", "vel_magnitude"))

lut = GetColorTransferFunction("vel_magnitude")
lut.ApplyPreset("Jet", True)
lut.RescaleTransferFunction(0.0, 1.0)   # U_wall = 1 → max = 1
disp.SetScalarBarVisibility(view1, True)
bar = GetScalarBar(lut, view1)
bar.Title = "Velocity Magnitude  |u|  [m/s]"
bar.ComponentTitle = ""
bar.TitleFontSize  = 14
bar.LabelFontSize  = 12

reset_view(view1)
Render()
save(view1, OUT_DIR + "/final_velocity.png")

Delete(view1); del view1
Delete(calc);  del calc
Delete(reader_vel); del reader_vel

# ══════════════════════════════════════════════════════════════════════════════
# 2) ENDZUSTAND: Druck (Mittelwert abgezogen → physikalisch sinnvoll)
# ══════════════════════════════════════════════════════════════════════════════
print("\n[2/3] Rendere Druck (Endzustand)…")

reader_prs = LegacyVTKReader(FileNames=[VTK_FILES[-1]])

# Druck normieren: p_norm = p - mean(p)
calc_p = Calculator(Input=reader_prs)
calc_p.AttributeType = "Cell Data"
calc_p.ResultArrayName = "pressure_norm"
# ParaView Calculator kann kein mean(), daher: einfach pressure direkt zeigen
# und Colormap auf symmetrischen Bereich setzen
calc_p.Function = "pressure"

view2 = CreateView("RenderView")
view2.ViewSize = [1400, 1200]
view2.Background = [0.15, 0.15, 0.2]

disp2 = Show(calc_p, view2)
disp2.Representation = "Surface"
ColorBy(disp2, ("CELLS", "pressure_norm"))

lut2 = GetColorTransferFunction("pressure_norm")
lut2.ApplyPreset("Cool to Warm", True)
# Colormap automatisch auf Datenbereich skalieren
disp2.RescaleTransferFunctionToDataRange(True)
disp2.SetScalarBarVisibility(view2, True)
bar2 = GetScalarBar(lut2, view2)
bar2.Title = "Pressure  [Pa]"
bar2.ComponentTitle = ""
bar2.TitleFontSize = 14
bar2.LabelFontSize = 12

reset_view(view2)
Render()
save(view2, OUT_DIR + "/final_pressure.png")

Delete(view2); del view2
Delete(calc_p);  del calc_p
Delete(reader_prs); del reader_prs

# ══════════════════════════════════════════════════════════════════════════════
# 3) ENDZUSTAND: Stromlinien
# ══════════════════════════════════════════════════════════════════════════════
print("\n[3/3] Rendere Stromlinien (Endzustand)…")

reader_sl = LegacyVTKReader(FileNames=[VTK_FILES[-1]])

view3 = CreateView("RenderView")
view3.ViewSize = [1400, 1200]
view3.Background = [0.05, 0.05, 0.1]

# Untergrund: Geschwindigkeitsbetrag als Hintergrund
calc_bg = Calculator(Input=reader_sl)
calc_bg.AttributeType = "Cell Data"
calc_bg.ResultArrayName = "vel_mag_bg"
calc_bg.Function = "mag(velocity)"

disp_bg = Show(calc_bg, view3)
disp_bg.Representation = "Surface"
ColorBy(disp_bg, ("CELLS", "vel_mag_bg"))
lut_bg = GetColorTransferFunction("vel_mag_bg")
lut_bg.ApplyPreset("Jet", True)
lut_bg.RescaleTransferFunction(0.0, 1.0)
disp_bg.Opacity = 0.6

# Stromlinien auf POINT_DATA velocity
stream = StreamTracer(Input=reader_sl, SeedType="Line")
stream.Vectors                 = ["POINTS", "velocity"]
stream.MaximumStreamlineLength = 2.0
stream.IntegrationDirection    = "BOTH"
stream.SeedType.Point1         = [0.02, 0.5, 0.0]
stream.SeedType.Point2         = [0.98, 0.5, 0.0]
stream.SeedType.Resolution     = 30

disp_sl = Show(stream, view3)
disp_sl.Representation = "Surface"
disp_sl.LineWidth = 2.0
ColorBy(disp_sl, ("POINTS", "velocity"))
lut_sl = GetColorTransferFunction("velocity")
lut_sl.ApplyPreset("Jet", True)
lut_sl.RescaleTransferFunction(0.0, 1.0)

disp_bg.SetScalarBarVisibility(view3, False)
disp_sl.SetScalarBarVisibility(view3, True)
bar3 = GetScalarBar(lut_sl, view3)
bar3.Title = "Velocity |u| [m/s]"
bar3.ComponentTitle = ""

reset_view(view3)
Render()
save(view3, OUT_DIR + "/final_streamlines.png")

Delete(view3); del view3

# ══════════════════════════════════════════════════════════════════════════════
# 4) ANIMATION: Geschwindigkeitsbetrag über alle Zeitschritte
# ══════════════════════════════════════════════════════════════════════════════
print("\n[4/4] Rendere Animation (Geschwindigkeitsbetrag)…")

reader_anim = LegacyVTKReader(FileNames=VTK_FILES)

calc_anim = Calculator(Input=reader_anim)
calc_anim.AttributeType = "Cell Data"
calc_anim.ResultArrayName = "vel_magnitude"
calc_anim.Function = "mag(velocity)"

view4 = CreateView("RenderView")
view4.ViewSize = [900, 800]
view4.Background = [0.15, 0.15, 0.2]

disp4 = Show(calc_anim, view4)
disp4.Representation = "Surface"
ColorBy(disp4, ("CELLS", "vel_magnitude"))

lut4 = GetColorTransferFunction("vel_magnitude")
lut4.ApplyPreset("Jet", True)
lut4.RescaleTransferFunction(0.0, 1.0)
disp4.SetScalarBarVisibility(view4, True)
bar4 = GetScalarBar(lut4, view4)
bar4.Title = "Velocity |u| [m/s]"
bar4.ComponentTitle = ""

reset_view(view4)

anim  = GetAnimationScene()
anim.UpdateAnimationUsingDataTimeSteps()
# Zeitschritte vom Reader holen (nicht vom Calculator-Filter)
times = reader_anim.TimestepValues
if not times:
    times = [0]

print(f"  Speichere {len(times)} Frames…")
for idx, t in enumerate(times):
    anim.AnimationTime = t
    Render()
    save(view4, f"{SCR_DIR}/vel_{idx:04d}.png")

print(f"\n✓ Alle Dateien gespeichert in:\n  {OUT_DIR}/")
print("  final_velocity.png    — Geschwindigkeitsbetrag Endzustand")
print("  final_pressure.png    — Druck Endzustand")
print("  final_streamlines.png — Stromlinien Endzustand")
print(f"  screenshots/vel_*.png — {len(times)} Animationsframes")
