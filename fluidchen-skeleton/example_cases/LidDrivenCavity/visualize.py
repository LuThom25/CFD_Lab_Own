"""
Automatische ParaView-Visualisierung der Lid-Driven Cavity Simulation.
Ladet alle VTK-Zeitschritte, setzt Velocity-Farben und speichert Screenshots + Animation.
"""

from paraview.simple import *
import os

# Pfad zum Output-Ordner
output_dir = os.path.dirname(os.path.abspath(__file__)) + "/LidDrivenCavity_Output"
screenshots_dir = output_dir + "/screenshots"
os.makedirs(screenshots_dir, exist_ok=True)

print("Lade VTK-Dateien...")

# Alle VTK-Dateien als Zeitreihe laden
import glob
vtk_files = sorted(glob.glob(output_dir + "/*.vtk"))
print(f"  {len(vtk_files)} Zeitschritte gefunden.")

reader = LegacyVTKReader(FileNames=vtk_files)
RenderAllViews()

# View einrichten
view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1200, 900]
view.Background = [0.1, 0.1, 0.1]

# Display-Objekt holen
display = Show(reader, view)
display.Representation = "Surface"

# --- Velocity-Farbe einstellen ---
print("Stelle Velocity-Farbe ein...")
ColorBy(display, ("CELLS", "velocity"))

# Farb-Skala anpassen (Regenbogen-Colormap)
velocityLUT = GetColorTransferFunction("velocity")
velocityLUT.ApplyPreset("Cool to Warm", True)

# Farbbalken anzeigen
display.SetScalarBarVisibility(view, True)
colorBar = GetScalarBar(velocityLUT, view)
colorBar.Title = "Velocity Magnitude"
colorBar.ComponentTitle = ""

# Kamera auf 2D-Ansicht setzen (von oben)
view.ResetCamera()
view.CameraPosition    = [0.5, 0.5, 3.0]
view.CameraFocalPoint  = [0.5, 0.5, 0.0]
view.CameraViewUp      = [0.0, 1.0, 0.0]
view.CameraParallelProjection = 1

Render()

# --- Screenshots für jeden Zeitschritt speichern ---
print("Speichere Screenshots...")
animScene = GetAnimationScene()
animScene.UpdateAnimationUsingDataTimeSteps()

timekeeper = animScene.TimeKeeper
times = reader.TimestepValues
if not times:
    times = [0]

for i, t in enumerate(times):
    timekeeper.Time = t
    Render()
    screenshot_path = f"{screenshots_dir}/frame_{i:04d}.png"
    SaveScreenshot(screenshot_path, view, ImageResolution=[1200, 900])

print(f"  {len(times)} Screenshots gespeichert in:\n  {screenshots_dir}/")

# --- Letzten Zeitschritt als finales Bild speichern ---
final_path = output_dir + "/final_state_velocity.png"
timekeeper.Time = times[-1]
Render()
SaveScreenshot(final_path, view, ImageResolution=[1200, 900])
print(f"\nFinales Bild (t={times[-1]:.2f}) gespeichert:\n  {final_path}")

print("\nFertig!")
