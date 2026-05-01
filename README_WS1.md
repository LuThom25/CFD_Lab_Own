# Worksheet 1 — Lid-Driven Cavity

## Requirements

- CMake + GCC (available in WSL/Linux)
- ParaView 6.0.1 (with `pvpython`) — Windows or Linux
- Python 3 with `matplotlib` and `pillow`
- `ffmpeg` for video generation

Install Python dependencies (WSL):
```bash
sudo apt install python3-matplotlib python3-pil ffmpeg
# OR inside a venv:
pip install matplotlib pillow
```

---

## 1. Compile

```bash
cd ~/CFD_Lab/cfd-lab-group-c/build
cmake ..
make -j$(nproc)
```

This produces the binary at `build/fluidchen`.

> Recompile with `make -j$(nproc)` every time you change C++ source files.

---

## 2. Run everything

From anywhere in the terminal:

```bash
cd ~/CFD_Lab/cfd-lab-group-c/example_cases/LidDrivenCavity && \
../../build/fluidchen LidDrivenCavity.dat && \
"/mnt/c/Program Files/ParaView 6.0.1/bin/pvpython.exe" visualize.py && \
bash make_videos.sh && \
python3 run_studies.py 5 6 7 8
```

### What each step does

| Command | What it does |
|---|---|
| `fluidchen LidDrivenCavity.dat` | Runs the CFD simulation, writes VTK files to `LidDrivenCavity_Output/` |
| `pvpython visualize.py` | Generates static images and animation frames from the VTK output |
| `bash make_videos.sh` | Assembles frames into `.mp4` videos using ffmpeg |
| `python3 run_studies.py 5 6 7 8` | Runs parameter studies for Tasks 5–8 and saves plots |

### To run only specific tasks

```bash
python3 run_studies.py 4        # Task 4 only
python3 run_studies.py 5 6      # Tasks 5 and 6
python3 run_studies.py 5 6 7 8  # Tasks 5 to 8
```

---

## 3. ParaView path

The `pvpython` path is hardcoded in the run command above. If your ParaView is installed elsewhere, find it with:

```bash
find /mnt/c -maxdepth 4 -iname "pvpython.exe" 2>/dev/null
```

Then replace `/mnt/c/Program Files/ParaView 6.0.1/bin/pvpython.exe` in the command accordingly.

On native Linux, `pvpython` is usually just:
```bash
pvpython visualize.py
```

---

## 4. Output

All results are saved to:
```
example_cases/LidDrivenCavity/LidDrivenCavity_Output/
├── task4/
│   ├── final_u.png
│   ├── final_v.png
│   ├── final_pressure.png
│   ├── final_velocity.png
│   ├── final_glyphs.png
│   ├── final_streamlines.png
│   ├── frames/         # animation frames
│   └── videos/         # .mp4 files
└── study_plots/        # plots for tasks 5–8
```

---

## 5. Simulation parameters (Task 4)

Defined in `LidDrivenCavity.dat`:

| Parameter | Value |
|---|---|
| imax, jmax | 50 × 50 |
| xlength, ylength | 1.0 |
| dt | 0.05 (adaptive, tau=0.5) |
| t_end | 50.0 |
| dt_value | 0.5 |
| eps | 0.001 |
| omg | 1.7 |
| gamma | 0.5 |
| itermax | 100 |
| nu | 0.01 (Re ≈ 100) |

---

## 6. Notes

- Task 5 runs 9 full simulations (one per omega value) — expect ~30 min total.
- Tasks 6–8 are faster since they use shorter `t_end` or fewer grid points.
- If `make_videos.sh` fails, check that `ffmpeg` is installed (`sudo apt install ffmpeg`).
- If `matplotlib` is missing: `pip install matplotlib pillow` inside your venv, or `sudo apt install python3-matplotlib python3-pil` system-wide.
