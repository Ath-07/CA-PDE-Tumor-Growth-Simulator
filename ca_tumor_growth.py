import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import imageio
import os

# cell states 
EMPTY = 0
TUMOR = 1
QUIESCENT = 2
NECROTIC = 3

# parameters 
N = 200            # grid size
T = 3000           # timesteps
capture_every = 2  # capture frame every k steps for GIF size control

# nutrient/PDE params (if use_pde==False, a static gradient is used)
use_pde = True
D = 1.2            # diffusion coefficient (increase so nutrients reach center faster)
dt = 0.25          # PDE timestep (slightly larger for faster diffusion)
lam = 0.04         # consumption rate (reduce so tumor doesn't starve center immediately)

# division thresholds / probabilities
div_threshold = 0.25        # lower threshold so division is possible at moderate nutrient
quiescent_threshold = 0.15  # lower quiescence threshold so cells stay proliferative longer
necrosis_threshold = 0.06
p_div = 0.95                # increase max division probability (more aggressive)
base_div_chance = 0.08      # larger baseline chance for stochastic growth


# neighborhood offsets (Moore model)
neighbors = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

# ---- helper functions ----
def nutrient_field_static(N):
    x, y = np.meshgrid(np.linspace(-1,1,N), np.linspace(-1,1,N))
    dist = np.sqrt(x*x + y*y)

    # nutrient high at boundary, low at center
    nut = dist
    nut = np.clip(nut, 0.0, 1.0)

    return nut


def laplacian(A):
    up    = np.roll(A, -1, axis=0)
    down  = np.roll(A, 1, axis=0)
    left  = np.roll(A, -1, axis=1)
    right = np.roll(A, 1, axis=1)
    return (up + down + left + right - 4*A)

def update_nutrient_pde(nutrient, grid, D, lam, dt):
    """Explicit Euler update of reaction-diffusion PDE:
       ∂N/∂t = D ∇^2 N - lam * (tumor|quiescent) * N + source
       Here source is only boundary (kept high).
    """
    lap = laplacian(nutrient)
    consumer = ((grid == TUMOR) | (grid == QUIESCENT)).astype(float)
    consumption = lam * consumer * nutrient
    # simple source: boundary values act as blood supply -> enforce after step
    N_new = nutrient + dt * (D * lap - consumption)
    # clamp
    N_new = np.clip(N_new, 0.0, 1.0)
    # enforce boundary supply
    N_new[0, :] = np.maximum(N_new[0, :], 1.0)
    N_new[-1, :] = np.maximum(N_new[-1, :], 1.0)
    N_new[:, 0] = np.maximum(N_new[:, 0], 1.0)
    N_new[:, -1] = np.maximum(N_new[:, -1], 1.0)
    return N_new

def attempt_divisions(grid, nutrient):
    """For each proliferating tumor cell, attempt to divide into a random empty neighbor."""
    births = np.zeros_like(grid, dtype=np.uint8)
    rng = np.random.default_rng()
    tumor_positions = np.argwhere(grid == TUMOR)
    rng.shuffle(tumor_positions)
    for (i, j) in tumor_positions:
        empties = []
        for di, dj in neighbors:
            x, y = i + di, j + dj
            if 0 <= x < grid.shape[0] and 0 <= y < grid.shape[1] and grid[x, y] == EMPTY:
                empties.append((x, y))
        if not empties:
            continue
        # scaled chance based on nutrient
        # allow nutrient below div_threshold to still contribute (so moderate nutrient isn't ignored)
        scaled = (nutrient[i,j] - div_threshold) / (1.0 - div_threshold + 1e-9)
        # soften negative scaled values so they still add a little probability (but not too much)
        scaled = np.clip(scaled, -0.5, 1.0)    # allow negative but not too negative
        prob = p_div * np.clip(base_div_chance + max(0.0, scaled), 0.0, 1.0)

        if rng.random() < prob:
            x, y = empties[rng.integers(len(empties))]
            births[x, y] = 1
    return births

# ---- main simulate function ----
def simulate_and_animate(N=N, T=T, save_path="tumor_ca_animation.gif"):
    # initialize CA grid
    grid = np.zeros((N, N), dtype=np.uint8)
    center = (N//2, N//2)
    grid[center] = TUMOR

    # initialize nutrient
    if use_pde:
        # start interior at a higher baseline so diffusion can feed the seed faster
        nutrient = nutrient_field_static(N) * 0.6
        # reinforce boundaries to be high (blood supply)
        nutrient[0, :] = 1.0
        nutrient[-1, :] = 1.0
        nutrient[:, 0] = 1.0
        nutrient[:, -1] = 1.0
    else:
        nutrient = nutrient_field_static(N)

    # very important: boost seed nutrient so initial tumor can proliferate
    nutrient[center] = max(nutrient[center], 0.6)

    # prepare figure (but using Agg backend - will capture frames via buffer_rgba())
    fig, ax = plt.subplots(figsize=(5,5))
    cmap = colors.ListedColormap(['white', 'tab:red', 'tab:orange', 'black'])
    norm = colors.BoundaryNorm([0,1,2,3,4], cmap.N)

    frames = []
    for t in range(T):
        if t % 50 == 0:
            print(f"t={t} center_nutrient={nutrient[center]:.3f} tumor_cells={np.sum(grid==TUMOR)}")

        # 1) nutrient update (PDE) if enabled
        if use_pde:
            nutrient = update_nutrient_pde(nutrient, grid, D, lam, dt)

        # 2) state transitions based on nutrient
        # Tumor -> Quiescent if nutrient low
        grid[(grid == TUMOR) & (nutrient < quiescent_threshold)] = QUIESCENT
        # Quiescent -> Necrotic if very low
        grid[(grid == QUIESCENT) & (nutrient < necrosis_threshold)] = NECROTIC

        # 3) proliferation attempts
        births = attempt_divisions(grid, nutrient)
        grid[births == 1] = TUMOR

        # 4) small random invasion (mechanical pushing)
        tumor_mask = (grid == TUMOR)
        takeover_candidates = np.zeros_like(grid, dtype=bool)
        for di, dj in neighbors:
            shifted = np.roll(np.roll(tumor_mask, di, axis=0), dj, axis=1)
            takeover_candidates |= (shifted & (grid == EMPTY))
        rng = np.random.default_rng()
        takeover = (rng.random(grid.shape) < 0.0005) & takeover_candidates
        grid[takeover] = TUMOR

        # 5) capture frame
        if (t % capture_every) == 0 or t == T-1:
            ax.clear()
            ax.axis('off')
            ax.set_title(f"timestep {t}")
            ax.imshow(grid, cmap=cmap, norm=norm, interpolation='nearest')

            # draw and capture via buffer_rgba (safe for headless)
            fig.canvas.draw()
            buf = fig.canvas.buffer_rgba()
            frame = np.asarray(buf)        # RGBA uint8
            frame_rgb = frame[..., :3].copy()
            frames.append(frame_rgb)

        # optional: small feedback printed occasionally
        if t % (T//5) == 0:
            print(f"Progress {t}/{T}  center_nutrient={nutrient[center]:.3f}  tumor_cells={np.sum(grid==TUMOR)}")

    plt.close(fig)

    # save GIF
    out_dir = os.getcwd()
    gif_full = os.path.join(out_dir, save_path)
    imageio.mimsave(gif_full, frames, fps=12)
    print("Saved animation to:", gif_full)
    return gif_full, grid, nutrient

# ---- run if main ----
if __name__ == "__main__":
    out_gif, final_grid, final_nutrient = simulate_and_animate()
    print(nutrient_field_static(5))
