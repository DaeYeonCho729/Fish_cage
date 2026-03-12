import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# ── Parameters ──────────────────────────────────────────────
Sn   = 0.19
dw   = 0.00242
Urel = 0.5
eps  = 1e-10
nu   = 1.05e-6

angles_deg = np.linspace(0, 90, 500)
angles_rad = np.deg2rad(angles_deg)
phi = np.pi / 2 - angles_rad


# ── Drag Models ─────────────────────────────────────────────
def bessonneau_cd(phi):
    Cn, Ct = 1.2, 0.1
    return Cn * np.cos(phi) + Ct * np.sin(phi)

def lee_cd_raw(deg):
    cd = np.zeros_like(deg, dtype=float)
    cd = np.where((deg >= 0)  & (deg <= 10), 0.00185  * deg + 0.033352, cd)
    cd = np.where((deg > 10)  & (deg <= 20), 0.005986 * deg - 0.01081,  cd)
    cd = np.where((deg > 20)  & (deg <= 40), 0.015606 * deg - 0.2101,   cd)
    cd = np.where((deg > 40)  & (deg <= 60), 0.0239   * deg - 0.54333,  cd)
    cd = np.where((deg > 60)  & (deg <= 70), 0.011907 * deg + 0.183297, cd)
    cd = np.where((deg > 70)  & (deg <= 90), 0.008774 * deg + 0.402622, cd)
    return cd

def lee_cd_smooth(angles):
    bp = np.array([0, 10, 20, 40, 60, 70, 90], dtype=float)
    spline = make_interp_spline(bp, lee_cd_raw(bp), k=3)
    return spline(angles)

def loland_cd(phi):
    return 0.04 + (-0.04 + 0.33*Sn + 6.54*Sn**2 - 4.88*Sn**3) * np.cos(phi)

def kristiansen_cd(phi):
    Re      = (dw * Urel) / (nu * (1.0 - Sn) + eps)
    x       = np.log10(Re)
    cd_circ = (-78.46675 + 254.73873*x - 327.88640*x**2 + 223.64577*x**3
               - 87.92234*x**4 + 20.00769*x**5 - 2.44894*x**6 + 0.12479*x**7)
    cd      = cd_circ * Sn * (2.0 - Sn) / (((1.0 - Sn)**2) * 2)
    return cd * np.cos(phi)


# ── Lift Models ─────────────────────────────────────────────
def bessonneau_cl(phi):
    Cn, Ct = 1.2, 0.1
    return Cn * np.sin(phi) - Ct * np.cos(phi)

def lee_cl_raw(deg):
    cl = np.zeros_like(deg, dtype=float)
    cl = np.where((deg >= 0)  & (deg <= 5),  0.001831 * deg + 0.020675,  cl)
    cl = np.where((deg > 5)   & (deg <= 10), 0.004817 * deg + 0.004268,  cl)
    cl = np.where((deg > 10)  & (deg <= 20), 0.008873 * deg - 0.03891,   cl)
    cl = np.where((deg > 20)  & (deg <= 30), 0.008873 * deg - 0.11956,   cl)
    cl = np.where((deg > 30)  & (deg <= 43), 0.013325 * deg - 0.13140,   cl)
    cl = np.where((deg > 43)  & (deg <= 50), 0.008639 * deg + 0.070035,  cl)
    cl = np.where((deg > 50)  & (deg <= 55), 0.002428 * deg + 0.383586,  cl)
    cl = np.where((deg > 55)  & (deg <= 60), -0.00503 * deg + 0.797889,  cl)
    cl = np.where((deg > 60)  & (deg <= 65), -0.01468 * deg + 1.38167,   cl)
    cl = np.where((deg > 65)  & (deg <= 72), -0.02828 * deg + 2.274163,  cl)
    cl = np.where((deg > 72)  & (deg <= 85), -0.01127 * deg + 1.046885,  cl)
    cl = np.where(deg > 85,                  -0.00703 * deg + 0.689754,  cl)
    return cl

def lee_cl_smooth(angles):
    bp = np.array([0, 5, 10, 20, 30, 43, 50, 55, 60, 65, 72, 85, 90], dtype=float)
    spline = make_interp_spline(bp, lee_cl_raw(bp), k=3)
    return spline(angles)

def loland_cl(phi):
    return (-0.05*Sn + 2.3*Sn**2 - 1.76*Sn**3) * np.sin(2 * phi)

def kristiansen_cl(phi):
    Re      = (dw * Urel) / (nu * (1.0 - Sn) + eps)
    x       = np.log10(Re)
    cd_circ = (-78.46675 + 254.73873*x - 327.88640*x**2 + 223.64577*x**3
               - 87.92234*x**4 + 20.00769*x**5 - 2.44894*x**6 + 0.12479*x**7)
    cd      = cd_circ * Sn * (2.0 - Sn) / (((1.0 - Sn)**2) * 2)
    cn_45   = 0.5 * cd
    y_45    = 0.25 * np.pi
    ct_45   = y_45 * (4.0 * cn_45 / (8.0 + cn_45))
    cl      = (cn_45 - ct_45) / np.sqrt(2)
    return cl * np.sin(phi) * 2


# ── Compute ──────────────────────────────────────────────────
cd_data = [bessonneau_cd(phi), lee_cd_smooth(angles_deg),
           loland_cd(phi),     kristiansen_cd(phi)]
cl_data = [bessonneau_cl(phi), lee_cl_smooth(angles_deg),
           loland_cl(phi),     kristiansen_cl(phi)]

styles = [
    {'marker': 's', 'label': '1D model (Bessonneau 1998)'},
    {'marker': 'D', 'label': '1D model (Lee et al. 2005)'},
    {'marker': 'o', 'label': '2D model (Loland 1991)'},
    {'marker': '^', 'label': '2D model (Kristiansen 2012)'},
]


# ── Plot ─────────────────────────────────────────────────────
plt.rcParams['font.family']      = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

for ax, data, ylabel, title in [
    (ax1, cd_data, 'Drag Coefficient',  'Drag Coefficient'),
    (ax2, cl_data, 'Lift Coefficient',  'Lift Coefficient'),
]:
    for cd, style in zip(data, styles):
        ax.plot(angles_deg, cd,
                color='black', linestyle='-',
                marker=style['marker'], markevery=50,
                markersize=5, linewidth=1.2,
                label=style['label'])

    ax.set_xlabel('Attack Angle (°)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlim(0, 90)
    ax.tick_params(direction='in', labelsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)

ax1.legend(fontsize=9.5, frameon=True, edgecolor='black', loc='upper left')
ax2.legend(fontsize=9.5, frameon=True, edgecolor='black', loc='upper right')

plt.tight_layout()
plt.savefig('cd_cl_comparison.png', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.show()