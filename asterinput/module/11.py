import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# ── Parameters ──────────────────────────────────────────────
Sn   = 0.19
dw   = 0.00141
Urel = 0.5
eps  = 1e-10
nu   = 1.05e-6

angles_deg = np.linspace(0, 90, 500)
angles_rad = np.deg2rad(angles_deg)
phi = np.pi / 2 - angles_rad


# ── Models ──────────────────────────────────────────────────
def bessonneau(phi):
    Cn, Ct = 1.2, 0.1
    return Cn * np.cos(phi) + Ct * np.sin(phi)

def lee_raw(deg):
    cd = np.zeros_like(deg, dtype=float)
    cd = np.where((deg >= 0)  & (deg <= 10), 0.00185  * deg + 0.033352, cd)
    cd = np.where((deg > 10)  & (deg <= 20), 0.005986 * deg - 0.01081,  cd)
    cd = np.where((deg > 20)  & (deg <= 40), 0.015606 * deg - 0.2101,   cd)
    cd = np.where((deg > 40)  & (deg <= 60), 0.0239   * deg - 0.54333,  cd)
    cd = np.where((deg > 60)  & (deg <= 70), 0.011907 * deg + 0.183297, cd)
    cd = np.where((deg > 70)  & (deg <= 90), 0.008774 * deg + 0.402622, cd)
    return cd

def lee_smooth(angles):
    breakpoints = np.array([0, 10, 20, 40, 60, 70, 90], dtype=float)
    cd_breaks   = lee_raw(breakpoints)
    spline      = make_interp_spline(breakpoints, cd_breaks, k=3)
    return spline(angles)

def loland(phi):
    return 0.04 + (-0.04 + 0.33*Sn + 6.54*Sn**2 - 4.88*Sn**3) * np.cos(phi)

def kristiansen(phi):
    Re      = (dw * Urel) / (nu * (1.0 - Sn) + eps)
    x       = np.log10(Re)
    cd_circ = (-78.46675 + 254.73873*x - 327.88640*x**2 + 223.64577*x**3
               - 87.92234*x**4 + 20.00769*x**5 - 2.44894*x**6 + 0.12479*x**7)
    cd      = cd_circ * Sn * (2.0 - Sn) / (((1.0 - Sn)**2) * 2)
    cn_45   = 0.5 * cd
    y_45    = 0.25 * np.pi
    ct_45   = y_45 * (4.0 * cn_45 / (8.0 + cn_45))
    return cd * np.cos(phi)

cd_bessonneau  = bessonneau(phi)
cd_lee         = lee_smooth(angles_deg)
cd_loland      = loland(phi)
cd_kristiansen = kristiansen(phi)
# ── Difference & Mean ───────────────────────────────────────
cd_loland_mean = np.mean(cd_loland)
cd_kristiansen_mean = np.mean(cd_kristiansen)

percent_higher = (cd_loland_mean - cd_kristiansen_mean) / cd_kristiansen_mean * 100

print(f"Loland mean Cd        : {cd_loland_mean:.6f}")
print(f"Kristiansen mean Cd   : {cd_kristiansen_mean:.6f}")
print(f"Mean difference       : {cd_loland_mean - cd_kristiansen_mean:.6f}")
print(f"Loland is {percent_higher:.2f}% higher than Kristiansen")

# ── Plot ────────────────────────────────────────────────────
plt.rcParams['font.family']      = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots(figsize=(7, 5))

styles = [
    {'marker': 's', 'label': '1D model (Bessonneau 1998)'},
    {'marker': 'D', 'label': '1D model (Lee et al. 2005)'},
    {'marker': 'o', 'label': '2D model (Loland 1991)'},
    {'marker': '^', 'label': '2D model (Kristiansen 2012)'},
]
cds = [cd_bessonneau, cd_lee, cd_loland, cd_kristiansen]

for cd, style in zip(cds, styles):
    ax.plot(angles_deg, cd,
            color='black',
            linestyle='-',
            marker=style['marker'],
            markevery=50,
            markersize=5,
            linewidth=1.2,
            label=style['label'])

ax.set_xlabel('Attack Angle (°)', fontsize=12)
ax.set_ylabel('Drag Coefficient', fontsize=12)
ax.set_xlim(0, 90)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g°'))
ax.tick_params(direction='in', labelsize=11)
ax.spines['top'].set_visible(False)    # 추가
ax.spines['right'].set_visible(False)  # 추가

ax.legend(fontsize=9.5, frameon=True, edgecolor='black', loc='upper left')

ax.grid(False)
plt.tight_layout()
plt.savefig('drag_coefficient_paper_style.png', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.show()