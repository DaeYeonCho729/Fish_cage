import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# ── Parameters ──────────────────────────────────────────────
Sn   = 0.35         # 019
dw   = 0.00141       # 0.00242
Urel = 0.5
eps  = 1e-10
nu   = 1.05e-6

angles_deg = np.linspace(0, 90, 500)
angles_rad = np.deg2rad(angles_deg)
phi = np.pi / 2 - angles_rad


# ── Models ──────────────────────────────────────────────────
def bessonneau_lift(phi):
    Cn, Ct = 1.2, 0.1
    # CL = CN sin(phi) - CT cos(phi)
    return Cn * np.sin(phi) - Ct * np.cos(phi)


def lee_raw_lift(deg):
    cl = np.zeros_like(deg, dtype=float)

    cl = np.where((deg >= 0)  & (deg <= 5),  0.001831 * deg + 0.020675, cl)
    cl = np.where((deg > 5)   & (deg <= 10), 0.004817 * deg + 0.004268, cl)
    cl = np.where((deg > 10)  & (deg <= 20), 0.008873 * deg - 0.03891,  cl)
    cl = np.where((deg > 20)  & (deg <= 30), 0.008873 * deg - 0.11956,  cl)
    cl = np.where((deg > 30)  & (deg <= 43), 0.013325 * deg - 0.13140,  cl)
    cl = np.where((deg > 43)  & (deg <= 50), 0.008639 * deg + 0.070035, cl)
    cl = np.where((deg > 50)  & (deg <= 55), 0.002428 * deg + 0.383586, cl)
    cl = np.where((deg > 55)  & (deg <= 60), -0.00503 * deg + 0.797889, cl)
    cl = np.where((deg > 60)  & (deg <= 65), -0.01468 * deg + 1.38167,  cl)
    cl = np.where((deg > 65)  & (deg <= 72), -0.02828 * deg + 2.274163, cl)
    cl = np.where((deg > 72)  & (deg <= 85), -0.01127 * deg + 1.046885, cl)
    cl = np.where((deg > 85),               -0.00703 * deg + 0.689754, cl)

    return cl


def lee_smooth_lift(angles):
    breakpoints = np.array([0, 5, 10, 20, 30, 43, 50, 55, 60, 65, 72, 85, 90], dtype=float)
    cl_breaks   = lee_raw_lift(breakpoints)
    spline      = make_interp_spline(breakpoints, cl_breaks, k=3)
    return spline(angles)


def loland_lift(phi):
    # 네 drag 코드 스타일에 맞춘 대응형
    # 필요하면 네가 쓰는 정확한 Loland CL 식으로 여기만 교체
    A = (-0.05*Sn + 2.3*Sn**2 - 1.76*Sn**3)
    return A * np.sin(2 * phi)


def kristiansen_lift(phi):
    # 네 drag 코드 스타일 기준의 대응형
    # 현재는 Re 기반 drag 계산 후 각도항만 lift 형태로 분리
    Re      = (dw * Urel) / (nu * (1.0 - Sn) + eps)
    x       = np.log10(Re)

    cd_circ = (
        -78.46675 + 254.73873*x - 327.88640*x**2 + 223.64577*x**3
        - 87.92234*x**4 + 20.00769*x**5 - 2.44894*x**6 + 0.12479*x**7
    )

    cd = cd_circ * Sn * (2.0 - Sn) / (((1.0 - Sn)**2) * 2)

    cn_45 = 0.5 * cd
    y_45  = 0.25 * np.pi
    ct_45 = y_45 * (4.0 * cn_45 / (8.0 + cn_45))

    # 단순 대응형 lift 표현
    return ct_45 * np.sin(2 * phi)


cl_bessonneau  = bessonneau_lift(phi)
cl_lee         = lee_smooth_lift(angles_deg)
cl_loland      = loland_lift(phi)
cl_kristiansen = kristiansen_lift(phi)


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
cls = [cl_bessonneau, cl_lee, cl_loland, cl_kristiansen]

for cl, style in zip(cls, styles):
    ax.plot(
        angles_deg, cl,
        color='black',
        linestyle='-',
        marker=style['marker'],
        markevery=50,
        markersize=5,
        linewidth=1.2,
        label=style['label']
    )

ax.set_xlabel('Attack Angle (°)', fontsize=12)
ax.set_ylabel('Lift Coefficient', fontsize=12)
ax.set_xlim(0, 90)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g°'))
ax.tick_params(direction='in', labelsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.legend(fontsize=9.5, frameon=True, edgecolor='black', loc='upper right')

ax.grid(False)
plt.tight_layout()
plt.savefig('lift_coefficient_paper_style.png', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.show()