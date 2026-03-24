import numpy as np

gravity = 9.81


class LinearRegularWave:
    # Airy wave theory
    # 입력:
    # H : wave height
    # T : wave period
    # h : water depth
    # wavelength : wave length
    # z0 : still water level

    def __init__(self, H, T, h, wavelength, direction_deg=0.0, phase_deg=0.0, z0=0.0):
        self.H = float(H)
        self.T = float(T)
        self.h = float(h)
        self.wavelength = float(wavelength)
        self.direction_deg = float(direction_deg)
        self.phase_deg = float(phase_deg)
        self.z0 = float(z0)

        self.a = 0.5 * self.H

        # 교과서 정의 그대로
        self.omega = 2.0 * np.pi / self.T
        self.k = 2.0 * np.pi / self.wavelength

        self.phase = np.deg2rad(self.phase_deg)

        beta = np.deg2rad(self.direction_deg)
        self.ex = np.array([np.cos(beta), np.sin(beta), 0.0], dtype=float)

    def eta(self, points, t):
        """
        자유수면 변위
        반환 shape: (N,)
        """
        pts = np.asarray(points, dtype=float)

        if pts.ndim == 1:
            pts = pts.reshape(1, 3)

        x_proj = pts[:, 0] * self.ex[0] + pts[:, 1] * self.ex[1]
        theta = self.k * x_proj - self.omega * float(t) + self.phase

        eta = self.a * np.cos(theta)
        eta = np.where(np.isfinite(eta), eta, 0.0)
        return eta

    def velocity(self, points, t):
        """
        Airy wave particle velocity
        반환 shape: (N, 3)

        Hui 방식 반영:
        - 자유수면 위(z > eta) 노드는 velocity = 0
        """
        pts = np.asarray(points, dtype=float)

        if pts.ndim == 1:
            pts = pts.reshape(1, 3)

        x_proj = pts[:, 0] * self.ex[0] + pts[:, 1] * self.ex[1]
        z = pts[:, 2] - self.z0

        theta = self.k * x_proj - self.omega * float(t) + self.phase
        kh = self.k * self.h

        # overflow 방어
        arg = self.k * (z + self.h)
        arg = np.clip(arg, -50.0, 50.0)
        kh_clip = np.clip(kh, -50.0, 50.0)

        denom = np.sinh(kh_clip)
        if abs(denom) < 1e-12:
            denom = 1e-12

        c1 = self.a * self.omega * np.cosh(arg) / denom
        c2 = self.a * self.omega * np.sinh(arg) / denom

        ux = c1 * np.cos(theta)
        uz = c2 * np.sin(theta)

        vel = np.zeros((len(pts), 3), dtype=float)
        vel[:, 0] = ux * self.ex[0]
        vel[:, 1] = ux * self.ex[1]
        vel[:, 2] = uz

        # NumPy 1.16.1 호환
        vel = np.where(np.isfinite(vel), vel, 0.0)

        # 자유수면 위는 0
        eta = self.eta(pts, t)
        above_surface = pts[:, 2] > eta
        vel[above_surface] = 0.0

        return vel