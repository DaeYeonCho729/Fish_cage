import numpy as np

gravity = 9.81


class LinearRegularWave:
    """
    DNVGL-RP-C205 기반 regular linear wave model

    입력:
        H : wave height [m]
        T : wave period [s]
        h : water depth [m]
        direction_deg : wave propagation direction [deg]
        phase_deg : initial phase [deg]
        z0 : still water level [m]

    내부 계산:
        omega = 2*pi/T
        L = DNV approximation(T, h)
        k = 2*pi/L
    """

    def __init__(self, H, T, h, direction_deg=0.0, phase_deg=0.0, z0=0.0):
        self.H = float(H)
        self.T = float(T)
        self.h = float(h)
        self.direction_deg = float(direction_deg)
        self.phase_deg = float(phase_deg)
        self.z0 = float(z0)

        if self.T <= 0.0:
            raise ValueError("wave period T must be positive.")
        if self.h <= 0.0:
            raise ValueError("water depth h must be positive.")

        self.a = 0.5 * self.H
        self.omega = 2.0 * np.pi / self.T
        self.phase = np.deg2rad(self.phase_deg)

        beta = np.deg2rad(self.direction_deg)
        self.ex = np.array([np.cos(beta), np.sin(beta), 0.0], dtype=float)

        # 1) T, h 로부터 DNV 근사식으로 wavelength 계산
        self.wavelength = self._dnv_wavelength(self.T, self.h)

        # 2) 계산된 wavelength 로부터 wave number 계산
        self.k = 2.0 * np.pi / self.wavelength

        # 3) phase speed
        self.c = self.wavelength / self.T

    def _dnv_wavelength(self, T, h):
        """
        DNVGL-RP-C205의 finite depth wavelength approximation 형태 사용

        deep water wavelength:
            L0 = g T^2 / (2*pi)

        finite depth correction:
            c^2/(g h) = [ x + 1/(1 + a1 x + a2 x^2 + a3 x^3 + a4 x^4) ]^-1
            where x = omega^2 h / g

        이후
            c = sqrt(g h * c^2/(g h))
            L = c T
        """
        omega = 2.0 * np.pi / T
        x = (omega ** 2) * h / gravity

        a1 = 0.666
        a2 = 0.445
        a3 = -0.105
        a4 = 0.272

        poly = 1.0 + a1 * x + a2 * (x ** 2) + a3 * (x ** 3) + a4 * (x ** 4)
        poly = max(poly, 1e-12)

        c2_over_gh = (x + 1.0 / poly) ** (-1)
        c = np.sqrt(gravity * h * c2_over_gh)
        L = c * T

        return float(max(L, 1e-8))

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

        eta = self.z0 + self.a * np.cos(theta)
        eta = np.where(np.isfinite(eta), eta, self.z0)
        return eta

    def velocity(self, points, t):
        """
        Airy wave particle velocity
        반환 shape: (N, 3)

        자유수면 위(z > eta) 노드는 velocity = 0
        """
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts.reshape(1, 3)

        x_proj = pts[:, 0] * self.ex[0] + pts[:, 1] * self.ex[1]
        z = pts[:, 2] - self.z0

        theta = self.k * x_proj - self.omega * float(t) + self.phase
        kh = np.clip(self.k * self.h, -50.0, 50.0)

        arg = np.clip(self.k * (z + self.h), -50.0, 50.0)

        denom = np.sinh(kh)
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

        vel = np.where(np.isfinite(vel), vel, 0.0)

        eta = self.eta(pts, t)
        above_surface = pts[:, 2] > eta
        vel[above_surface] = 0.0

        return vel