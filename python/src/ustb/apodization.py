"""Apodization computation matching MATLAB uff.apodization.

Computes pixel-dependent apodization weights for receive elements or
transmit waves, implementing Eqs. (22)-(23) of the Generalized Beamformer.
The geometry and window definitions are ported from +uff/apodization.m.
"""

import numpy as np
from ustb.enums import Wavefront, Window


class Apodization:
    """Apodization weights for receive or transmit channels.

    Mirrors MATLAB uff.apodization. For receive apodization, set probe and
    focus. For transmit apodization, set sequence and focus.
    """

    def __init__(self):
        self.probe = None
        self.focus = None
        self.sequence = None
        self.f_number = np.array([1.0, 1.0])
        self.window = Window.none
        self.MLA = np.array([1, 1])
        self.MLA_overlap = np.array([0, 0])
        self.tilt = np.array([0.0, 0.0])
        self.minimum_aperture = np.array([1e-3, 1e-3])
        self.maximum_aperture = np.array([10.0, 10.0])
        self.grating_lobe_angle = np.array([np.pi, np.pi])
        self.apodization_vector = None
        self.origin = None

    @property
    def N_elements(self):
        if self.probe is not None:
            return self.probe.N_elements
        elif self.sequence is not None:
            return len(self.sequence)
        return 0

    @property
    def data(self):
        """Compute apodization matrix [N_pixels x N_elements_or_waves]."""
        if self.focus is None:
            raise ValueError("Apodization requires focus (scan) to be set")

        if self.probe is not None:
            return self._receive_apodization()
        elif self.sequence is not None:
            return self._transmit_apodization()
        else:
            return np.ones((self._N_pixels, 1), dtype=np.float32)

    @property
    def _N_pixels(self):
        return getattr(self.focus, "N_pixels", None) or np.asarray(self.focus.x).size

    def _focus_xyz(self):
        x = np.asarray(self.focus.x, dtype=np.float64).ravel()
        z = np.asarray(self.focus.z, dtype=np.float64).ravel()
        y = getattr(self.focus, "y", None)
        y = np.zeros_like(x) if y is None else np.asarray(y, dtype=np.float64).ravel()
        return x, y, z

    @staticmethod
    def _pair(value):
        """Expand a scalar setting to [x, y], as the MATLAB setters do."""
        value = np.atleast_1d(np.asarray(value, dtype=np.float64)).ravel()
        if value.size == 1:
            return np.array([value[0], value[0]])
        return value[:2]

    # ------------------------------------------------------------------
    # Receive (aperture) apodization
    # ------------------------------------------------------------------
    def _receive_apodization(self):
        """Compute receive apodization based on probe geometry and f-number."""
        N_pixels = self._N_pixels
        N_elements = self.probe.N_elements

        if self.apodization_vector is not None:
            vector = np.asarray(self.apodization_vector)
            if vector.size == N_elements:
                apo = np.tile(vector.ravel(), (N_pixels, 1))
            else:
                apo = vector.reshape(N_pixels, N_elements)
            return apo.astype(np.float32)

        if self.window == Window.none:
            return np.ones((N_pixels, N_elements), dtype=np.float32)

        if self.window == Window.sta:
            if self.origin is None:
                raise ValueError("origin must be set to use STA apodization")
            dist = np.sqrt(
                (np.asarray(self.probe.x) - self.origin.x) ** 2
                + (np.asarray(self.probe.y) - self.origin.y) ** 2
                + (np.asarray(self.probe.z) - self.origin.z) ** 2
            ).ravel()
            closest = (dist == dist.min()).astype(np.float32)
            return np.tile(closest, (N_pixels, 1))

        tan_theta, tan_phi = self._incidence_aperture()
        f_number = self._pair(self.f_number)
        return self._apply_window(np.abs(f_number[0] * tan_theta),
                                  np.abs(f_number[1] * tan_phi))

    def _origin_per_pixel(self, origin):
        """Origin coordinates broadcastable against the pixels.

        A list of origins (one per scanline of a sector scan) is mapped to the
        pixels by azimuth index; pyuff_ustb orders sector-scan pixels with
        azimuth varying fastest.
        """
        if not isinstance(origin, (list, tuple)):
            return origin.x, origin.y, origin.z
        n_azimuth = len(origin)
        if self._N_pixels % n_azimuth != 0:
            raise ValueError("Number of origins does not match the scan azimuth axis")
        index = np.arange(self._N_pixels) % n_azimuth
        ox = np.array([o.x for o in origin])[index]
        oy = np.array([o.y for o in origin])[index]
        oz = np.array([o.z for o in origin])[index]
        return ox[:, None], oy[:, None], oz[:, None]

    def _incidence_aperture(self):
        """Tangents of the pixel-element angles (MATLAB incidence_aperture)."""
        px, py, pz = self._focus_xyz()
        ex = np.asarray(self.probe.x, dtype=np.float64).ravel()[None, :]
        ey = np.asarray(self.probe.y, dtype=np.float64).ravel()[None, :]
        ez = np.asarray(self.probe.z, dtype=np.float64).ravel()[None, :]
        px, py, pz = px[:, None], py[:, None], pz[:, None]

        curvilinear = hasattr(self.probe, "radius") or hasattr(self.probe, "radius_x")
        origin = self.origin
        if origin is None:
            if hasattr(self.probe, "radius_x"):
                origin = _Point(0.0, 0.0, -float(self.probe.radius_x))
            elif hasattr(self.probe, "radius"):
                origin = _Point(0.0, 0.0, -float(self.probe.radius))
            elif hasattr(self.focus, "azimuth_axis"):
                origin = self.focus.origin

        if curvilinear:
            element_azimuth = np.arctan2(ex - origin.x, ez - origin.z)
            pixel_azimuth = np.arctan2(px - origin.x, pz - origin.z)
            pixel_distance = np.sqrt((px - origin.x) ** 2 + (pz - origin.z) ** 2)
            x_dist = origin.z * (pixel_azimuth - element_azimuth)
            y_dist = origin.y - ey
            z_dist = pixel_distance - origin.z + 0.0 * ex
        elif origin is not None:
            # Sector scan, or a user-set origin: aperture centred at the origin
            x0, y0, z0 = self._origin_per_pixel(origin)
            pixel_distance = np.sqrt((px - x0) ** 2 + (py - y0) ** 2 + (pz - z0) ** 2)
            x_dist = ex - x0
            y_dist = ey - y0
            z_dist = pixel_distance + 0.0 * ex
        else:
            # Flat probe and linear scan: aperture centred at each pixel's [x, y]
            x_dist = px - ex
            y_dist = py - ey
            z_dist = pz - ez

        x_dist, y_dist, z_dist = np.broadcast_arrays(x_dist, y_dist, z_dist)
        tilt = self._pair(self.tilt)
        x_dist, y_dist, z_dist = _rotate_points(x_dist, y_dist, z_dist, tilt[0], tilt[1])
        zx_dist, zy_dist = self._limit_aperture(z_dist)
        with np.errstate(divide="ignore", invalid="ignore"):
            return x_dist / zx_dist, y_dist / zy_dist

    def _limit_aperture(self, z_dist):
        """Clamp the depth used for the aperture to [minimum, maximum] * f-number."""
        f_number = self._pair(self.f_number)
        min_ap = self._pair(self.minimum_aperture) * f_number
        max_ap = self._pair(self.maximum_aperture) * f_number
        limited = []
        for axis in range(2):
            zd = np.array(z_dist, dtype=np.float64, copy=True)
            zd[np.abs(z_dist) <= min_ap[axis]] = np.sign(zd[np.abs(z_dist) <= min_ap[axis]]) * min_ap[axis]
            zd[np.abs(z_dist) >= max_ap[axis]] = np.sign(zd[np.abs(z_dist) >= max_ap[axis]]) * max_ap[axis]
            limited.append(zd)
        return limited

    # ------------------------------------------------------------------
    # Transmit (wave) apodization
    # ------------------------------------------------------------------
    def _transmit_apodization(self):
        """Compute transmit apodization based on the wave sequence."""
        N_pixels = self._N_pixels
        N_waves = len(self.sequence)

        if self.apodization_vector is not None:
            vector = np.asarray(self.apodization_vector).ravel()
            if vector.size != N_waves:
                raise ValueError("apodization_vector must have one value per wave")
            return np.tile(vector, (N_pixels, 1)).astype(np.float32)

        if self.window == Window.none:
            return np.ones((N_pixels, N_waves), dtype=np.float32)

        if self.window == Window.scanline:
            return self._scanline_apodization()

        tan_theta, tan_phi = self._incidence_wave()
        f_number = self._pair(self.f_number)
        return self._apply_window(np.abs(f_number[0] * tan_theta),
                                  np.abs(f_number[1] * tan_phi))

    def _incidence_wave(self):
        """Tangents of the pixel angles seen from each wave (MATLAB incidence_wave)."""
        px, py, pz = self._focus_xyz()
        N_pixels, N_waves = px.size, len(self.sequence)
        tilt = self._pair(self.tilt)
        f_number = self._pair(self.f_number)
        tan_theta = np.zeros((N_pixels, N_waves))
        tan_phi = np.zeros((N_pixels, N_waves))

        for n, wave in enumerate(self.sequence):
            source = wave.source
            wavefront = getattr(wave.wavefront, "value", wave.wavefront)
            if int(wavefront) == int(Wavefront.plane) or np.isinf(source.distance):
                # Plane wave: one weight per wave, the same for every pixel
                tan_theta[:, n] = np.tan(source.azimuth - tilt[0])
                tan_phi[:, n] = np.tan(source.elevation - tilt[1])
                continue

            # Diverging or converging wave
            origin = wave.origin
            SP = np.column_stack([px - source.x, py - source.y, pz - source.z])
            OP = np.column_stack([px - origin.x, py - origin.y, pz - origin.z])
            OS = np.array([source.x - origin.x, source.y - origin.y, source.z - origin.z])
            if np.linalg.norm(OS) <= np.finfo(float).eps:
                # STA: the source is at the origin. MATLAB only uses the probe
                # radius here when a probe is set, which it never is for waves.
                OS = np.array([0.0, 0.0, 1.0])

            zu = OS / np.sqrt(np.sum(OS ** 2))
            yu = np.cross(zu, [1.0, 0.0, 0.0])
            xu = np.cross(zu, yu)

            z_dist = SP @ zu
            x_dist = SP @ xu
            y_dist = SP @ yu
            zx_dist, zy_dist = self._limit_aperture(z_dist)

            # Grating lobe masking
            gl = self._pair(self.grating_lobe_angle)
            eps = np.finfo(float).eps
            zx_dist[np.abs(np.arctan2(OP[:, 0], OP[:, 2]) - np.arctan2(OS[0], OS[2])) > 0.5 * gl[0]] = eps
            zy_dist[np.abs(np.arctan2(OP[:, 1], OP[:, 2]) - np.arctan2(OS[1], OS[2])) > 0.5 * gl[1]] = eps

            with np.errstate(divide="ignore", invalid="ignore"):
                tan_theta[:, n] = x_dist / zx_dist
                tan_phi[:, n] = y_dist / zy_dist

        return tan_theta, tan_phi

    def _scanline_apodization(self):
        """Scanline apodization for focused imaging."""
        N_pixels = getattr(self.focus, "N_pixels", None) or self.focus.x.size
        N_waves = len(self.sequence)

        scan_x = self.focus.x.ravel()
        scan_y = self.focus.y.ravel() if self.focus.y is not None else np.zeros_like(scan_x)

        source_az = np.array([w.source.azimuth for w in self.sequence])

        apo = np.zeros((N_pixels, N_waves), dtype=np.float32)

        if N_waves > 1:
            d_az = np.abs(np.diff(source_az)).mean()
        else:
            d_az = 1.0

        for n_wave in range(N_waves):
            az = source_az[n_wave]
            scan_az = np.arctan2(scan_x, self.focus.z.ravel() + 1e-20)
            diff = np.abs(scan_az - az)
            mask = diff <= d_az * self.MLA[0] / 2.0
            apo[mask, n_wave] = 1.0

        return apo

    # ------------------------------------------------------------------
    # Windows (MATLAB uff.apodization window methods)
    # ------------------------------------------------------------------
    def _apply_window(self, ratio_theta, ratio_phi):
        """Separable window of the ratios |F * tan(angle)|; support is ratio <= 1/2."""
        window = _WINDOWS.get(Window(int(getattr(self.window, "value", self.window))))
        if window is None:
            raise ValueError(f"Unknown apodization window: {self.window!r}")
        return (window(ratio_theta) * window(ratio_phi)).astype(np.float32)


class _Point:
    """Minimal point with x, y, z coordinates."""

    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z


def _rotate_points(x, y, z, theta, phi):
    """Port of MATLAB tools.rotate_points."""
    if abs(theta) > 0:
        x, z = x * np.cos(theta) - z * np.sin(theta), x * np.sin(theta) + z * np.cos(theta)
    if abs(phi) > 0:
        y, z = y * np.cos(phi) - z * np.sin(phi), y * np.sin(phi) + z * np.cos(phi)
    return x, y, z


def _windowed(profile):
    """Evaluate ``profile`` inside the support (ratio <= 1/2) and 0 elsewhere.

    MATLAB multiplies a 0/1 support mask with the profile, which gives NaN where
    the ratio is infinite or NaN (e.g. pixels at zero depth); here those are 0.
    """
    def window(ratio):
        ratio = np.asarray(ratio, dtype=np.float64)
        out = np.zeros_like(ratio)
        inside = ratio <= 0.5
        out[inside] = profile(ratio[inside])
        return out
    return window


def _tukey(roll):
    def profile(r):
        taper = r > 0.5 * (1 - roll)
        out = np.ones_like(r)
        out[taper] = 0.5 * (1 + np.cos(2 * np.pi / roll * (r[taper] - roll / 2 - 0.5)))
        return out
    return _windowed(profile)


_WINDOWS = {
    Window.boxcar: _windowed(np.ones_like),
    Window.hanning: _windowed(lambda r: 0.5 + 0.5 * np.cos(2 * np.pi * r)),
    Window.hamming: _windowed(lambda r: 0.53836 + 0.46164 * np.cos(2 * np.pi * r)),
    Window.tukey25: _tukey(0.25),
    Window.tukey50: _tukey(0.50),
    Window.tukey75: _tukey(0.75),
    Window.triangle: _windowed(lambda r: 1 - 2 * r),
}
