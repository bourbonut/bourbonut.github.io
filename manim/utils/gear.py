from manim import ParametricFunction, VGroup, Arc, Line, RED, YELLOW
from glm import mat3, vec3, inverse
from math import sqrt, radians, pi, cos, sin, tan, acos, asin, atan2

from .maths import u, v, O, X, Y, Z, anglebt


def involute(t, r, t0=0):
    return r * (u(t - t0) - t * v(t - t0))


def interference_curve(t, r, x, y, t0):
    return involute(t, r, t0) - x * u(t - t0) + y * v(t - t0)


def derived_involute(t, r, t0):
    return r * t * u(t - t0)


def derived_interference_curve(t, r, x, y, t0):
    return derived_involute(t, r, t0) - x * v(t - t0) - y * u(t - t0)


def jacobian_involute(rb, rp, x, y, t0):
    return lambda t1, t2: mat3(
        derived_involute(t1, rb, t0), -derived_interference_curve(t2, rp, x, y, t0), Z
    )


def profile(m, z, alpha=radians(20), ka=1, kf=1.25, interference=False):
    ha = m * ka  # adendum height
    hf = m * kf  # deddendum height
    p = pi * m  # step
    rp = m * z * 0.5
    ra = rp + ha
    rf = rp - hf
    rb = rp * cos(alpha)

    ta = sqrt(ra * ra / (rb * rb) - 1)
    tp = sqrt(rp * rp / (rb * rb) - 1)

    ha = ka * m
    la = p * 0.5 - 2 * tan(alpha) * ha  # addendum length of tooth

    if interference:
        ts = tp - atan2(tp, 1)
        phase = pi / z + 2 * (tp - atan2(tp, 1))
        phase_empty = 2 * pi / z - phase
        sl = atan2(la * 0.5, (rp - ka * m))
        f = lambda t1, t2: involute(t1, rb) - interference_curve(
            t2, rp, ka * m, 0.5 * la, phase_empty * 0.5
        )
        J = jacobian_involute(rb, rp, ka * m, 0.5 * la, phase_empty * 0.5)
        t1, t2, t3 = 0.5 * ta, -0.5 * ta, 0
        for i in range(8):
            t1, t2, t3 = vec3(t1, t2, t3) - inverse(J(t1, t2)) * f(t1, t2)

        side1 = ParametricFunction(lambda x: involute(x, rb), t_range=[t1, ta])
        interference1 = ParametricFunction(
            lambda t: interference_curve(t, rp, ka * m, 0.5 * la, phase_empty * 0.5),
            t_range=[t2, 0],
        )
        side2 = side1.copy().apply_matrix(mat3(X, -Y, Z)).rotate_about_origin(phase)
        interference2 = (
            interference1.copy().apply_matrix(mat3(X, -Y, Z)).rotate_about_origin(phase)
        )
        angle_a = ta - atan2(ta, 1)
        top = Arc(ra, angle_a, phase - 2 * angle_a)
        M = interference_curve(0, rp, ka * m, 0.5 * la, phase_empty * 0.5)
        angle_M = anglebt(M, u(0.5 * phase)) * 2
        bottom = Line(
            M, mat3(u(-2 * pi / z + angle_M), v(-2 * pi / z + angle_M), Z) * M
        )

        i = ParametricFunction(
            lambda t: interference_curve(t, rp, ka * m, 0.5 * la, phase_empty * 0.5),
            t_range=[-ta, ta],
            color=RED,
        )
        h = ParametricFunction(
            lambda t: involute(t, rp, 0.5 * phase_empty - sl),
            t_range=[-0.5 * ta, 0.5 * ta],
            color=YELLOW,
        )
        return VGroup(
            side1,
            interference1,
            side2,
            interference2,
            top,
            bottom,
            i,
            # Line(O, (rp - ka * m) * u(-0.5 * phase_empty + sl)),
            h,
        ).rotate_about_origin(-phase * 0.5)
        # return VGroup(
        #     side1, interference1, side2, interference2, top, bottom,
        # ).rotate_about_origin(-phase * 0.5)
    else:
        tb = 0
        phase = pi / z + 2 * (tp - atan2(tp, 1))
        side1 = ParametricFunction(lambda x: involute(x, rb), t_range=[tb, ta])
        side2 = side1.copy().apply_matrix(mat3(X, -Y, Z)).rotate_about_origin(phase)
        angle_a = ta - atan2(ta, 1)
        top = ParametricFunction(
            lambda t: ra * u(t), t_range=[0.96 * angle_a, 1.02 * (phase - angle_a)]
        )

        small_radius = 0.5 * (rb - rf)
        segment_after = Line((rf + small_radius) * u(phase), rb * 1.01 * u(phase))
        segment_before = Line((rf + small_radius) * u(0), rb * 1.01 * u(0))

        small_angle = asin(small_radius / rf) * 0.92
        bottom = ParametricFunction(
            lambda t: rf * u(t),
            t_range=[0.98 * (phase + small_angle), 1.02 * (2 * pi / z - small_angle)],
        )
        small_center = (rf + small_radius) * X
        arc_after = ParametricFunction(
            lambda t: small_center + small_radius * u(t), t_range=[-pi, -pi / 2]
        ).rotate_about_origin(phase + small_angle)
        arc_before = ParametricFunction(
            lambda t: small_center + small_radius * u(t), t_range=[pi / 2, pi]
        ).rotate_about_origin(-small_angle)

        return VGroup(
            side1,
            side2,
            top,
            segment_before,
            segment_after,
            arc_after,
            arc_before,
            bottom,
        ).rotate_about_origin(-phase * 0.5)
