from manim import ParametricFunction, VGroup, Arc, Line, RED, YELLOW, Circle
from glm import mat3, vec3, inverse
from math import sqrt, radians, pi, cos, sin, tan, acos, asin, atan2
from functools import partial
from collections import namedtuple

from .maths import u, v, O, X, Y, Z, anglebt, rotation
from . import rack


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


def pitch_radius(m, z):
    return m * z * 0.5


def base_radius(m, z, alpha):
    return pitch_radius(m, z) * cos(alpha)


def angle_involute(r, rb):
    return sqrt(r * r / (rb * rb) - 1)


def profile(m, z, alpha=radians(20), ka=1, kf=1.25, interference=False):
    ha = m * ka  # adendum height
    hf = m * kf  # deddendum height
    p = pi * m  # step
    rp = pitch_radius(m, z)
    ra = rp + ha
    rf = rp - hf
    rb = base_radius(m, z, alpha)

    ta = angle_involute(ra, rb)
    tp = angle_involute(rp, rb)

    duplicate = (
        lambda obj, angle: obj.copy()
        .apply_matrix(mat3(X, -Y, Z))
        .rotate_about_origin(angle)
    )

    if interference:
        la = rack.addendum_length(m, alpha, ka)
        ts = tp - atan2(tp, 1)
        phase = pi / z + 2 * (tp - atan2(tp, 1))
        phase_empty = 2 * pi / z - phase
        angle_top = ta - atan2(ta, 1)
        tmin = la * 0.5 / rp

        Functions = namedtuple("Functions", ["involute", "interference"])
        functions = Functions(
            partial(involute, r=rb),
            partial(interference_curve, r=rp, x=ha, y=0.5 * la, t0=phase_empty * 0.5),
        )

        # Newton method
        f = lambda t1, t2: functions.involute(t1) - functions.interference(t2)
        J = jacobian_involute(rb, rp, ha, 0.5 * la, phase_empty * 0.5)
        t1, t2, t3 = 0.5 * ta, -0.5 * ta, 0
        for i in range(8):
            t1, t2, t3 = vec3(t1, t2, t3) - inverse(J(t1, t2)) * f(t1, t2)

        # Involute and interference curve
        side = ParametricFunction(functions.involute, t_range=[t1, ta])
        interference = ParametricFunction(functions.interference, t_range=[t2, tmin])
        duplicated_objs = map(partial(duplicate, angle=phase), (side, interference))

        # Top and bottom
        top = Arc(ra, angle_top, phase - 2 * angle_top)
        M = functions.interference(tmin)
        angle_bottom = anglebt(M, u(0.5 * phase)) * 2
        bottom = Line(M, rotation(-2 * pi / z + angle_bottom) * M)

        return VGroup(
            side, interference, *duplicated_objs, top, bottom
        ).rotate_about_origin(-phase * 0.5)

    else:
        phase = pi / z + 2 * (tp - atan2(tp, 1))

        # Involute
        side = ParametricFunction(partial(involute, r=rb), t_range=[0, ta])

        # Arc parameters
        ArcParameters = namedtuple("ArcParameters", ["center", "radius", "angle"])
        r = 0.5 * (rb - rf)
        t = -atan2(r, rf + r)
        arcp = ArcParameters((rf + r) * u(t), r, t)
        arc = Arc(arcp.radius, -pi + arcp.angle, -pi / 2, arc_center=arcp.center)

        # Joint, top and bottom
        angle_top = ta - atan2(ta, 1)
        top = Arc(ra, angle_top, phase - 2 * angle_top)
        joint = Line(arcp.center + arcp.radius * u(-3 * pi / 2 + arcp.angle), rb * X)
        M = arcp.center + arcp.radius * u(-pi + arcp.angle)
        angle_bottom = anglebt(M, u(0.5 * phase)) * 2
        bottom = Line(M, rotation(-2 * pi / z + angle_bottom) * M)

        # Duplicated objects
        duplicated_objs = map(partial(duplicate, angle=phase), (side, arc, joint))

        return VGroup(
            side, top, arc, joint, bottom, *duplicated_objs
        ).rotate_about_origin(-phase * 0.5)
