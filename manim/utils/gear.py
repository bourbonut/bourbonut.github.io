from manim import ParametricFunction, VGroup, Arc, Line, RED, YELLOW, Dot
from glm import mat3, vec3, inverse
from math import sqrt, radians, pi, cos, atan2
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


def angle_step(z):
    return 2 * pi / z


def pitch_radius(m, z):
    return m * z * 0.5


def base_radius(m, z, alpha):
    return pitch_radius(m, z) * cos(alpha)


def addendum_radius(m, z, ka):
    return pitch_radius(m, z) + ka * m


def dedendum_radius(m, z, kf):
    return pitch_radius(m, z) - kf * m


def angle_involute(r, rb):
    return sqrt(r * r / (rb * rb) - 1)


def profile(m, z, alpha=radians(20), ka=1, kf=1.25, interference=True):
    # Parameters
    ha = m * ka  #                  addendum height
    hf = m * kf  #                  dedendum height
    p = pi * m  #                   step
    rp = pitch_radius(m, z) #       pitch radius
    ra = rp + ha #                  addendum radius
    rf = rp - hf #                  dedendum radius
    rb = base_radius(m, z, alpha) # base radius

    ta = angle_involute(ra, rb) #   addendum angle
    tp = angle_involute(rp, rb) #   pitch angle

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

        # Top and bottom
        top = Arc(ra, angle_top, phase - 2 * angle_top)
        M = functions.interference(tmin)
        angle_bottom = anglebt(M, u(0.5 * phase)) * 2
        bottom = Line(M, rotation(-2 * pi / z + angle_bottom) * M)

        # Patches
        top_dot = Dot(ra * u(angle_top), radius=0.02)
        interference_dot = Dot(functions.interference(t2), radius=0.02)
        involute_dot = Dot(functions.involute(t1), radius=0.02)
        bottom_dot = Dot(M, radius=0.02)
        dots = (top_dot, interference_dot, involute_dot, bottom_dot)

        # Duplicated objects
        duplicated_objs = map(
            partial(duplicate, angle=phase), (side, interference) + dots
        )

        return VGroup(
            side, interference, *duplicated_objs, top, bottom, *dots
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

        # Patches
        top_dot = Dot(ra * u(angle_top), radius=0.02)
        side_dot = Dot(rb * X, radius=0.02)
        point = arcp.center + arcp.radius * u(-3 * pi / 2 + arcp.angle)
        joint_dot = Dot(point, radius=0.02)
        bottom_dot = Dot(M, radius=0.02)
        dots = (top_dot, side_dot, joint_dot, bottom_dot)

        # Duplicated objects
        duplicated_objs = map(
            partial(duplicate, angle=phase), (side, arc, joint) + dots
        )

        return VGroup(
            side, top, arc, joint, bottom, *duplicated_objs, *dots
        ).rotate_about_origin(-phase * 0.5)
