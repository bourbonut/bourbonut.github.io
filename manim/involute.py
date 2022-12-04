from manim import *
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from itertools import pairwise, starmap, chain

O = np.array([0, 0, 0])
X = np.array([1, 0, 0])
Y = np.array([0, 1, 0])
Z = np.array([0, 0, 1])

rotation = lambda t: np.array([[cos(t), -sin(t), 0], [sin(t), cos(t), 0], [0, 0, 1]])


def multilines(points, closed=False):
    return VGroup(
        *chain(
            (starmap(Line, pairwise(points))),
            (Line(points[0], points[-1])) if closed else (),
        )
    )


def repeat(vgroup, times, step, direction=X):
    final_vgroup = vgroup.copy()
    for i in range(times):
        translated_vgroup = vgroup.copy().shift(direction * i * step)
        for line in translated_vgroup:
            final_vgroup += line
    return final_vgroup


def default_height(step, alpha):
    return 0.5 * 2.25 * step / pi * cos(alpha) / cos(radians(20))


def rackprofile(step, h=None, offset=0, alpha=radians(20)):
    if h is None:
        h = default_height(step, alpha)
    e = offset  # change name for convenience

    # fraction of the tooth above the primitive circle
    x = 0.5 + 2 * e / step * tan(alpha)

    return multilines(
        [
            np.array([step * x / 2 - tan(alpha) * (h - e), h - e, 0]),
            np.array([step * x / 2 - tan(alpha) * (-h - e), -h - e, 0]),
            np.array([step * (2 - x) / 2 + tan(alpha) * (-h - e), -h - e, 0]),
            np.array([step * (2 - x) / 2 + tan(alpha) * (h - e), h - e, 0]),
            np.array([step * (2 + x) / 2 - tan(alpha) * (h - e), h - e, 0]),
        ],
    )


def involute(c, t0, t):
    x = vec2(cos(t), sin(t))
    y = vec2(-x[1], x[0])
    return c * (x + y * (t0 - t))


def involuteof(c, t0, d, t):
    x = vec2(cos(t), sin(t))
    y = vec2(-x[1], x[0])
    return (c + d) * x + c * y * (t0 - t)


def involuteat(c, r):
    return sqrt((r / c) ** 2 - 1)


def angle(p):
    return atan2(p[1], p[0])


def interpol1(a, b, x):
    return (1 - x) * a + x * b


def gearprofile(step, z, h=None, offset=0, alpha=radians(20)):
    if h is None:
        h = default_height(step, alpha)
    p = step * z / (2 * pi)  # primitive circle
    c = p * cos(alpha)  # tangent circle to gliding axis

    e = offset  # change name for convenience
    # fraction of the tooth above the primitive circle
    x = 0.5 + 2 * offset / step * tan(alpha)

    o0 = angle(involute(c, 0, tan(alpha)))  # offset of contact curve
    oi = atan((h - e) / p * tan(alpha))  # offset of interference curve

    l0 = involuteat(c, p + h + e)  # interval size of contact curve

    # if the tooth height is unreachable, find the intersection between the two contact curves
    t0 = -step / (2 * p) * x - o0
    t = t0 + l0
    if involute(c, t0, t)[1] > 0:
        # Newton solver method
        for i in range(5):
            f = sin(t) + cos(t) * (t0 - t)
            J = -sin(t) * (t0 - t)
            t -= f / J
        l0 = t - t0

    # if there is an interference curve
    interference = c > p - h + e
    if interference:
        # Newton solver method
        # to compute the parameters (t1, t2) of the intersection between contact line and interference line
        t0, ti = o0, oi
        # initial state
        t1 = t0 - tan(alpha)  # put contact line on primitive
        # put interference point on base circle
        t2 = ti + sqrt(c ** 2 - (p - h + e) ** 2) / p
        for i in range(8):
            ct1, ct2 = cos(t1), cos(t2)
            st1, st2 = sin(t1), sin(t2)
            # function value
            f = (
                c * vec2(ct1, st1)
                + c * (t0 - t1) * vec2(-st1, ct1)
                + (h - e - p) * vec2(ct2, st2)
                - p * (ti - t2) * vec2(-st2, ct2)
            )
            # jacobian matrix (f partial derivatives)
            J = mat2(
                -c * (t0 - t1) * vec2(ct1, st1),
                p * (ti - t2) * vec2(ct2, st2) + (h - e) * vec2(-st2, ct2),
            )
            # iteration
            t1, t2 = vec2(t1, t2) - inverse(J) * f
        li = t2 - ti  # interval size of interference curve
        s0 = t0 - t1  # generation start of contact curve
    else:
        s0 = involuteat(c, p - h + e)

    pts = []
    # number of points to place
    n = 20

    # parameter for first side
    place = step / (2 * p) * x  # place of intersection with the primitive circle
    t0 = place + o0  # start of contact curve
    ti = place + oi  # start of interference curve
    # contact line
    for i in range(n + 1):
        t = interpol1(t0 - l0, t0 - s0, i / n)
        v = involute(c, t0, t)
        pts.append(vec3(v, 0))
    # interference line
    if interference:
        for i in range(n + 1):
            t = interpol1(ti + li, ti, i / n)
            v = involuteof(p, ti, -h + e, t)
            pts.append(vec3(v, 0))

    # parameters for second side
    place = step / (2 * p) * (2 - x)
    t0 = place - o0
    ti = place - oi
    # interference line
    if interference:
        for i in range(n + 1):
            t = interpol1(ti, ti - li, i / n)
            v = involuteof(p, ti, -h + e, t)
            pts.append(vec3(v, 0))
    # contact line
    for i in range(n + 1):
        t = interpol1(t0 + s0, t0 + l0, i / n)
        v = involute(c, t0, t)
        pts.append(vec3(v, 0))

    pts.append(angleAxis(step / p, vec3(0, 0, 1)) * pts[0])

    return multilines(map(np.array, pts))


class Test(Scene):
    def construct(self):
        wire = rackprofile(2)
        self.add(repeat(wire, 3, 2))
        # self.add(wire)


class InvoluteFunction(Scene):
    def involute(self, t):
        x = np.array([cos(t), sin(t), 0])
        y = np.array([-x[1], x[0], 0])
        return self.radius * (x + y * (self.t0 - t))

    def construct(self):
        # Parameters
        self.t0 = 0
        self.t_offset = 0
        self.radius = 1
        rate = 0.5
        origin = O

        axes = Axes(axis_config={"include_ticks": False})
        axes_labels = axes.get_axis_labels()
        dotc = Dot(X * self.radius, color=RED)
        doti = Dot(X * self.radius, color=BLUE)
        circle = Circle(self.radius)
        rot90 = rotation(pi * 0.5)

        self.curve = VGroup()
        self.curve.add(Line(doti.get_center(), doti.get_center()))

        def go_around_circle(mob, dt):
            self.t_offset += dt * rate
            mob.move_to(np.array([cos(self.t_offset), sin(self.t_offset), 0]))

        def go_line_circle():
            return Line(origin, dotc.get_center())

        def go_around_involute(mob, dt):
            mob.move_to(self.involute(self.t_offset))

        def go_line_involute():
            p1, p2 = dotc.get_center(), doti.get_center()
            line = Line(p1, p2)
            # if self.t_offset > pi * 0.25:
            #     middle = (p1 + p2) * 0.5
            #     v = p2 - p1
            #     v = v / np.linalg.norm(v)
            #     d = 0.4 * rot90 @ v
            #     return VGroup(line, MathTex(r"t \times R", font_size=36).move_to(middle + d))
            # else:
            return line

        def draw_right_angle():
            line1, line2 = Line(origin, dotc.get_center()), Line(
                dotc.get_center(), doti.get_center()
            )
            inter = np.cross(
                line1.get_end() - line1.get_start(), line2.get_end() - line2.get_start()
            )
            class_ = (
                Line
                if inter[2] == 0.0
                else lambda l1, l2: RightAngle(
                    l1, l2, length=min(0.2, 0.6 * self.t_offset), quadrant=(-1, 1)
                )
            )
            return class_(line1, line2)

        def draw_involute():
            last_line = self.curve[-1]
            new_line = Line(
                last_line.get_end(), self.involute(self.t_offset), color=BLUE
            )
            self.curve.add(new_line)
            return self.curve

        def draw_angle():
            line1, line2 = Line(origin, origin + X), Line(origin, dotc.get_center())
            inter = np.cross(
                line1.get_end() - line1.get_start(), line2.get_end() - line2.get_start()
            )
            class_ = Line if inter[2] == 0.0 else Angle
            # if self.t_offset > pi * 0.25:
            #     angle = class_(line1, line2, radius = 0.5 + 3 * SMALL_BUFF, other_angle=False)
            #     return VGroup(class_(line1, line2), MathTex(r"t", font_size=36).move_to(angle.point_from_proportion(0.5)))
            # else:
            return class_(line1, line2)

        dotc.add_updater(go_around_circle)
        doti.add_updater(go_around_involute)

        origin2circle = always_redraw(go_line_circle)
        circle2involute = always_redraw(go_line_involute)
        rightangle = always_redraw(draw_right_angle)
        angle = always_redraw(draw_angle)
        involute = always_redraw(draw_involute)

        self.add(axes, axes_labels, circle)
        self.add(dotc, doti)
        self.add(origin2circle, circle2involute, angle, rightangle, involute)
        self.wait(10)
        # func = ParametricFunction(self.involute, t_range=np.array([0., 2.*pi]), fill_opacity=0)
        # self.add(func)


class Gear(Scene):

    def construct(self):
        m = 0.5
        z = 12
        p = pi * m
        d_pitch = m * z
        d_addendum = d_pitch + 2 * m
        d_dedendum = d_pitch - 2.5 * m
        r_pitch = d_pitch * 0.5
        r_addendum = d_addendum * 0.5
        r_dedendum = d_dedendum * 0.5
        z = 12
        profile = gearprofile(p, z)
        full_profile = profile.copy()
        for i in range(12 - 1):
            full_profile += profile.copy().rotate_about_origin((i + 1) * 2 * pi / z, Z)

        circle_pitch = Circle(r_pitch, color=RED)
        circle_addendum = DashedVMobject(Circle(r_addendum, color=BLUE))
        circle_dedendum = DashedVMobject(Circle(r_dedendum, color=BLUE))
        self.add(full_profile, circle_pitch, circle_addendum, circle_dedendum)
