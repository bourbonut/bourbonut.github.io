from manim import *
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from itertools import pairwise, starmap, chain

O = vec3(0, 0, 0)
X = vec3(1, 0, 0)
Y = vec3(0, 1, 0)
Z = vec3(0, 0, 1)

rotation = lambda t: mat3([cos(t), -sin(t), 0, [sin(t), cos(t), 0], [0, 0, 1]])


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


def revolution(vgroup, angle, times):
    final_vgroup = vgroup.copy()
    for i in range(times):
        rotated_vgroup = vgroup.copy().rotate_about_origin(angle * i)
        for line in rotated_vgroup:
            final_vgroup += line
    return final_vgroup


def rackprofile(m, alpha=radians(20), ka=1, kf=1.25):
    ha = m * ka  # addendum height
    hf = m * kf  # dedendum height
    p = pi * m  # step
    la = p * 0.5 - 2 * tan(alpha) * ha  # addendum length of tooth
    lf = p * 0.5 - 2 * tan(alpha) * hf  # dedendum length of tooth

    #    B ____ C
    #     /    \
    #    /      \
    # A /      D \____ E

    return multilines(
        [
            vec3(-tan(alpha) * hf, -hf, 0),  # A
            vec3(tan(alpha) * ha, ha, 0),  # B
            vec3(tan(alpha) * ha + la, ha, 0),  # C
            vec3(p * 0.5 + tan(alpha) * hf, -hf, 0),  # D
            vec3(p * 0.5 + lf + tan(alpha) * hf, -hf, 0),  # E
        ],
    )


def anglebt(a, b):
    return acos(dot(a, b) / (length(a) * length(b)))


def u(t):
    return vec3(cos(t), sin(t), 0)


def v(t):
    return vec3(-sin(t), cos(t), 0)

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


def gearprofile(m, z, alpha=radians(20), ka=1, kf=1.25, interference=False):
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
        phase = PI / z + 2 * (tp - atan2(tp, 1))
        phase_empty = 2 * PI / z - phase
        # sl = asin(la * 0.5 / rb)
        # sl = atan2(la * 0.5 , rb)
        # sl = atan2(la * 0.5 , rp) 
        sl = atan2(la * 0.5, (rp - ka * m))
        print(rb, rp, 0.5 * phase, la * 0.5)
        f = lambda t1, t2: involute(t1, rb) - interference_curve(t2, rp, ka * m, 0.5 * la, phase_empty * 0.5)
        J = jacobian_involute(rb, rp, ka * m, 0.5 * la, phase_empty * 0.5)
        t1, t2, t3 = 0.5 * ta, -0.5 * ta, 0
        for i in range(2):
            t1, t2, t3 = vec3(t1, t2, t3) - inverse(J(t1, t2)) * f(t1, t2)

        side1 = ParametricFunction(lambda x: involute(x, rb), t_range=[t1, ta])
        interference1 = ParametricFunction(
            lambda t: interference_curve(t, rp, ka * m, 0.5 * la, phase_empty * 0.5), t_range=[t2, 0]
        )
        side2 = side1.copy().apply_matrix(mat3(X, -Y, Z)).rotate_about_origin(phase)
        interference2 = (
            interference1.copy().apply_matrix(mat3(X, -Y, Z)).rotate_about_origin(phase)
        )
        angle_a = ta - atan2(ta, 1)
        top = Arc(ra, angle_a, phase - 2 * angle_a)
        M = interference_curve(0, rp, ka * m, 0.5 * la, phase_empty * 0.5)
        angle_M = anglebt(M, u(0.5 * phase)) * 2
        bottom = Line(M, mat3(u(-2 * PI / z + angle_M), v(-2 * PI / z + angle_M), Z) * M)

        i = ParametricFunction(lambda t: interference_curve(t, rp, ka * m, 0.5 * la, phase_empty * 0.5), t_range=[-ta, ta], color=RED)
        h = ParametricFunction(lambda t: involute(t, rp, 0.5 * phase_empty - sl), t_range=[-0.5 * ta, 0.5 * ta], color=YELLOW)
        return VGroup(
            side1, interference1, side2, interference2, top, bottom, i, Line(O, (rp - ka * m) * u(-0.5 * phase_empty + sl)), h, 
        ).rotate_about_origin(-phase * 0.5)
        # return VGroup(
        #     side1, interference1, side2, interference2, top, bottom,
        # ).rotate_about_origin(-phase * 0.5)
    else:
        tb = 0
        phase = PI / z + 2 * (tp - atan2(tp, 1))
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
            t_range=[0.98 * (phase + small_angle), 1.02 * (2 * PI / z - small_angle)],
        )
        small_center = (rf + small_radius) * X
        arc_after = ParametricFunction(
            lambda t: small_center + small_radius * u(t), t_range=[-PI, -PI / 2]
        ).rotate_about_origin(phase + small_angle)
        arc_before = ParametricFunction(
            lambda t: small_center + small_radius * u(t), t_range=[PI / 2, PI]
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


class Test(Scene):
    def construct(self):
        m = 1
        z = 12
        rp = m * z * 0.5
        l = -tan(radians(20)) * (1.25 * m)
        la = pi * m * 0.5 - 2 * tan(radians(20)) * 1 * m  # addendum length of tooth
        lf = pi * m * 0.5 - 2 * tan(radians(20)) * 1.25 * m  # dedendum length of tooth
        wire = rackprofile(m)
        g = VGroup(repeat(wire, 3, PI * m)) 
        # self.add(g.rotate_about_origin(PI / 2).shift(rp * X + (- PI * m - l + 0.5 * lf) * Y).shift(-rp * X))
        self.add(revolution(gearprofile(m, z, interference=True), 2 * PI / z, z).rotate_about_origin(PI / z).shift(-rp * X))
        self.add(g.rotate_about_origin(PI / 2).shift(rp * X + (- PI * m - l + 0.5 * lf) * Y).shift(-rp * X - (PI * m / 2 * Y)))
        # self.add(revolution(gearprofile(m, z, interference=True), 2 * PI / z, z).shift(-rp * X))
        # self.add(Line(-10 * X, 10 *X))


class InvoluteFunction(Scene):
    def involute(self, t):
        x = vec3(cos(t), sin(t), 0)
        y = vec3(-x[1, x[0], 0])
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
            mob.move_to(vec3(cos(self.t_offset), sin(self.t_offset), 0))

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
        # func = ParametricFunction(self.involute, t_range=vec3(0., 2.*pi), fill_opacity=0)
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


class Rack(Scene):
    def construct(self):
        m = 0.5
        self.add(rackprofile(m))
