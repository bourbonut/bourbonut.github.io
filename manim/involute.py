from re import I
from manim import *
from utils import *
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from utils import rack, gear, maths, fast

X = maths.X
Y = maths.Y


class Test(Scene):
    def construct(self):
        m = 1
        z = 12
        rp = m * z * 0.5
        alpha = radians(20)
        ka = 1
        kf = 1.25
        la = rack.addendum_length(m, alpha, ka)  # addendum length of tooth
        lf = rack.dedendum_length(m, alpha, kf)  # dedendum length of tooth
        racklines = VGroup(fast.repeat(rack.profile(m), 3, PI * m))
        # self.add(g.rotate_about_origin(PI / 2).shift(rp * X + (- PI * m - l + 0.5 * lf) * Y).shift(-rp * X))
        self.add(
            fast.revolution(gear.profile(m, z, interference=True), 2 * PI / z, z)
            .rotate_about_origin(PI / z)
            .shift(-rp * X)
        )
        self.add(
            racklines.rotate_about_origin(PI / 2)
            .shift(rp * X + (-PI * m + 0.5 * lf) * Y)
            .shift(-rp * X - (PI * m / 2 * Y))
        )
        # self.add(revolution(gearprofile(m, z, interference=True), 2 * PI / z, z).shift(-rp * X))


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
