from manim import *
from utils import *
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from utils import rack, gear, maths, fast

X = maths.X
Y = maths.Y


class RackAndGear(Scene):
    def construct(self):
        m = 1
        z = 8
        alpha = radians(20)
        ka = 1
        kf = 1.25
        rp = gear.pitch_radius(m, z)
        la = rack.addendum_length(m, alpha, ka)  # addendum length of tooth
        lf = rack.dedendum_length(m, alpha, kf)  # dedendum length of tooth
        step = rack.step(m)
        angle_step = gear.angle_step(z)
        racklines = VGroup(fast.repeat(rack.profile(m), 3, step))

        # Visualization 1
        self.add(
            fast.revolution(gear.profile(m, z, interference=True), angle_step, z)
            .rotate_about_origin(angle_step * 0.5)
            .shift(-rp * X)
        )
        self.add(
            racklines.rotate_about_origin(pi * 0.5)
            .shift(rp * X + (-step + 0.5 * lf) * Y)
            .shift(-rp * X - (step * 0.5 * Y))
        )

        # Visualization 2
        # self.add(
        #     racklines.rotate_about_origin(pi / 2)
        #     .shift(rp * X + (-step + 0.5 * lf) * Y)
        #     .shift(-rp * X)
        # )
        # self.add(
        #     fast.revolution(gear.profile(m, z, interference=False), angle_step, z).shift(
        #         -rp * X
        #     )
        # )


class Gear(Scene):
    def construct(self):
        m = 0.5
        z = 12
        ka = 1
        kf = 1.25
        alpha = radians(20)
        angle_step = gear.angle_step(z)

        profile = fast.revolution(gear.profile(m, z, alpha, ka, kf), angle_step, z)

        # Radius
        rp = gear.pitch_radius(m, z)
        rb = gear.base_radius(m, z, alpha)
        ra = gear.addendum_radius(m, z, ka)
        rf = gear.dedendum_radius(m, z, kf)

        # Circles
        circle_pitch = Circle(rp, color=RED, stroke_width=2)
        circle_base = Circle(rb, color=YELLOW, stroke_width=2)
        circle_addendum = DashedVMobject(Circle(ra, color=BLUE, stroke_width=2))
        circle_dedendum = DashedVMobject(Circle(rf, color=BLUE, stroke_width=2))

        self.add(profile, circle_pitch, circle_base, circle_addendum, circle_dedendum)


class Rack(Scene):
    def construct(self):
        m = 1
        step = rack.step(m)
        ka, kf = 1, 1.25
        ha = rack.addendum_height(m, ka)
        hf = rack.addendum_height(m, kf)
        line = lambda y, color: Line(
            y * Y, 3 * step * X + y * Y, color=color, stroke_width=2
        )
        objs = VGroup(
            fast.repeat(rack.profile(m), 3, step),
            line(0, YELLOW),
            DashedVMobject(line(ha, BLUE)),
            DashedVMobject(line(-hf, BLUE)),
        )
        self.add(objs.shift(-1.5 * step * X))
