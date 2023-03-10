from manim import *
from manim.mobject.geometry.tips import ArrowTriangleFilledTip
from utils import *
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from utils import rack, gear, maths, fast

O = maths.O
X = maths.X
Y = maths.Y


def CustomDoubleArrow(start_point, end_point, **kwargs):
    line1 = Line(start_point, end_point)
    line2 = Line(end_point, start_point)
    line1.add_tip(**kwargs)
    line2.add_tip(**kwargs)
    return VGroup(line1, line2)


class RackAndGear(Scene):
    def construct(self):
        # Parameters
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

        racklines = fast.repeat(rack.profile(m), 3, step)

        # Pitch characteristics
        pitch_circle = DashedVMobject(Circle(rp, color=YELLOW))
        pitch_line = DashedLine(O, 3 * step * X, color=YELLOW, dash_length=1)

        self.add(
            fast.revolution(gear.profile(m, z, interference=False), angle_step, z)
            .rotate_about_origin(angle_step * 0.5)
            .shift(-0.5 * rp * X),
            pitch_circle.shift(-0.5 * rp * X),
            Dot(-0.5 * rp * X),
        )
        offset_y = -1.5 * step + 0.5 * lf
        self.add(
            racklines.rotate_about_origin(pi * 0.5).shift(offset_y * Y + 0.5 * rp * X),
            pitch_line.rotate_about_origin(pi * 0.5).shift(offset_y * Y + 0.5 * rp * X),
        )


def gear_information(interference=False):
    m = 1
    z = 12
    ka = 1
    kf = 1.25
    alpha = radians(20)
    angle_step = gear.angle_step(z)

    profile = fast.revolution(
        gear.profile(m, z, alpha, ka, kf, interference), angle_step, z
    )

    # Radius
    rp = gear.pitch_radius(m, z)
    rb = gear.base_radius(m, z, alpha)
    ra = gear.addendum_radius(m, z, ka)
    rf = gear.dedendum_radius(m, z, kf)

    # Circles
    circle_pitch = Circle(rp, color=YELLOW, stroke_width=2)
    circle_base = Circle(rb, color=RED, stroke_width=2)
    circle_addendum = DashedVMobject(Circle(ra, color=BLUE, stroke_width=2))
    circle_dedendum = DashedVMobject(Circle(rf, color=BLUE, stroke_width=2))

    # Center
    center = Dot(O)

    # Arrows
    pitch_arrow = Line(O, rp * maths.u(PI / 3), color=YELLOW)
    base_arrow = Line(O, rb * maths.u(2 * PI / 3), color=RED)
    addendum_arrow = Line(O, ra * maths.u(PI / 6), color=BLUE)
    dedendum_arrow = Line(O, rf * maths.u(PI - PI / 6), color=BLUE)

    pitch_arrow.add_tip(tip_shape=ArrowTriangleFilledTip)
    base_arrow.add_tip(tip_shape=ArrowTriangleFilledTip)
    addendum_arrow.add_tip(tip_shape=ArrowTriangleFilledTip)
    dedendum_arrow.add_tip(tip_shape=ArrowTriangleFilledTip)

    # Text
    rp_text = MathTex("r_p", color=YELLOW).move_to(
        pitch_arrow.get_center() + 0.1 * ra * maths.v(PI / 3)
    )
    rb_text = MathTex("r_b", color=RED).move_to(
        base_arrow.get_center() + 0.1 * ra * maths.v(2 * PI / 3)
    )
    ra_text = MathTex("r_a", color=BLUE).move_to(
        addendum_arrow.get_center() + 0.1 * ra * maths.v(PI / 6)
    )
    rf_text = MathTex("r_f", color=BLUE).move_to(
        dedendum_arrow.get_center() + 0.1 * ra * maths.v(PI - PI / 6)
    )

    objs = VGroup(
        profile,
        circle_pitch,
        circle_base,
        circle_addendum,
        circle_dedendum,
        pitch_arrow,
        base_arrow,
        addendum_arrow,
        dedendum_arrow,
        rp_text,
        rb_text,
        ra_text,
        rf_text,
        center,
    )
    return objs.shift(-0.55 * rp * Y)


class GearWithInterference(Scene):
    def construct(self):
        self.add(gear_information(True))


class GearWithoutInterference(Scene):
    def construct(self):
        self.add(gear_information(False))


class Rack(Scene):
    def construct(self):
        m = 1
        step = rack.step(m)
        alpha = radians(20)
        ka, kf = 1, 1.25
        ha = rack.addendum_height(m, ka)
        hf = rack.addendum_height(m, kf)
        la = rack.addendum_length(m, alpha, ka)
        line = lambda y, color: Line(
            y * Y, 3 * step * X + y * Y, color=color, stroke_width=2
        )
        A = step * X - hf * Y
        C = (step + (ha + hf) * tan(alpha) + 0.5 * la) * X
        M = C + ha * Y
        P = M + 0.5 * la * X
        extra = -0.5 * step * X

        # Rack
        rackprofile = fast.repeat(rack.profile(m), 3, step)

        # Lines
        pitch_line = line(0, YELLOW)
        addedumn_line = DashedVMobject(line(ha, BLUE))
        dedendum_line = DashedVMobject(line(-hf, BLUE))
        extra_lines = [
            DashedVMobject(Line(extra, O)),
            DashedVMobject(Line(extra + ha * Y, ha * Y)),
            DashedVMobject(Line(extra - hf * Y, -hf * Y)),
        ]
        alpha_line1 = DashedVMobject(
            Line(
                A + (ha + hf) * maths.u(PI / 2 - alpha),
                A + (3 * ha + hf) * maths.u(PI / 2 - alpha),
            )
        )
        alpha_line2 = DashedVMobject(
            Line(A, A + (3 * ha + hf) * sin(PI / 2 - alpha) * Y)
        )

        # Dots
        center = Dot(C)
        middle = Dot(M)
        corner = Dot(P)

        # Double arrows
        # Custom class because Manim has tips bugged
        ha_arrow = CustomDoubleArrow(0.5 * extra, 0.5 * extra + ha * Y, tip_length=0.2, tip_width=0.2)
        hf_arrow = CustomDoubleArrow(0.5 * extra, 0.5 * extra - hf * Y, tip_length=0.2, tip_width=0.2)
        alpha_arrow = CurvedDoubleArrow(
            A + (2 * ha + hf) * maths.u(PI / 2 - alpha),
            A + (2 * ha + hf) * Y,
            tip_length=0.2,
            radius=2 * ha + hf,
        )
        X_arrow = CustomDoubleArrow(C, M, tip_width=0.15, tip_length=0.15)
        Y_arrow = CustomDoubleArrow(M, P, tip_width=0.15, tip_length=0.15)

        # Text
        ha_text = MathTex("h_a").move_to(-0.4 * step * X + 0.5 * ha * Y)
        hf_text = MathTex("h_f").move_to(-0.4 * step * X - 0.5 * hf * Y)
        alpha_text = MathTex("\\alpha").move_to(
            A + (2.5 * ha + hf) * maths.u(PI / 2 - 0.5 * alpha)
        )
        X_text = MathTex("X", font_size=28).move_to((C + P) * 0.5)
        Y_text = MathTex("Y", font_size=28).move_to((C + P) * 0.5 * (X + 2.5 * Y))

        objs = VGroup(
            rackprofile,
            pitch_line,
            addedumn_line,
            dedendum_line,
            alpha_line1,
            alpha_line2,
            *extra_lines,
            center,
            middle,
            corner,
            ha_arrow,
            hf_arrow,
            alpha_arrow,
            X_arrow,
            Y_arrow,
            ha_text,
            hf_text,
            alpha_text,
            X_text,
            Y_text,
        )
        self.add(objs.shift(-1.25 * step * X))
