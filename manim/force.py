from manim import *
from utils import *
from math import radians, tan

X = maths.X
Y = maths.Y

class CurvedArrowAntiClockwise(ArcBetweenPoints):
    def __init__(self, start_point, end_point, **kwargs):
        if "tip_shape_end" in kwargs:
            kwargs["tip_shape"] = kwargs.pop("tip_shape_end")
        from manim.mobject.geometry.tips import ArrowTriangleFilledTip

        tip_shape_end = kwargs.pop("tip_shape_end", ArrowTriangleFilledTip)
        super().__init__(start_point, end_point, **kwargs)
        self.add_tip(at_start=True, tip_shape=tip_shape_end) 


class Gear2Gear(Scene):
    def construct(self):
        # Parameters
        m = 0.85
        alpha = radians(20)
        z_pinion = 8
        z_wheel = 10
        rp_pinion = gear.pitch_radius(m, z_pinion)
        rb_pinion = gear.base_radius(m, z_pinion, alpha)
        rp_wheel = gear.pitch_radius(m, z_wheel)
        rb_wheel = gear.base_radius(m, z_wheel, alpha)
        angle_step_pinion = gear.angle_step(z_pinion)
        angle_step_wheel = gear.angle_step(z_wheel)

        # Profiles
        pinion = (
            fast.revolution(
                gear.profile(m, z_pinion, alpha), angle_step_pinion, z_pinion
            )
            .shift(rp_pinion * Y)
            .rotate(-angle_step_pinion * 0.25)
        )
        wheel = (
            fast.revolution(gear.profile(m, z_wheel, alpha), angle_step_wheel, z_wheel)
            .shift(-rp_wheel * Y)
            .rotate(angle_step_wheel * 0.25)
        )

        # Dots
        center = Dot(maths.O)
        center_pinion = Dot(rp_pinion * Y, color=RED)
        center_wheel = Dot(-rp_wheel * Y, color=BLUE)

        t_pinion = Dot(rp_pinion * Y + rb_pinion * maths.u(-alpha - PI / 2), color=RED)
        t_wheel = Dot(-rp_wheel * Y + rb_wheel * maths.u(-alpha + PI / 2), color=BLUE)

        # Circles
        circle_base_pinion = Circle(rb_pinion).shift(rp_pinion * Y)
        circle_pitch_pinion = DashedVMobject(Circle(rp_pinion)).shift(rp_pinion * Y)

        circle_base_wheel = Circle(rb_wheel, color=BLUE).shift(-rp_wheel * Y)
        circle_pitch_wheel = DashedVMobject(Circle(rp_wheel, color=BLUE)).shift(
            -rp_wheel * Y
        )

        # Names of circles
        cbp_text = Text("Base circle", color=RED, font_size=24).shift(rp_pinion * Y + rb_pinion * 0.6 * X)
        cpp_text = Text("Pitch circle", color=RED, font_size=24).shift(rp_pinion * Y + rp_pinion * 1.35 * X)
        cbw_text = Text("Base circle", color=BLUE, font_size=24).shift(-rp_wheel * Y - rb_wheel * 0.6 * X)
        cpw_text = Text("Pitch circle", color=BLUE, font_size=24).shift(-rp_wheel * Y - rp_wheel * 1.35 * X)

        # Line of action
        p1 = -rp_pinion / tan(alpha) * X + rp_pinion * Y
        p2 = rp_wheel / tan(alpha) * X - rp_wheel * Y
        line_action = DashedLine(p1, p2, dash_length=1.5)

        # Alpha angle
        line_alpha = DashedLine(-rp_pinion * X, maths.O)
        alpha_angle = CurvedDoubleArrow(
            0.9 * rp_pinion * maths.u(PI - alpha),
            0.9 * rp_pinion * maths.u(-PI),
            radius=0.9 * rp_pinion,
            tip_length=0.2,
        )
        alpha_text = MathTex("\\alpha").move_to(
            alpha_angle.get_center() + rp_pinion * 0.15 * maths.u(-PI - 0.5 * alpha)
        )

        # Names of dots
        Op = MathTex("O_1", font_size=32, color=RED).next_to(center_pinion, DOWN)
        Ow = MathTex("O_2", font_size=32, color=BLUE).next_to(center_wheel, UP)
        Tp = (
            MathTex("T_1", font_size=32, color=RED)
            .next_to(t_pinion, UP)
            .shift(-0.05 * rp_pinion * X)
        )
        Tw = (
            MathTex("T_2", font_size=32, color=BLUE)
            .next_to(t_wheel, DOWN)
            .shift(0.05 * rp_wheel * X)
        )
        P = MathTex("P", font_size=32).next_to(center, UP).shift(-0.05 * rp_pinion * X)

        # Force
        F = Arrow(start=-0.02 * p1, end=0.5 * p1, stroke_width=8)
        F_text = MathTex("\\overrightarrow{F_{2 / 1}}").move_to(0.5 * p1 + 0.15 * rp_pinion * Y)

        # Movement orientation
        pinion_m = CurvedArrow(1.25 * rp_pinion * maths.u(-PI * 0.25), 1.25 * rp_pinion * maths.u(-PI * 0.25 + alpha), radius = 1.25 * rp_pinion, color=RED).shift(rp_pinion * Y)
        wheel_m = CurvedArrowAntiClockwise(1.25 * rp_wheel * maths.u(PI * 0.25 - alpha), 1.25 * rp_wheel * maths.u(PI * 0.25), radius = 1.25 * rp_wheel, color=BLUE).shift(-rp_wheel * Y)

        vg = VGroup(
            pinion,
            wheel,
            circle_base_pinion,
            circle_pitch_pinion,
            circle_base_wheel,
            circle_pitch_wheel,
            line_action,
            line_alpha,
            alpha_angle,
            alpha_text,
            center_pinion,
            center_wheel,
            cbp_text,
            cpp_text,
            cbw_text,
            cpw_text,
            Tp,
            Tw,
            P,
            Op,
            Ow,
            F,
            F_text,
            pinion_m,
            wheel_m,
            t_pinion,
            t_wheel,
            center,
        ).shift(-(rp_pinion - rp_wheel) * 0.5 * Y)
        self.add(vg)
