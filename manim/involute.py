from manim import *
from utils import *
from utils.maths import O, X, Y, Z
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from utils import rack, gear, maths, fast
from functools import partial


class InvoluteFunction(Scene):
    def construct(self):
        # Parameters
        m, z, alpha = 1, 12, radians(20)
        rb = gear.base_radius(m, z, alpha)
        self.t = 0
        rate = 0.5
        offset = -rb * X - 0.55 * rb * Y

        # Graph configuration
        axes = Axes(
            x_length=14, y_length=14, axis_config={"include_ticks": False}
        ).shift(offset)
        axes_labels = axes.get_axis_labels()

        # Dots
        on_circle = Dot(rb * maths.u(self.t), color=RED).shift(offset)
        on_involute = Dot(gear.involute(self.t, rb), color=BLUE).shift(offset)

        # Base circle
        circle = Circle(rb).shift(offset)

        # Involute curve
        self.curve = VGroup(Line(on_involute.get_center(), on_involute.get_center()))

        def go_around_circle(mob, dt):
            self.t += dt * rate
            mob.move_to(rb * maths.u(self.t)).shift(offset)

        def go_line_circle():
            return Line(O + offset, on_circle.get_center())

        def go_around_involute(mob, dt):
            mob.move_to(gear.involute(self.t, rb)).shift(offset)

        def go_line_involute():
            p1, p2 = on_circle.get_center(), on_involute.get_center()
            return Line(p1, p2)

        def draw_right_angle():
            line1, line2 = Line(O + offset, on_circle.get_center()), Line(
                on_circle.get_center(), on_involute.get_center()
            )
            inter = cross(
                line1.get_end() - line1.get_start(),
                line2.get_end() - line2.get_start(),
            )
            rightangle_function = lambda l1, l2: RightAngle(
                l1, l2, length=min(0.2, 0.6 * self.t), quadrant=(-1, 1)
            )
            class_ = Line if inter[2] == 0.0 else rightangle_function
            return class_(line1, line2)

        def draw_involute():
            last_line = self.curve[-1]
            new_line = Line(
                last_line.get_end(), gear.involute(self.t, rb) + offset, color=BLUE
            )
            self.curve.add(new_line)
            return self.curve

        on_circle.add_updater(go_around_circle)
        on_involute.add_updater(go_around_involute)

        origin2circle = always_redraw(go_line_circle)
        circle2involute = always_redraw(go_line_involute)
        rightangle = always_redraw(draw_right_angle)
        involute = always_redraw(draw_involute)

        self.add(axes, axes_labels, circle)
        self.add(on_circle, on_involute)
        self.add(origin2circle, circle2involute, rightangle, involute)
        self.wait(3.5)


class Static(Scene):
    def construct(self):
        # Parameters
        m, z, alpha = 1, 12, radians(20)
        ka = 1
        rb = gear.base_radius(m, z, alpha)
        offset = -rb * X - 0.55 * rb * Y
        tmax = gear.angle_involute(gear.addendum_radius(m, z, ka), rb) * 1.3

        # Axes
        axes = Axes(
            x_length=25, y_length=14, axis_config={"include_ticks": False}
        ).shift(offset)
        axes_labels = axes.get_axis_labels()

        # Involute, dots and circle
        involute = ParametricFunction(
            partial(gear.involute, r=rb), t_range=[0, 2 * tmax], color=BLUE
        ).shift(offset)
        pos_circle = rb * maths.u(tmax)
        pos_involute = gear.involute(tmax, rb)
        on_circle = Dot(pos_circle, color=RED).shift(offset)
        on_involute = Dot(pos_involute, color=BLUE).shift(offset)
        circle = Circle(rb, color=RED).shift(offset)

        # Lines
        origin2circle = Line(offset, on_circle.get_center())
        circle2involute = Line(on_circle.get_center(), on_involute.get_center())
        origin2involute = Line(offset, on_involute.get_center())
        extra = DashedLine(pos_circle, 1.6 * pos_circle).shift(offset)

        # Angles
        var_angle = CurvedDoubleArrow(
            1.5 * rb * X, 1.5 * pos_circle, tip_length=0.2, radius=1.5 * rb
        ).shift(offset)
        involute_angle = CurvedDoubleArrow(
            1.3 * rb * X,
            1.3 * rb * normalize(pos_involute),
            radius=1.3 * rb,
            tip_length=0.2,
        ).shift(offset)
        delta_angle = CurvedDoubleArrow(
            1.3 * rb * normalize(pos_involute),
            1.3 * pos_circle,
            radius=1.3 * rb,
            tip_length=0.2,
        ).shift(offset)
        right_angle = RightAngle(origin2circle, circle2involute, quadrant=(-1, 1))

        # Text
        r_text = MathTex("r").move_to(origin2circle.get_center() - 0.2 * (X - Y))
        rt_text = MathTex("r t").move_to(circle2involute.get_center() + 0.2 * (X + 1.5 * Y))
        t_text = MathTex("t").move_to(1.6 * rb * maths.u(0.5 * tmax) + offset)
        gamma_text = MathTex("\gamma").move_to(1.4 * rb * maths.u(0.5 * (tmax - atan2(tmax, 1))) + offset)
        delta_text = MathTex("\delta").move_to(1.4 * rb * maths.u(tmax - 0.5 * atan2(tmax, 1)) + offset)

        # Vectors
        ut = Arrow(start=pos_circle - 0.2 * maths.u(tmax), end=pos_circle + 1.5 * maths.u(tmax)).shift(offset)
        vt = Arrow(start=pos_circle - 0.2 * maths.v(tmax), end=pos_circle + 1.5 * maths.v(tmax)).shift(offset)
        ut_text = MathTex("\\overrightarrow{u(t)}").move_to(pos_circle + 1.5 * maths.u(tmax) + offset - 0.5 * (X - 0.5 * Y))
        vt_text = MathTex("\\overrightarrow{v(t)}").move_to(pos_circle + 1.7 * maths.v(tmax) + offset + 0.3 * (X + Y))

        self.add(
            axes,
            axes_labels,
            involute,
            circle,
            origin2circle,
            circle2involute,
            origin2involute,
            extra,
            var_angle,
            involute_angle,
            delta_angle,
            right_angle,
            r_text,
            rt_text,
            t_text,
            gamma_text,
            delta_text,
            ut,
            vt,
            ut_text,
            vt_text,
            on_circle,
            on_involute,
        )
