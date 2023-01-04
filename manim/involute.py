from manim import *
from utils import *
from utils.maths import O, X, Y, Z
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from utils import rack, gear, maths, fast


class InvoluteFunction(Scene):
    def construct(self):
        # Parameters
        m, z, alpha = 1, 12, radians(20)
        rb = gear.base_radius(m, z, alpha)
        self.t = 0
        rate = 0.5
        offset = -rb * X - 0.55 * rb * Y

        # Graph configuration
        axes = Axes(x_length=14, y_length=14, axis_config={"include_ticks": False}).shift(offset)
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
            new_line = Line(last_line.get_end(), gear.involute(self.t, rb) + offset, color=BLUE)
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
