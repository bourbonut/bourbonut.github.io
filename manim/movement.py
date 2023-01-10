from functools import partial
from manim import *
from utils import *
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from utils import rack, gear, maths, fast
import numpy as np

X = maths.X
Y = maths.Y


class CustomMoveAlongPath(Animation):
    def __init__(
        self,
        mobject: Mobject,
        path: callable,
        parameter: float,
        suspend_mobject_updating: bool | None = False,
        **kwargs,
    ) -> None:
        super().__init__(
            mobject, suspend_mobject_updating=suspend_mobject_updating, **kwargs
        )
        self.path = path
        self.original = self.mobject.copy()
        self.parameter = parameter

    def interpolate_mobject(self, alpha: float) -> None:
        t = self.rate_func(alpha) * self.parameter
        angle = maths.anglebt(self.path(t), maths.X)
        self.mobject.become(self.original.copy().rotate(t, about_point=ORIGIN).shift(self.path(t)))


class GearOnRack(Scene):
    def construct(self):
        m, z, alpha = 1, 8, radians(20)
        ka = 1
        kf = 1.25
        rp = gear.pitch_radius(m, z)
        la = rack.addendum_length(m, alpha, ka)  # addendum length of tooth
        lf = rack.dedendum_length(m, alpha, kf)  # dedendum length of tooth
        step = rack.step(m)
        angle_step = gear.angle_step(z)

        gearprofile = fast.revolution(
            gear.profile(m, z, alpha, interference=True), angle_step, z
        )
        rackprofile = fast.repeat(rack.profile(m, alpha), 2 * z, step)
        rackprofile = rackprofile.rotate_about_origin(pi * 0.5).shift(
            (-(0.5 + z) * step + 0.5 * lf) * Y
        )
        gearprofile = gearprofile.rotate_about_origin(angle_step * 0.5).shift(-rp * X)

        self.play(
            Rotate(gearprofile, about_point=-rp * X, angle=angle_step),
            rackprofile.animate.shift(step * Y),
            run_time=10,
            run_func=linear,
        )


class RackOnGear(Scene):
    def construct(self):
        m, z, alpha = 0.5, 8, radians(20)
        ka = 1
        kf = 1.25
        rp = gear.pitch_radius(m, z)
        la = rack.addendum_length(m, alpha, ka)  # addendum length of tooth
        lf = rack.dedendum_length(m, alpha, kf)  # dedendum length of tooth
        step = rack.step(m)
        angle_step = gear.angle_step(z)

        rp = gear.pitch_radius(m, z)
        inv = ParametricFunction(partial(gear.involute, r=rp), t_range=[0, pi], color=BLUE)

        gearprofile = fast.revolution(
            gear.profile(m, z, alpha, interference=True), angle_step, z
        )
        gearprofile = gearprofile.rotate_about_origin(angle_step * 0.5)

        rackprofile = fast.repeat(rack.profile(m, alpha), 2 * z, step)
        rackprofile = rackprofile.rotate_about_origin(pi * 0.5).shift(
            (-(0.5 + z) * step + 0.5 * lf) * Y
        )
        vg = VGroup(Dot(ORIGIN), rackprofile)

        self.add(gearprofile)
        self.add(inv)
        self.play(CustomMoveAlongPath(vg, partial(gear.involute, r=rp), pi), run_time=5)
