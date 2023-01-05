from manim import *
from utils import *
from glm import *
import numpy as np
from math import radians, pi, cos, sin, tan, acos, asin, atan2
from utils import rack, gear, maths, fast
import numpy as np

X = maths.X
Y = maths.Y


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
        )
