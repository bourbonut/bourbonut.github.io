from manim import *
from manim.utils.color import Colors
from itertools import starmap
from functools import reduce
import random
import glm
from utils.maths import *

TARGET = 3 * X + 2 * Y + 3 * Z


def get_rect(w, h, origin, posx, posy, color="#ffd36a"):
    rec = Rectangle("#ffd36a", width=w, height=h)
    rec.stroke_width = 1
    rec.set_fill(color, opacity=0.5)
    rec.set_stroke(color)
    return rec.shift(origin + (0.5 * w + posx) * X + 0.5 * posy * Y)


class Context(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)
        center = Dot3D(point=O, radius=0.1, color=RED)
        dirs = [X, -X, Y, -Y, Z, -Z]
        arrows = [Arrow3D(start=d, end=2 * d, color=YELLOW, resolution=8) for d in dirs]
        self.add(axes, center, *arrows)


class Eat(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)
        dot1 = Dot3D(point=O, radius=0.1, color=RED)
        dot2 = Dot3D(point=axes.coords_to_point(*TARGET), radius=0.1, color=BLUE)
        self.add(axes, dot2, dot1)
        weights = [3, 2, 3]
        dirs = [axes.coords_to_point(*d) for d in [X, Y, Z]]
        while any(map(lambda x: x!=0, weights)):
            i = random.choices([0, 1, 2], weights=weights)[0]
            weights[i] -= 1
            direction = dirs[i]
            self.play(dot1.animate.shift(direction))
        self.wait(0.5)


class Simplified(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)
        dot1 = Dot3D(point=ORIGIN, radius=0.1, color=RED)
        dot2 = Dot3D(point=axes.coords_to_point(*X), radius=0.1, color=BLUE)
        self.add(axes, dot1, dot2)


class RewardsDot(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        axes.add(axes.get_axis_labels())
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)
        dot1 = Dot3D(point=ORIGIN, radius=0.1, color=RED)
        dot2 = Dot3D(point=axes.coords_to_point(*X), radius=0.1, color=BLUE)
        self.add(axes, dot2, dot1)
        rewards = [-1, -1, -1, 10, -1, -1]
        dirs = [axes.coords_to_point(*d) for d in [-X, -Y, -Z, X, Y, Z]]
        rots = [False, True, True, False, True, True]
        shifts = [Z, Z, -Y, Z, Z, Y]
        ratios = [1, 1, 0.5, 1, 1, 0.5]
        indices = list(range(6)) * 2
        random.shuffle(indices)
        for idx in indices:
            d = dirs[idx]
            r = rewards[idx]
            shift = shifts[idx]
            rot = rots[idx]
            a, b = ratios[idx], 0.5
            color = GREEN if r > 0 else RED
            text = (
                Text(f"{r}", color=color)
                .shift(a * d + b * shift)
                .rotate(PI / 2, axis=X)
            )
            if rot:
                text = text.rotate(PI / 2, axis=Z)
            text_line = Line(a * d + b * shift, a * d + shift)
            text.move_to(text_line.get_start())
            self.play(dot1.animate.move_to(d), MoveAlongPath(text, text_line))
            self.remove(text)
            self.play(FadeOut(dot1))
            dot1 = Dot3D(point=ORIGIN, radius=0.1, color=RED)
            self.play(FadeIn(dot1))


class Dice(Scene):
    def make_face(self, n):
        frame_width = config["frame_width"]
        frame_height = config["frame_height"]
        s = frame_height / 12
        square = RoundedRectangle(
            corner_radius=s * 0.25, width=s, height=s, fill_opacity=1, fill_color="#170343",
        )
        make_point = lambda: Circle(s * 0.05, color="#ffd36a", fill_opacity=1)
        v = 0.25 * s * (X + Y)
        w = 0.25 * s * (X - Y)
        u = 0.25 * s * X
        if n == 1:
            circle = make_point()
            return VGroup(square, circle)
        elif n == 2:
            circle1 = make_point().shift(v)
            circle2 = make_point().shift(-v)
            return VGroup(square, circle1, circle2)
        elif n == 3:
            circle1 = make_point()
            circle2 = make_point().shift(v)
            circle3 = make_point().shift(-v)
            return VGroup(square, circle1, circle2, circle3)
        elif n == 4:
            circle1 = make_point().shift(v)
            circle2 = make_point().shift(-v)
            circle3 = make_point().shift(w)
            circle4 = make_point().shift(-w)
            return VGroup(square, circle1, circle2, circle3, circle4)
        elif n == 5:
            circle1 = make_point()
            circle2 = make_point().shift(v)
            circle3 = make_point().shift(-v)
            circle4 = make_point().shift(w)
            circle5 = make_point().shift(-w)
            return VGroup(square, circle1, circle2, circle3, circle4, circle5)
        elif n == 6:
            circle1 = make_point().shift(u)
            circle2 = make_point().shift(-u)
            circle3 = make_point().shift(v)
            circle4 = make_point().shift(-v)
            circle5 = make_point().shift(w)
            circle6 = make_point().shift(-w)
            return VGroup(square, circle1, circle2, circle3, circle4, circle5, circle6)

    def construct(self):
        frame_width = config["frame_width"]
        frame_height = config["frame_height"]
        R = [10, -1, -1, -1, -1, -1]

        lineh = Line(-frame_height * Y, frame_height * Y)
        linev = Line(-frame_width * X, frame_width * X).shift(
            0.8 * 0.5 * frame_height * Y
        )
        action = Text("Actions").shift(
            -0.25 * frame_width * X + 0.9 * 0.5 * frame_height * Y
        )
        reward = Text("Rewards").shift(
            0.25 * frame_width * X + 0.9 * 0.5 * frame_height * Y
        )
        r = (0.8 + 1) / 6
        va = -0.5 * 0.8 * frame_width * X
        vr = 0.5 * 0.5 * frame_width * X
        faces = []
        action_tex = []
        proba_tex = []
        reward_tex = []
        shifta = lambda x: x.shift(va + (0.65 - i * r) * 0.5 * frame_height * Y)
        shiftr = lambda x: x.shift(vr + (0.65 - i * r) * 0.5 * frame_height * Y)
        for i in range(6):
            atex = shifta(Tex(f"$A_{i + 1} =$"))
            ptex = Tex(f"$P(A_{i + 1}) = \\frac 1 6$")
            rtex = shiftr(Tex(f"$R(A_{i + 1}) = {R[i]}$"))
            face = self.make_face(i + 1)
            face.next_to(atex, RIGHT, buff=0.25)
            ptex.next_to(face, RIGHT, buff=0.45)
            faces.append(face)
            action_tex.append(atex)
            proba_tex.append(ptex)
            reward_tex.append(rtex)
        self.add(linev, lineh)
        self.add(action, reward)
        for i in range(6):
            self.add(action_tex[i], *faces[i], proba_tex[i])
        for i in range(6):
            self.add(reward_tex[i])



class ExpectedRewardInit(Scene):
    def f1(self, x):
        if 0 <= x < 1:
            return 10 / 6
        elif 1 <= x < 6:
            return -1 / 6
        else:
            return 0

    def f2(self, x):
        if 0 <= x < 1:
            return 10 / 5
        elif 2 <= x < 6:
            return -1 / 5
        else:
            return 0

    def f3(self, x):
        if 0 <= x < 1:
            return 10 / 4
        elif 2 <= x < 5:
            return -1 / 4
        else:
            return 0

    def f4(self, x):
        if 0 <= x < 1:
            return 10 / 3
        elif 2 <= x < 3:
            return -1 / 3
        elif 4 <= x < 5:
            return -1 / 3
        else:
            return 0

    def f5(self, x):
        if 0 <= x < 1:
            return 10 / 2
        elif 2 <= x < 3:
            return -1 / 2
        else:
            return 0

    def f6(self, x):
        if 0 <= x < 1:
            return 10
        else:
            return 0

    def reduce(self, x, y):
        wunit = self.ax.get_x_unit_size()
        hunit = self.ax.get_y_unit_size()
        origin = self.ax.get_origin()
        values = [self.f[self.idx](0.5 + i) for i in range(6)]
        self.idx += 1
        y = VGroup(
            *[
                get_rect(wunit, values[i] * hunit, origin, i * wunit, values[i] * hunit)
                for i in range(6)
            ]
        )
        rects, vlabels = x
        new_vlabels = self.generate_vlabels(y, values)
        self.play(FadeOut(vlabels), ReplacementTransform(rects, y))
        self.play(FadeIn(new_vlabels))
        self.wait(0.5)
        return (y, new_vlabels)

    def generate_vlabels(self, rects, values):
        vlabels = VGroup()
        labels = self.a(6 - self.idx + 1)
        for i in range(self.idx - 1):
            labels[self.r[i]] = "$-1 \\times 0$"
        if self.idx == 5:
            labels[0] = "$10 \\times 1$"
        for rect, v, ac in zip(rects, values, labels):
            label = Tex(str(ac))
            label.font_size = 28

            pos = UP if (v > 0) else DOWN
            label.next_to(rect, pos, buff=MED_SMALL_BUFF)
            vlabels.add(label)

        return vlabels

    def construct(self):
        self.f = [self.f1, self.f2, self.f3, self.f4, self.f5, self.f6]
        self.r = [1, 5, 3, 4, 2]
        self.a = (
            lambda i: [f"$10 \\times \\frac 1 {i}$"]
            + [f"$-1 \\times \\frac 1 {i}$"] * 5
        )
        self.idx = 1
        ax = Axes(
            x_range=[-1, 7],
            y_range=[-1, 10],
            x_axis_config={"numbers_to_include": [1, 2, 3, 4, 5, 6]},
            y_axis_config={"numbers_to_include": [-1] + [i + 1 for i in range(10)]},
            tips=False,
        )
        self.ax = ax
        labels = ax.get_axis_labels()
        wunit = ax.get_x_unit_size()
        hunit = ax.get_y_unit_size()
        origin = ax.get_origin()
        values = [self.f1(0.5 + i) for i in range(6)]
        rects = VGroup(
            *[
                get_rect(wunit, values[i] * hunit, origin, i * wunit, values[i] * hunit)
                for i in range(6)
            ]
        )
        va_labels = self.a(6)
        ac_labels = [f"$A_{i + 1}$" for i in range(6)]
        alabels = VGroup()
        vlabels = VGroup()
        for rect, v, va, ac in zip(rects, values, va_labels, ac_labels):
            label = Tex(str(va))
            label.font_size = 28

            pos = UP if (v > 0) else DOWN
            label.next_to(rect, pos, buff=MED_SMALL_BUFF)
            vlabels.add(label)

            label = Tex(str(ac))
            label.font_size = 28

            pos = DOWN if (v > 0) else UP
            label.next_to(rect, pos, buff=MED_SMALL_BUFF)
            alabels.add(label)

        self.add(ax)
        self.add(rects)
        self.add(vlabels, alabels)


class ExpectedRewardAnimation(Scene):
    def f1(self, x):
        if 0 <= x < 1:
            return 10 / 6
        elif 1 <= x < 6:
            return -1 / 6
        else:
            return 0

    def f2(self, x):
        if 0 <= x < 1:
            return 10 / 5
        elif 2 <= x < 6:
            return -1 / 5
        else:
            return 0

    def f3(self, x):
        if 0 <= x < 1:
            return 10 / 4
        elif 2 <= x < 5:
            return -1 / 4
        else:
            return 0

    def f4(self, x):
        if 0 <= x < 1:
            return 10 / 3
        elif 2 <= x < 3:
            return -1 / 3
        elif 4 <= x < 5:
            return -1 / 3
        else:
            return 0

    def f5(self, x):
        if 0 <= x < 1:
            return 10 / 2
        elif 2 <= x < 3:
            return -1 / 2
        else:
            return 0

    def f6(self, x):
        if 0 <= x < 1:
            return 10
        else:
            return 0

    def reduce(self, x, y):
        wunit = self.ax.get_x_unit_size()
        hunit = self.ax.get_y_unit_size()
        origin = self.ax.get_origin()
        values = [self.f[self.idx](0.5 + i) for i in range(6)]
        self.idx += 1
        y = VGroup(
            *[
                get_rect(wunit, values[i] * hunit, origin, i * wunit, values[i] * hunit)
                for i in range(6)
            ]
        )
        rects, vlabels = x
        new_vlabels = self.generate_vlabels(y, values)
        self.play(FadeOut(vlabels), ReplacementTransform(rects, y))
        self.play(FadeIn(new_vlabels))
        self.wait(0.5)
        return (y, new_vlabels)

    def generate_vlabels(self, rects, values):
        vlabels = VGroup()
        labels = self.a(6 - self.idx + 1)
        for i in range(self.idx - 1):
            labels[self.r[i]] = "$-1 \\times 0$"
        if self.idx == 5:
            labels[0] = "$10 \\times 1$"
        for rect, v, ac in zip(rects, values, labels):
            label = Tex(str(ac))
            label.font_size = 28

            pos = UP if (v > 0) else DOWN
            label.next_to(rect, pos, buff=MED_SMALL_BUFF)
            vlabels.add(label)

        return vlabels

    def construct(self):
        self.f = [self.f1, self.f2, self.f3, self.f4, self.f5, self.f6]
        self.r = [1, 5, 3, 4, 2]
        self.a = (
            lambda i: [f"$10 \\times \\frac 1 {i}$"]
            + [f"$-1 \\times \\frac 1 {i}$"] * 5
        )
        self.idx = 1
        ax = Axes(
            x_range=[-1, 7],
            y_range=[-1, 10],
            x_axis_config={"numbers_to_include": [1, 2, 3, 4, 5, 6]},
            y_axis_config={"numbers_to_include": [-1] + [i + 1 for i in range(10)]},
            tips=False,
        )
        self.ax = ax
        labels = ax.get_axis_labels()
        wunit = ax.get_x_unit_size()
        hunit = ax.get_y_unit_size()
        origin = ax.get_origin()
        values = [self.f1(0.5 + i) for i in range(6)]
        rects = VGroup(
            *[
                get_rect(wunit, values[i] * hunit, origin, i * wunit, values[i] * hunit)
                for i in range(6)
            ]
        )
        va_labels = self.a(6)
        ac_labels = [f"$A_{i + 1}$" for i in range(6)]
        alabels = VGroup()
        vlabels = VGroup()
        for rect, v, va, ac in zip(rects, values, va_labels, ac_labels):
            label = Tex(str(va))
            label.font_size = 28

            pos = UP if (v > 0) else DOWN
            label.next_to(rect, pos, buff=MED_SMALL_BUFF)
            vlabels.add(label)

            label = Tex(str(ac))
            label.font_size = 28

            pos = DOWN if (v > 0) else UP
            label.next_to(rect, pos, buff=MED_SMALL_BUFF)
            alabels.add(label)

        self.add(ax, rects, vlabels, alabels)
        reduce(self.reduce, range(5), (rects, vlabels))


class RandomPath(ThreeDScene):
    def generate_path(self, axes):
        weights = [3, 2, 3]
        dirs = [X, Y, Z]
        lines = []
        M = O
        while any(map(lambda x: x != 0, weights)):
            i = random.choices([0, 1, 2], weights=weights)[0]
            weights[i] -= 1
            d = dirs[i]
            line = Line(
                axes.coords_to_point(*M), axes.coords_to_point(*(M + d)), color=GREEN
            )
            lines.append(line)
            M = M + d
        return lines

    def construct(self):
        axes = ThreeDAxes()
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)
        dot1 = Dot3D(point=O, radius=0.1, color=RED)
        dot2 = Dot3D(
            point=axes.coords_to_point(*(3 * X + 2 * Y + 3 * Z)), radius=0.1, color=BLUE
        )
        lines = self.generate_path(axes)
        self.play(Write(axes))
        self.play(Create(dot1), Create(dot2))
        self.begin_ambient_camera_rotation(rate=PI / 10, about="theta")
        for line in lines:
            self.play(Create(line))
        self.wait(0.5)
        self.remove(*lines)
        self.wait(0.5)
        for _ in range(8):
            lines = self.generate_path(axes)
            self.add(*lines)
            self.wait(0.5)
            self.remove(*lines)
            self.wait(0.5)
        self.stop_ambient_camera_rotation()

class Image(Scene):
    def construct(self):
        img = ImageMobject("figures/pacman.png").scale(5)
        self.add(img)
