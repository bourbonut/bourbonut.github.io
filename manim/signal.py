from manim import *
from glm import *
from math import exp
from random import random
from functools import reduce

frame_width = config["frame_width"]
frame_height = config["frame_height"]
MYYELLOW = "#ffd36a"

T_START = 2
T_END = 2.1

O = vec3(0)
X = vec3(1, 0, 0)
Y = vec3(0, 1, 0)
Z = vec3(0, 0, 1)


def neural_network(radius=0.1):
    r = radius
    circle_ref = lambda color: Circle(r, color=color, fill_opacity=0.5)
    dy = frame_height * r * 0.5
    dx = dy * 1.5
    dirs = [O, dy * Y, -dy * Y, dx * X + 0.5 * dy * Y, dx * X - 0.5 * dy * Y]
    circles = [
        circle_ref(color=(BLUE if i < 3 else RED)).shift(d) for i, d in enumerate(dirs)
    ]
    connections = [(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)]
    lines = [Line(dirs[a] + r * X, dirs[b] - r * X) for a, b in connections]
    barycenter = reduce(add, dirs) / len(dirs)
    return VGroup(*circles, *lines).shift(-barycenter)


class Impulse(Scene):
    def construct(self):
        axes = Axes(
            x_range=[0, 7, 0.2],
            y_range=[0, 6, 0.5],
            tips=False,
        ).shift(0.3 * X)
        f1 = lambda t: 4 * exp(10 * (t - T_START)) - 1
        f2 = lambda t: (-1 - exp(-(t - T_END) * 2))
        m = (f2(T_END) - f1(T_START)) / (T_END - T_START)
        p = f1(T_START) - m * T_START
        impulse_graph1 = FunctionGraph(f1, x_range=[0, T_START], color=RED)
        impulse_graph2 = FunctionGraph(
            lambda t: m * t + p, x_range=[T_START, T_END], color=RED
        )
        impulse_graph3 = FunctionGraph(f2, x_range=[T_END, 10], color=RED)

        threshold = DashedLine(-Y, 6 * X - Y, dash_length=0.5)
        resting_potential = Line(4.5 * X - 1.5 * Y, 6 * X - 1.5 * Y)

        mt1 = MathTex("-50 \\text{ mV}", font_size=28).shift(-6.4 * X - 1.5 * Y)
        mt2 = MathTex("+30 \\text{ mV}", font_size=28).shift(-6.4 * X + 2.5 * Y)
        mt3 = MathTex("0 \\text{ mV}", font_size=28).shift(-6.2 * X + Y)

        t1 = Text("Threshold Level", font_size=28).next_to(threshold, UP)
        t2 = Text("Resting potential", font_size=28).next_to(resting_potential, DOWN)
        t3 = Text("Time (ms)", font_size=28).next_to(axes, DOWN)

        impulse = VGroup(impulse_graph1, impulse_graph2, impulse_graph3).shift(
            -5.7 * X - 0.5 * Y
        )
        self.add(axes, mt1, mt2, mt3, threshold, resting_potential, t1, t2, t3)
        self.play(Create(impulse), run_time=3)


class RandomSignal(Scene):
    def construct(self):
        nn = neural_network(0.5)
        dy = frame_height * 0.6 * 0.5
        signal = VGroup(
            *[
                SVGMobject(f"figures/signal{i + 1}.svg").shift((1 - i) * dy * Y)
                for i in range(3)
            ]
        )
        dy_nn = frame_height * 0.5 * 0.5
        arrows = VGroup(*[Arrow(start= -2 * X + (1 - i) * dy * Y, end=0.35 * X + (1 - i) * dy_nn * Y) for i in range(3)])
        objs = VGroup(signal.shift(-4 * X), nn.shift(2 * X), arrows)
        self.add(objs.shift(-objs.get_center()))
        self.add(signal)
