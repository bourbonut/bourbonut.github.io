from manim import *
from math import cos, sin, pi
from itertools import starmap
from functools import reduce
from operator import add, itemgetter
import random, pickle, glm
from utils.maths import *

frame_width = config["frame_width"]
frame_height = config["frame_height"]
colors = [PURPLE, PURE_RED, YELLOW, PINK, ORANGE, PURE_BLUE, PURE_GREEN, DARK_BROWN]
times = [0.0972, 0.7141, 0.8308, 0.9018, 0.9234, 0.9387, 0.9473, 0.9527, 0.957, 0.9621]
starting_value = times[0]
times = times[1:] + [0.9682]

def neural_network(radius=0.1, color=BLUE, idx=None):
    idx = random.randint(0, 1234567) if idx is None else idx
    random.seed(42 + idx)
    r = radius
    fcolor = (lambda: random.choice(colors)) if color=="random" else (lambda: color)
    circle_ref = lambda: Circle(r, color=fcolor(), fill_opacity=0.5)
    dy = frame_height * r * 0.5
    dx = dy * 1.5
    dirs = [O, dy * Y, -dy * Y, dx * X + 0.5 * dy * Y, dx * X - 0.5 * dy * Y]
    circles = [circle_ref().shift(d) for d in dirs]
    connections = [(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)]
    lines = [Line(dirs[a] + r * X, dirs[b] - r * X) for a, b in connections]
    barycenter = reduce(add, dirs) / len(dirs)
    return VGroup(*circles, *lines).shift(-barycenter)

class ServerAndCellphones(Scene):
    def construct(self):
        cell = SVGMobject("figures/cellphone.svg")
        cells = [
            cell.copy().shift(2.75 * rotation(2 * pi / 5 * i) * Y) for i in range(5)
        ]
        server = SVGMobject("figures/server.svg", height=1.5)
        self.add(server, *cells)


class ServerAsNN(Scene):
    def construct(self):
        nn = neural_network(0.2, color=BLUE)
        circle = Circle(1.25, color="#ffd36a")
        circle_nn = VGroup(nn, circle)
        server = SVGMobject("figures/server.svg", height=1.5)
        arrow = Arrow(start=LEFT, end=RIGHT)
        self.add(server.next_to(arrow, LEFT), arrow, circle_nn.next_to(arrow, RIGHT))


class Selection(Scene):
    def construct(self):
        nn = neural_network(0.2, color=BLUE)
        nns = [nn.copy() for _ in range(5)]
        circle = Circle(1.25, color="#ffd36a")
        cell = SVGMobject("figures/cellphone.svg")
        cells = [
            cell.copy().shift(2.75 * rotation(2 * pi / 5 * i) * Y) for i in range(5)
        ]
        self.add(nn, circle, *cells)
        self.add(*[Cross().move_to(cells[i].get_center()) for i in (1, 4)])


class UploadToCellphones(Scene):
    def construct(self):
        nn = neural_network(0.2, color=BLUE)
        nns = [nn.copy() for _ in range(5)]
        circle = Circle(1.25, color="#ffd36a")
        cell = SVGMobject("figures/cellphone.svg")
        cells = [
            cell.copy().shift(2.75 * rotation(2 * pi / 5 * i) * Y) for i in range(5)
        ]
        self.add(nn, circle, *[cells[i] for i in (0, 2, 3)])
        self.play(
            *[
                nns[i].animate.scale(0.25).move_to(2.75 * rotation(2 * pi / 5 * i) * Y)
                for i in (0, 2, 3)
            ]
        )
        self.wait(0.3)
        self.wait(0.5)


class Training(Scene):
    def construct(self):
        cell = SVGMobject("figures/cellphone.svg")
        indices = (0, 2, 3)
        cells = [
            cell.copy().shift(2.75 * rotation(2 * pi / 5 * i) * Y) for i in indices
        ]
        dy_up = cell.height * 0.2
        dy_down = cell.height * 0.125
        nn = neural_network(0.2, color=BLUE)
        nns = [
            nn.copy().scale(0.25).move_to(2.75 * rotation(2 * pi / 5 * i) * Y)
            for i in indices
        ]
        rnns = [
            neural_network(0.2, color="random", idx=idx)
            .scale(0.25)
            .move_to(-dy_down * Y + 2.75 * rotation(2 * pi / 5 * i) * Y)
            for idx, i in enumerate(indices)
        ]

        circle = Circle(1.25, color="#ffd36a")
        self.add(nn, circle, *cells, *nns)
        self.play(*[x.animate.shift(-dy_down * Y) for x in nns])
        waiting = SVGMobject("figures/waiting.svg", height=0.5)
        waitings = [
            waiting.copy().shift(dy_up * Y + 2.75 * rotation(2 * pi / 5 * i) * Y)
            for i in indices
        ]
        self.play(*map(FadeIn, waitings))
        self.play(*[Rotate(x, angle=-2 * PI, rate_func=linear) for x in waitings])
        self.play(*[Rotate(x, angle=-2 * PI, rate_func=linear) for x in waitings])
        self.play(*[Rotate(x, angle=-2 * PI, rate_func=linear) for x in waitings])
        self.play(*[Transform(x, y) for x, y in zip(nns, rnns)])


class DownloadToServer(Scene):
    def construct(self):
        cell = SVGMobject("figures/cellphone.svg")
        indices = (0, 2, 3)
        cells = [
            cell.copy().shift(2.75 * rotation(2 * pi / 5 * i) * Y) for i in indices
        ]
        dy_up = cell.height * 0.2
        dy_down = cell.height * 0.125
        rnns = [
            neural_network(0.2, color="random", idx=idx)
            .scale(0.25)
            .move_to(-dy_down * Y + 2.75 * rotation(2 * pi / 5 * i) * Y)
            for idx, i in enumerate(indices)
        ]

        circle = Circle(1.25, color="#ffd36a")
        nn = neural_network(0.2, color=BLUE)
        waiting = SVGMobject("figures/waiting.svg", height=0.5)
        waitings = [
            waiting.copy().shift(dy_up * Y + 2.75 * rotation(2 * pi / 5 * i) * Y)
            for i in indices
        ]
        nn = neural_network(0.2, color=BLUE)
        dirs = [2.75 * Y - 0.25 * frame_width * X, 2.75 * Y + 0.25 * frame_width * X]
        rect = Rectangle(
            width=0.9 * frame_width,
            height=0.5 * 0.7 * frame_height,
        )
        text = Text(
            "Server",
        ).next_to(rect, DOWN)

        self.add(nn, circle, *cells, *rnns, *waitings)
        self.play(
            FadeOut(circle),
            FadeOut(nn),
            *map(FadeOut, waitings),
            *[x.animate.shift(dy_down * Y) for x in rnns],
        )
        lines = [
            Line(2.75 * rotation(2 * pi / 5 * i) * Y, dirs[idx])
            for idx, i in enumerate((2, 3))
        ]
        cells = [cells[0]] + [
            x.move_to(line.get_start()) for x, line in zip(cells[1:], lines)
        ]
        self.play(
            *[MoveAlongPath(cell, line) for cell, line in zip(cells[1:], lines)],
            *[x.animate.move_to(d) for d, x in zip(dirs, rnns[1:])],
        )
        self.play(Create(rect), Write(text))
        self.wait(1)
        self.play(*[rnn.animate.scale(4).shift(0.12 * X - 2.75 * Y) for rnn in rnns])
        self.play(*map(Uncreate, cells))
        self.wait(0.3)


class Aggregation(Scene):
    def construct(self):
        dirs = [-0.25 * frame_width * X, O, 0.25 * frame_width * X]
        a = (1, 0, 2)
        rnns = [
            neural_network(0.2, color="random", idx=j).shift(dirs[i]) for i, j in zip(range(3), a)
        ]
        rect = Rectangle(
            width=0.9 * frame_width,
            height=0.5 * 0.7 * frame_height,
        )
        text = Text(
            "Server",
        ).next_to(rect, DOWN)
        urnn = neural_network(0.2, color="random", idx=-1)
        self.add(text, rect, *rnns)
        self.play(*[rnns[i].animate.shift(-dirs[i]) for i in (0, 2)])
        self.play(*map(FadeOut, rnns), FadeIn(urnn))
        self.wait(1)


class Restart(Scene):
    def construct(self):
        rect = Rectangle(
            width=0.9 * frame_width,
            height=0.5 * 0.7 * frame_height,
        )
        text = Text(
            "Server",
        ).next_to(rect, DOWN)
        urnn = neural_network(0.2, color="random", idx=-1)
        circle = Circle(1.25, color="#ffd36a")
        cell = SVGMobject("figures/cellphone.svg")
        cells = [
            cell.copy().shift(2.75 * rotation(2 * pi / 5 * i) * Y) for i in range(5)
        ]
        nns = [urnn.copy() for _ in range(5)]
        self.add(urnn, text, rect)
        self.play(Unwrite(text), Transform(rect, circle))
        self.play(*map(Create, cells))
        self.play(
            *[
                nns[i].animate.scale(0.25).move_to(2.75 * rotation(2 * pi / 5 * i) * Y)
                for i in range(5)
            ]
        )
        self.wait(1)


class CompleteTraining(Scene):

    CELLINDICES = list(range(5))
    DX = -0.25 * frame_width * X

    def selection(self, cells):
        self.indices = set(random.sample(self.CELLINDICES, k=3))
        others = set(self.CELLINDICES) - self.indices
        self.play(*map(Uncreate, itemgetter(*others)(cells)), run_time=0.3)

    def upload(self):
        self.nns = [self.nn.copy() for _ in self.indices]
        pos = lambda i: self.DX - self.dy_down * Y + 2.75 * rotation(2 * pi / 5 * i) * Y
        movement = lambda i, x: x.animate.scale(0.25).move_to(pos(i))
        self.play(*starmap(movement, zip(self.indices, self.nns)), run_time=0.3)

    def train(self, waiting):
        pos = lambda i: self.DX + self.dy_up * Y + 2.75 * rotation(2 * pi / 5 * i) * Y
        waitings = [waiting.copy().shift(pos(i)) for i in self.indices]
        new = (
            lambda nn: neural_network(0.2, "random")
            .scale(0.25)
            .move_to(nn.get_center())
        )
        new_nns = [new(nn) for nn in self.nns]
        self.play(*map(FadeIn, waitings), run_time=0.3)
        self.play(
            *map(lambda x: Rotate(x, angle=-2 * PI, rate_func=linear), waitings),
            run_time=0.3,
        )
        self.play(*starmap(Transform, zip(self.nns, new_nns)), run_time=0.3)
        self.play(*map(FadeOut, waitings), run_time=0.3)

    def aggregate(self):
        new_nn = neural_network(0.2, "random").move_to(self.nn.get_center())
        movement = lambda x: x.animate.scale(4).move_to(self.nn.get_center())
        self.play(*map(movement, self.nns), run_time=0.3)
        self.play(*map(FadeOut, self.nns), Transform(self.nn, new_nn), run_time=0.3)

    def reset(self, cells, ref):
        pos = lambda i: self.DX + 2.75 * rotation(2 * pi / 5 * i) * Y
        others = set(self.CELLINDICES) - self.indices
        for i in others:
            cells[i] = ref.copy().shift(pos(i))
        self.play(*map(Create, itemgetter(*others)(cells)), run_time=0.3)

    def point(self, index, y):
        step = 10
        length = 6  # x_length, y_length
        r = length / step
        x = index * r
        y *= r / 0.1
        return self.origin + x * X + y * Y

    def construct(self):
        plane = NumberPlane(
            x_range=(0, 10, 1),
            y_range=[0, 1, 0.1],
            x_length=6,
            y_length=6,
            axis_config={"include_numbers": True},
            y_axis_config={"label_direction": LEFT},
            tips=False,
        )
        plane = plane.shift(0.25 * frame_width * X)
        x_label = Text("Rounds").scale(0.35).next_to(plane, DOWN)
        y_label = Text("Accuracy").scale(0.35).rotate(PI / 2).next_to(plane, LEFT)
        self.origin = plane.get_origin()
        self.add(plane, x_label, y_label)

        self.nn = neural_network(0.2, color=BLUE).shift(self.DX)
        circle = Circle(1.25, color="#ffd36a").shift(self.DX)
        cell = SVGMobject("figures/cellphone.svg")
        waiting = SVGMobject("figures/waiting.svg", height=0.5)
        self.dy_up = cell.height * 0.2
        self.dy_down = cell.height * 0.125
        pos = lambda i: self.DX + 2.75 * rotation(2 * pi / 5 * i) * Y
        cells = [cell.copy().shift(pos(i)) for i in range(5)]
        start = self.point(0, starting_value)
        dot = Dot(start, 0.08)
        self.add(circle, self.nn, *cells)
        for i, time in enumerate(times):
            self.selection(cells)
            self.upload()
            self.train(waiting)
            self.aggregate()
            next_ = self.point(i + 1, time)
            line = Line(start, next_)
            vg = VGroup(dot, line)
            self.play(Create(vg))
            dot = Dot(next_, 0.08)
            start = next_
            self.reset(cells, cell)
