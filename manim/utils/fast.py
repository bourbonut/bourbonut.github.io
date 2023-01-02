from manim import VGroup, Line, Dot
from itertools import pairwise, starmap, chain

from .maths import X


def multilines(points, closed=False):
    return VGroup(
        *chain(
            (starmap(Line, pairwise(points))),
            (map(lambda p: Dot(p, radius=0.02), points)),
            (Line(points[0], points[-1])) if closed else (),
        )
    )


def repeat(vgroup, times, step, direction=X):
    final_vgroup = vgroup.copy()
    for i in range(times):
        translated_vgroup = vgroup.copy().shift(direction * i * step)
        for line in translated_vgroup:
            final_vgroup += line
    return final_vgroup


def revolution(vgroup, angle, times):
    final_vgroup = vgroup.copy()
    for i in range(times):
        rotated_vgroup = vgroup.copy().rotate_about_origin(angle * i)
        for line in rotated_vgroup:
            final_vgroup += line
    return final_vgroup
