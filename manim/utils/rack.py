from glm import vec3
from math import radians, pi, tan

from .fast import multilines


def addendum_length(m, alpha, ka):
    return pi * m * 0.5 - 2 * tan(alpha) * m * ka


def dedendum_length(m, alpha, kf):
    return pi * m * 0.5 - 2 * tan(alpha) * m * kf

def up2down_length(m, alpha, ka, kf):
    return tan(alpha) * m * (ka + kf)


def profile(m, alpha=radians(20), ka=1, kf=1.25):
    # Parameters
    ha = m * ka  #                              addendum height
    hf = m * kf  #                              dedendum height
    p = pi * m  #                               step
    la = addendum_length(m, alpha, ka)  #       addendum length of tooth
    lf = dedendum_length(m, alpha, kf)  #       dedendum length of tooth
    lud = up2down_length(m, alpha, ka, kf) #    length up to down
    tan_a = tan(alpha)

    #    B ____ C
    #     /    \
    #    /      \
    # A /      D \____ E

    return multilines(
        [
            vec3(0, -hf, 0),  #                    A
            vec3(lud, ha, 0),  #                   B
            vec3(lud + la, ha, 0),  #              C
            vec3(2 * lud + la, -hf, 0),  #         D
            vec3(2 * lud + la + lf, -hf, 0),  #    E
        ],
    )
