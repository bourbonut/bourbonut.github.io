from glm import vec3
from math import radians, pi, tan 

from .fast import multilines


def addendum_length(m, alpha, ka):
    return pi * m * 0.5 - 2 * tan(alpha) * m * ka


def dedendum_length(m, alpha, kf):
    return pi * m * 0.5 - 2 * tan(alpha) * m * kf


def profile(m, alpha=radians(20), ka=1, kf=1.25):
    # Parameters
    ha = m * ka  #                          addendum height
    hf = m * kf  #                          dedendum height
    p = pi * m  #                           step
    la = addendum_length(m, alpha, ka)  #   addendum length of tooth
    lf = dedendum_length(m, alpha, kf)  #   dedendum length of tooth
    tan_a = tan(alpha)

    #    B ____ C
    #     /    \
    #    /      \
    # A /      D \____ E

    return multilines(
        [
            vec3(-tan_a * hf, -hf, 0),  #                  A
            vec3(tan_a * ha, ha, 0),  #                    B
            vec3(tan_a * ha + la, ha, 0),  #               C
            vec3(p * 0.5 + tan_a * hf, -hf, 0),  #         D
            vec3(p * 0.5 + lf + tan_a * hf, -hf, 0),  #    E
        ],
    )
