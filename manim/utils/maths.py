from glm import vec3, mat3, dot, length
from math import cos, sin, acos

O = vec3(0, 0, 0)
X = vec3(1, 0, 0)
Y = vec3(0, 1, 0)
Z = vec3(0, 0, 1)


def anglebt(a, b):
    return acos(dot(a, b) / (length(a) * length(b)))


def u(t):
    return vec3(cos(t), sin(t), 0)


def v(t):
    return vec3(-sin(t), cos(t), 0)


def rotation(t):
    return mat3(u(t), v(t), Z)
