from typing import Callable
from manim import *
from glm import *
from utils import *
from utils.maths import O, X, Y, Z, rotation, anglebt
from math import pi, atan2, sin, cos, asin, acos, tan


def get_pitch_cone_angle(
    z_pinion: int, z_wheel: int, shaft_angle: float = 0.5 * pi
) -> float:
    """
    Return the pitch cone angle of the pinion called `gamma_p`.
    The pitch cone angle of the wheel is equal to `shaft_angle - gamma_p`.

    Parameters
    ----------
    z_pinion : int
        The number of teeth on the bevel pinion
    z_wheel : int
        The number of teeth on the bevel wheel
    shaft_angle : float
        The shaft angle

    Returns
    -------
    float
        The pitch cone angle
    """
    return atan2(sin(shaft_angle), ((z_wheel / z_pinion) + cos(shaft_angle)))


def spherical_involute(cone_angle:float, t0:float, t:float) -> vec3:
    """
    Return spherical involute function

    Parameters
    ----------
    t : float
        The angular position
    t0 : float
        The difference phase
    cone_angle : float
        The cone angle

    Returns
    -------
    vec3
        A point on the spherical involute
    """
    cos_g, sin_g = cos(cone_angle), sin(cone_angle)
    return vec3(
        sin_g * cos(t * sin_g) * cos(t + t0) + sin(t * sin_g) * sin(t + t0),
        sin_g * cos(t * sin_g) * sin(t + t0) - sin(t * sin_g) * cos(t + t0),
        cos_g * cos(t * sin_g),
    )

class BevelGearSection(Scene):
    def construct(self):
        z = z_pinion = 10
        z_wheel = 20
        m = 0.5
        step = pi * m
        ka = 1
        kd = 1.25
        pressure_angle = pi / 9
        pitch_cone_angle = get_pitch_cone_angle(z_pinion, z_wheel)

        gamma_p = pitch_cone_angle  # for convenience
        gamma_b = asin(cos(pressure_angle) * sin(gamma_p))
        cos_b, sin_b = cos(gamma_b), sin(gamma_b)
        rp = z * step / (2 * pi)
        rho1 = rp / sin(gamma_p)
        rho0 = 2 * rho1 / 3
        k = sin(gamma_p) / z
        gamma_r = gamma_p - atan2(2 * kd * k, 1)
        gamma_f = gamma_p + atan2(2 * ka * k, 1)

        # Make teeth (up and down)
        perp = rotation(gamma_p) * Y

        points = [
            rho0 * rotation(gamma_p) * X - 2 / 3 * kd * m * perp,
            rho0 * rotation(gamma_p) * X + 2 / 3 * ka * m * perp,
            rho1 * rotation(gamma_p) * X + ka * m * perp,
            rho1 * rotation(gamma_p) * X - kd * m * perp,
        ]
        up_tooth = Polygon(*points, color=WHITE)

        # Make bevel gear body
        rev_points = [vec3(p.x, -p.y, 0) for p in points]
        down_tooth = Polygon(*rev_points, color=WHITE)

        body_points = [
            rho0 * rotation(gamma_p) * X - 2 * kd * m * perp,
            rho1 * rotation(gamma_p) * X - 2.2 * kd * m * perp,
        ]
        body_points += [vec3(p.x, -p.y, 0) for p in reversed(body_points)]

        body = Polygon(
            body_points[0],
            points[1],
            points[2],
            body_points[1],
            body_points[2],
            rev_points[2],
            rev_points[1],
            body_points[3],
            color=WHITE,
        )

        # Make dashed lines
        p1 = rho0 * rotation(gamma_p) * X - 2 / 3 * kd * m * perp
        p2 = rho0 * rotation(gamma_p) * X + 2 / 3 * ka * m * perp
        dashed_lines = [
            DashedLine(O, p1, color=BLUE),
            DashedLine(O, p2, color=BLUE),
            DashedLine(O, rho1 * rotation(gamma_p) * X, color=YELLOW),
            DashedLine(O, rho1 * rotation(gamma_b) * X, color=RED),
            DashedLine(O, vec3(p1.x, -p1.y, 0), color=BLUE),
            DashedLine(O, vec3(p2.x, -p2.y, 0), color=BLUE),
            DashedLine(O, rho1 * rotation(-gamma_p) * X, color=YELLOW),
            DashedLine(O, rho1 * rotation(-gamma_b) * X, color=RED),
            DashedLine(O, rho1 * X),
        ]

        # Rho representation
        rho_repr = VGroup(
            Line(perp * 0.2, perp * 0.9),
            Line(
                rho1 * rotation(gamma_p) * X + perp * 0.6,
                rho1 * rotation(gamma_p) * X + perp * 1.3,
            ),
            DoubleArrow(
                perp * 0.55,
                rho1 * rotation(gamma_p) * X + perp * 0.95,
                stroke_width=3,
                max_tip_length_to_length_ratio=0.04,
                buff=0.025,
            ),
            MathTex("\\rho_1").move_to(
                0.5 * rho1 * rotation(gamma_p) * X + perp
            ),
        )

        # Gamma representation
        gamma_repr = []
        kr = 0.4
        for tex, gamma, k, color in [
            ("\\gamma_r", gamma_r, 0.34, BLUE),
            ("\\gamma_b", gamma_b, 0.46, RED),
            ("\\gamma_p", gamma_p, 0.58, YELLOW),
            ("\\gamma_f", gamma_f, 0.7, BLUE),
        ]:
            gamma_repr.append(
                VGroup(
                    ArcBetweenPoints(
                        rho1 * k * X,
                        rho1 * k * rotation(gamma) * X,
                        radius=rho1 * k,
                        color=color,
                    )
                    .add_tip(tip_width=0.15, tip_length=0.15, at_start=True)
                    .add_tip(tip_width=0.15, tip_length=0.15),
                    MathTex(tex, color=color).move_to(
                        rho1 * (k + 0.05) * rotation(gamma * 0.5) * X
                    ),
                )
            )

        # ka and kd representations
        middle_line = Line(
            rho1 * 1.02 * rotation(gamma_p) * X, rho1 * 1.2 * rotation(gamma_p) * X
        )
        ka_repr = VGroup(
            Line(
                rho1 * 1.02 * rotation(gamma_p) * X + ka * m * perp,
                rho1 * 1.2 * rotation(gamma_p) * X + ka * m * perp,
            ),
            Line(
                rho1 * 1.08 * rotation(gamma_p) * X + ka * m * perp * 1.5,
                rho1 * 1.08 * rotation(gamma_p) * X - ka * m * perp * 0.5,
                color=BLUE,
            ),
            Cross(scale_factor=0.08, stroke_color=BLUE, stroke_width=3)
            .rotate(gamma_p)
            .move_to(rho1 * 1.08 * rotation(gamma_p) * X + ka * m * perp),
            Cross(scale_factor=0.08, stroke_color=BLUE, stroke_width=3)
            .rotate(gamma_p)
            .move_to(rho1 * 1.08 * rotation(gamma_p) * X),
            MathTex("k_a \\cdot m", color=BLUE)
            .move_to(rho1 * 1.15 * rotation(gamma_p) * X + ka * m * perp * 2.2),
        )
        kd_repr = VGroup(
            Line(
                rho1 * 1.02 * rotation(gamma_p) * X - kd * m * perp,
                rho1 * 1.2 * rotation(gamma_p) * X - kd * m * perp,
            ),
            Line(
                rho1 * 1.14 * rotation(gamma_p) * X - kd * m * perp * 1.5,
                rho1 * 1.14 * rotation(gamma_p) * X + kd * m * perp * 0.5,
                color=BLUE,
            ),
            Cross(scale_factor=0.08, stroke_color=BLUE, stroke_width=3)
            .rotate(gamma_p)
            .move_to(rho1 * 1.14 * rotation(gamma_p) * X - kd * m * perp),
            Cross(scale_factor=0.08, stroke_color=BLUE, stroke_width=3)
            .rotate(gamma_p)
            .move_to(rho1 * 1.14 * rotation(gamma_p) * X),
            MathTex("k_d \\cdot m", color=BLUE)
            .move_to(rho1 * 1.13 * rotation(gamma_p) * X - kd * m * perp * 2),
        )

        # rp representation
        rp_repr = VGroup(
            DoubleArrow(
                rp / tan(gamma_p) * X,
                rho1 * rotation(gamma_p) * X,
                stroke_width=3,
                max_tip_length_to_length_ratio=0.04,
                buff=0.025,
                color=YELLOW,
            ),
            MathTex("r_p", color=YELLOW).move_to(0.95 * rp / tan(gamma_p) * X + rp * 0.5 * Y)
        )

        group = VGroup(
            up_tooth,
            down_tooth,
            body,
            *dashed_lines,
            rho_repr,
            *gamma_repr,
            middle_line,
            ka_repr,
            kd_repr,
            rp_repr,
        ).center()
        self.add(group)

def arc(
    func: Callable[[float], vec3],
    n: int = 10,
    color: str | ManimColor = WHITE,
) -> list[Line]:
    """
    Returns an arc

    Parameters
    ----------
    func : Callable[[float], vec3]
        Function which generates the arc
    n : int
        Resolution of the arc
    color : str | ManimColor
        Color of the arc

    Returns
    -------
    list[Line]
        List of lines which composes the arc
    """
    return [
        Line(func(i / n), func((i + 1) / n), color=color) for i in range(n)
    ] + [
        Line(func((i + 0.5) / n), func((i + 1.5) / n), color=color) for i in range(n - 1)
    ]

def arc_with_arrows(
    func: Callable[[float], vec3],
    tex_content: str,
    direction: vec3 = Y,
    n: int = 10,
    height: float = 0.15,
    reverse: bool = False,
    font_size: int = 48,
    color: str | ManimColor = WHITE,
    phi_rot: float = 75 * DEGREES,
    theta_rot: float = pi / 2 - 20 * DEGREES,
) -> list[VMobject]:
    """
    Generates an arc with arrows

    Parameters
    ----------
    func : Callable[[float], vec3]
        Function which generates the arc
    tex_content : str
        Tex content
    direction : vec3
        Dirction for the tex content
    n : int
        Resolution of the arc
    height : float
        Height of arrows
    reverse : bool
        Reverse arrows
    font_size : int
        Font size of text
    color : str | ManimColor
        Color of the arc
    phi_rot : float
        Phi camera angle for tex content orientation
    theta_rot : float
        Theta camera angle for tex content orientation

    Returns
    -------
    list[VMobject]
        List of VMobject which composes the arc
    """
    s1, s2, e1, e2 = (3, 0, n - 3, n)
    arrow_1 = Arrow3D(func(s1 / n), func(s2 / n), resolution=0, thickness=0, height=height, color=color)
    arrow_2 = Arrow3D(func(e1 / n), func(e2 / n), resolution=0, thickness=0, height=height, color=color)
    if reverse:
        arrow_1 = arrow_1.rotate(pi, func(s2 / n), about_point=func(s2 / n))
        arrow_2 = arrow_2.rotate(pi, func(e2 / n), about_point=func(e2 / n))
    return [
        Line(func(i / n), func((i + 1) / n), color=color)
        for i in range(n)
    ] + [
        Line(func((i + 0.5) / n), func((i + 1.5) / n), color=color)
        for i in range(n - 1)
    ] + [
        arrow_1,
        arrow_2,
        MathTex(tex_content, color=color, font_size=font_size)
        .rotate(phi_rot, X, about_point=O)
        .rotate(theta_rot, Z, about_point=O)
        .next_to(func(0.5), direction=direction),
    ]

class SphericalRepr(ThreeDScene):

    def construct(self):
        z = z_pinion = 10
        z_wheel = 20
        m = 0.5
        step = pi * m
        ka = 1
        kd = 1.25
        pressure_angle = pi / 9
        pitch_cone_angle = get_pitch_cone_angle(z_pinion, z_wheel)

        gamma_p = pitch_cone_angle  # for convenience
        gamma_b = asin(cos(pressure_angle) * sin(gamma_p))
        cos_b, sin_b = cos(gamma_b), sin(gamma_b)
        rp = z * step / (2 * pi)
        rho1 = rp / sin(gamma_p)
        rho0 = 2 * rho1 / 3
        k = sin(gamma_p) / z
        gamma_r = gamma_p - atan2(2 * kd * k, 1)
        gamma_f = gamma_p + atan2(2 * ka * k, 1)

        epsilon = pi / 3
        phi = epsilon * sin(gamma_b)

        self.set_camera_orientation(phi=75 * DEGREES, theta=-20 * DEGREES, zoom=0.8, focal_distance=1000)

        rho1_circle = Circle(rho1, color=WHITE)

        rb = rho1 * sin(gamma_b)
        base_circle = (
            Circle(radius=rb, color=WHITE)
            .rotate(pi / 2, Y, about_point=O)
            .move_to(rho1 * cos(gamma_b) * X)
            .rotate(gamma_b, Y, about_point=O)
        )

        positions = {
            "O": rho1 * rotate(gamma_b, Y) * X,
            "P": rho1 * rotate(phi, Z) * X,
            "M": rho1 * X,
            "Q": rho1 * rotate(-epsilon, rotate(gamma_b, Y) * X) * X,
        }

        center = VGroup(
            Dot3D(positions["O"]),
            MathTex("O")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["O"], direction=-normalize(Z + Y)),
        )

        P = VGroup(
            Dot3D(positions["P"]),
            MathTex("P")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["P"], direction=normalize(Z + Y)),
        )

        M = VGroup(
            Dot3D(positions["M"]),
            MathTex("M")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["M"], direction=-normalize(Z + Y)),
        )

        Q = VGroup(
            Dot3D(positions["Q"]),
            MathTex("Q")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["Q"], direction=Y),
        )

        points = [center, P, M, Q]

        cone = Cone(rb, rho1 * cos(gamma_b), -rotate(gamma_b, Y) * X, fill_color=DARK_BLUE).set_opacity(0.05)
        sphere = Sphere(radius=rho1, u_range=(0, pi), fill_color=DARK_GRAY).rotate(-pi / 2, X, about_point=O).set_opacity(0.01)

        shapes = [sphere, cone]

        lines = [
            DashedLine(O, positions["O"], dash_length=0.2, color=GREEN),
            # DashedLine(O, rotate(-epsilon, positions["O"]) * rotate(pi / 2, Z) * positions["M"], dash_length=0.2, color=BLUE),
        ]

        n = 10
        # Manim is BUGGED ! If you set normal_vector, it does NOT work. I created my own arcs
        # Arc between O and M
        lines += arc(
            func=lambda t: rotate(gamma_b * t, cross(positions["O"], positions["M"])) * positions["O"],
            n=n,
            color=GREEN,
        )

        # Arc between O and Q
        lines += arc(
            func=lambda t: rotate(gamma_b * t, cross(positions["O"], positions["Q"])) * positions["O"],
            n=n,
            color=BLUE,
        )

        # epsilon (OM, OQ)
        arrows = arc_with_arrows(
            func=lambda t: rotate(-epsilon * t, rotate(gamma_b, Y) * X) * rho1 * rotate(gamma_b * 0.5, Y) * X,
            tex_content="\\epsilon",
            direction=0.5 * normalize(Z + Y),
            n=n,
            height=0.1,
            color=BLUE,
        )

        # gamma_b (O_S O, O_S M) where O_S is the origin
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rho1 * rotate(gamma_b * t, Y) * X,
            tex_content="\\gamma_b",
            direction=Y,
            n=n,
            height=0.1,
            color=GREEN,
        )

        # phi (O_S M, O_S P) where O_S is the origin
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rho1 * rotate(phi * t, Z) * X,
            tex_content="\\varphi",
            direction=0.5 * normalize(X - Z),
            n=n,
            height=0.1,
            color=RED,
        )

        # epsilon (y2, y3)
        arrows += arc_with_arrows(
            func=lambda t: rotate(-epsilon * t, rotate(gamma_b, Y) * X) * 0.7 * 0.5 * rho1 * rotate(gamma_b, Y) * Z,
            tex_content="\\epsilon",
            n=n,
            direction=0.5 * normalize(Y + Z),
            height=0.15,
            color=BLUE,
        )

        # gamma_b (y0, y2)
        arrows += arc_with_arrows(
            func=lambda t: 0.7 * 0.5 * rho1 * rotate(gamma_b * t, Y) * Z,
            tex_content="\\gamma_b",
            direction=-Y,
            n=n,
            height=0.05,
            color=GREEN,
        )

        # phi (x0, x1)
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rotate(phi * t, Z) * rho1 * Y,
            tex_content="\\varphi",
            direction=normalize(-X + 1.5 * Y),
            n=n,
            height=0.1,
            color=RED,
        )

        rho1_repr = VGroup(
            Arrow3D(O, rho1 * rotate(-pi / 4, Z) * X),
            MathTex("\\rho_1")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(rho1 * 0.5 * rotate(-pi / 4, Z) * X, direction=Z),
        )

        rb_repr = VGroup(
            Arrow3D(cos(gamma_b) * positions["O"], rho1 * rotate(epsilon, positions["O"]) * X),
            MathTex("r_b")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(
                rho1 * cos(gamma_b) / cos(gamma_b * 0.5) * rotate(epsilon, positions["O"]) * rotate(0.5 * gamma_b, Y) * X,
                direction=-Z
            ),
        )

        base_text = (
            Text("Base circle", font_size=32)
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(base_circle, direction=-Y)
        )

        # # Add curve
        rot = (
            rotate(pi / 2, positions["Q"]) *
            rotate(-epsilon, positions["O"]) *
            rotate(-pi / 2 - phi, Z) *
            rotate(-pi / 2, X) * rotate(pi, Z)
        )
        curve = ParametricFunction(
            lambda t: rho1 * rot * spherical_involute(gamma_b, 0, t),
            t_range=[0, epsilon],
            color=YELLOW
        )

        system0 = VGroup(
            Arrow3D(O, 1.1 * positions["P"]),
            Arrow3D(O, 1.1 * rotate(pi / 2, Z) * positions["P"]),
            Arrow3D(O, 0.5 * rho1 * Z),
            MathTex("\\overrightarrow{x_0}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(1.1 * rotate(pi / 2, Z) * positions["P"], direction=rotate(pi / 2, Z) * positions["P"] / rho1),
            MathTex("\\overrightarrow{z_0}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(1.1 * positions["P"], direction=Y),
            MathTex("\\overrightarrow{y_0}, \\overrightarrow{y_1}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(0.5 * rho1 * Z, direction=Y),
        )

        system1 = VGroup(
            Arrow3D(O, 1.1 * positions["M"]),
            Arrow3D(O, 1.1 * rotate(pi / 2, Z) * positions["M"]),
            MathTex("\\overrightarrow{x_1}, \\overrightarrow{x_2}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(1.1 * rotate(pi / 2, Z) * positions["M"], direction=rotate(pi / 2, Z) * positions["M"] / rho1),
            MathTex("\\overrightarrow{z_1}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(1.1 * positions["M"], direction=positions["M"] / rho1),
        ).set_color(RED)

        system2 = VGroup(
            Arrow3D(positions["O"], 1.1 * positions["O"]),
            Arrow3D(O, 0.5 * rotate(-pi / 2, Y) * positions["O"]),
            MathTex("\\overrightarrow{y_2}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(0.5 * rotate(-pi / 2, Y) * positions["O"], direction=Y),
            MathTex("\\overrightarrow{z_2}, \\overrightarrow{z_3}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(1.1 * positions["O"], direction=positions["O"] / rho1),
        ).set_color(GREEN)

        system3 = VGroup(
            Arrow3D(O, 0.5 * rotate(-epsilon, positions["O"]) * rotate(-pi / 2, Y) * positions["O"]),
            # Arrow3D(
            #     rotate(-epsilon, positions["O"]) * rotate(pi / 2, Z) * positions["M"],
            #     1.1 * rotate(-epsilon, positions["O"]) * rotate(pi / 2, Z) * positions["M"]
            # ),
            MathTex("\\overrightarrow{y_3}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(0.6 * rotate(-epsilon, positions["O"]) * rotate(-pi / 2, Y) * rho1 * cos(gamma_b) * rotate(gamma_b, Y) * X),
            # MathTex("\\overrightarrow{x_3}")
            # .rotate(75 * DEGREES, X, about_point=O)
            # .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            # .move_to(1.2 * rotate(-epsilon, positions["O"]) * rotate(pi / 2, Z) * rho1 * X),
        ).set_color(BLUE)

        group = VGroup(
            *shapes,
            *points,
            *lines,
            *arrows,
            rho1_circle,
            base_circle,
            rho1_repr,
            rb_repr,
            base_text,
            curve,
            system0,
            system1,
            system2,
            system3,
        )

        self.add(group.move_to(0.2 * Z + Y))

class SphericalRepr2(ThreeDScene):

    def construct(self):
        z = z_pinion = 10
        z_wheel = 20
        m = 0.5
        step = pi * m
        ka = 1
        kd = 1.25
        pressure_angle = pi / 9
        pitch_cone_angle = get_pitch_cone_angle(z_pinion, z_wheel)

        gamma_p = pitch_cone_angle  # for convenience
        gamma_b = asin(cos(pressure_angle) * sin(gamma_p))
        cos_b, sin_b = cos(gamma_b), sin(gamma_b)
        rp = z * step / (2 * pi)
        rho1 = rp / sin(gamma_p)
        rho0 = 2 * rho1 / 3
        k = sin(gamma_p) / z
        gamma_r = gamma_p - atan2(2 * kd * k, 1)
        gamma_f = gamma_p + atan2(2 * ka * k, 1)

        epsilon = pi / 3
        phi = epsilon * sin(gamma_b)

        self.set_camera_orientation(phi=75 * DEGREES, theta=-20 * DEGREES, zoom=0.8, focal_distance=1000)

        rho1_circle = Circle(rho1, color=WHITE)

        rb = rho1 * sin(gamma_b)
        base_circle = (
            Circle(radius=rb, color=WHITE)
            .rotate(pi / 2, Y, about_point=O)
            .move_to(rho1 * cos(gamma_b) * X)
            .rotate(gamma_b, Y, about_point=O)
        )

        positions = {
            "O": rho1 * rotate(gamma_b, Y) * X,
            "P": rho1 * rotate(phi, Z) * X,
            "M": rho1 * X,
            "Q": rho1 * rotate(-epsilon, rotate(gamma_b, Y) * X) * X,
        }

        center = VGroup(
            Dot3D(positions["O"]),
            MathTex("O")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["O"], direction=-normalize(Z + Y)),
        )

        P = VGroup(
            Dot3D(positions["P"]),
            MathTex("P")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["P"], direction=normalize(Z + Y)),
        )

        M = VGroup(
            Dot3D(positions["M"]),
            MathTex("M")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["M"], direction=-normalize(Z + Y)),
        )

        Q = VGroup(
            Dot3D(positions["Q"]),
            MathTex("Q")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(positions["Q"], direction=Y),
        )

        points = [center, P, M, Q]

        cone = Cone(rb, rho1 * cos(gamma_b), -rotate(gamma_b, Y) * X, fill_color=DARK_BLUE).set_opacity(0.05)
        sphere = Sphere(radius=rho1, u_range=(0, pi), fill_color=DARK_GRAY).rotate(-pi / 2, X, about_point=O).set_opacity(0.01)

        shapes = [sphere, cone]

        lines = [
            DashedLine(O, positions["O"], dash_length=0.2, color=GREEN),
            Line3D(O, positions["P"], color=WHITE),
            Line3D(O, positions["M"], color=RED),
        ]

        n = 10
        # Manim is BUGGED ! If you set normal_vector, it does NOT work. I created my own arcs
        # Arc between O and M
        lines += arc(
            func=lambda t: rotate(gamma_b * t, cross(positions["O"], positions["M"])) * positions["O"],
            n=n,
            color=GREEN,
        )

        # Arc between O and Q
        lines += arc(
            func=lambda t: rotate(gamma_b * t, cross(positions["O"], positions["Q"])) * positions["O"],
            n=n,
            color=BLUE,
        )

        # Missing angles
        gamma = anglebt(positions["O"], positions["P"])
        phi2 = anglebt(positions["M"] - positions["O"], positions["P"] - positions["O"])
        theta = anglebt(positions["P"] - positions["O"], positions["Q"] - positions["O"])
        eta = anglebt(positions["M"] - positions["P"], positions["O"] - positions["P"])

        # Arc between O and P
        lines += arc(
            func=lambda t: rotate(gamma * t, cross(positions["O"], positions["P"])) * positions["O"],
            n=n,
            color=WHITE,
        )

        # epsilon (OM, OQ)
        arrows = arc_with_arrows(
            func=lambda t: rotate(-epsilon * t, positions["O"]) * rho1 * rotate(gamma_b * 0.5, Y) * X,
            tex_content="\\epsilon",
            direction=0.3 * normalize(Z + Y),
            n=n,
            height=0.1,
            color=BLUE,
        )

        # gamma_b (O_S O, O_S M) where O_S is the origin
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rho1 * rotate(gamma_b * t, Y) * X,
            tex_content="\\gamma_b",
            direction=-0.5 * Y,
            n=n,
            height=0.1,
            color=GREEN,
        )

        # phi (O_S M, O_S P) where O_S is the origin
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rho1 * rotate(phi * t, Z) * X,
            tex_content="\\varphi",
            direction=0.3 * normalize(X - Z),
            n=n,
            height=0.1,
            color=RED,
        )

        # phi2 (OM, OP)
        arrows += arc_with_arrows(
            func=lambda t: rho1 * rotate(-phi2 * t, positions["O"]) * rotate(gamma_b * 0.2, Y) * X,
            tex_content="\\phi",
            direction=0.3 * normalize(Z + Y),
            n=n,
            height=0.1,
            color=WHITE,
        )

        # theta (OP, OQ)
        arrows += arc_with_arrows(
            func=lambda t: rho1 * rotate(-theta * t - phi2, positions["O"]) * rotate(gamma_b * 0.2, Y) * X,
            tex_content="\\theta",
            direction=0.3 * normalize(0.5 * Z + Y),
            n=n,
            height=0.05,
            color=WHITE,
        )

        # eta (PM, PO)
        arrows += arc_with_arrows(
            func=lambda t: rho1 * rotate(eta * t, positions["P"]) * rotate(phi * 0.7, Z) * X,
            tex_content="\\eta",
            direction=-0.1 * rotate(eta * 0.5, X) * Y,
            n=n,
            height=0.05,
            color=WHITE,
        )

        # gamma (O_S O, O_S P) where O_S is the origin
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rotate(gamma * t, cross(positions["O"], positions["P"])) * positions["O"],
            tex_content="\\gamma",
            direction=Y,
            n=n,
            height=0.1,
            color=WHITE,
        )

        rho1_repr = VGroup(
            Arrow3D(O, rho1 * rotate(-pi / 4, Z) * X),
            MathTex("\\rho_1")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(rho1 * 0.5 * rotate(-pi / 4, Z) * X, direction=Z),
        )

        rb_repr = VGroup(
            Arrow3D(cos(gamma_b) * positions["O"], rho1 * rotate(epsilon, positions["O"]) * X),
            MathTex("r_b")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(
                rho1 * cos(gamma_b) / cos(gamma_b * 0.5) * rotate(epsilon, positions["O"]) * rotate(0.5 * gamma_b, Y) * X,
                direction=-Z
            ),
        )

        base_text = (
            Text("Base circle", font_size=32)
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .next_to(base_circle, direction=-Y)
        )

        # Add curve
        rot = (
            rotate(pi / 2, positions["Q"]) *
            rotate(-epsilon, positions["O"]) *
            rotate(-pi / 2 - phi, Z) *
            rotate(-pi / 2, X) * rotate(pi, Z)
        )
        curve = ParametricFunction(
            lambda t: rho1 * rot * spherical_involute(gamma_b, 0, t),
            t_range=[0, epsilon],
            color=YELLOW
        )

        group = VGroup(
            *shapes,
            *points,
            *lines,
            *arrows,
            rho1_circle,
            base_circle,
            rho1_repr,
            rb_repr,
            base_text,
            curve,
        ).move_to(0.2 * Z)
        self.add(group)


class BevelFlat(ThreeDScene):

    def construct(self):
        z_pinion = 10
        z_wheel = 20
        m = 0.3
        step = pi * m
        ka = 1
        kd = 1.25
        pressure_angle = pi / 9
        pitch_cone_angle = get_pitch_cone_angle(z_pinion, z_wheel)

        # Pinion
        gamma_pp = pitch_cone_angle
        gamma_bp = asin(cos(pressure_angle) * sin(gamma_pp))
        rpp = z_pinion * step / (2 * pi)
        rho1 = rho1p = rpp / sin(gamma_pp)
        rbp = rho1p * sin(gamma_bp)
        phip = acos(tan(gamma_bp) / tan(gamma_pp))

        # Wheel
        gamma_pw = pi * 0.5 - pitch_cone_angle
        gamma_bw = asin(cos(pressure_angle) * sin(gamma_pw))
        rpw = z_wheel * step / (2 * pi)
        rho1w = rpw / sin(gamma_pw)
        rbw = rho1w * sin(gamma_bw)
        phiw = acos(tan(gamma_bw) / tan(gamma_pw))

        self.set_camera_orientation(phi=90 * DEGREES, theta=00 * DEGREES, zoom=1.1, focal_distance=1000)

        O1 = rho1 * rotate(gamma_pp, Y) * X
        O2 = rho1 * rotate(-gamma_pw, Y) * X
        P = rho1 * X
        M1 = rotate(-phip, O1) * rotate(gamma_pp, Y) * (rho1 * cos(gamma_bp) * X + rbp * Z)
        M2 = rotate(-phiw, O2) * rotate(-gamma_pw, Y) * (rho1 * cos(gamma_bw) * X - rbw * Z)

        sphere = Sphere(radius=rho1, fill_color=DARK_GRAY).set_opacity(0.01)

        circles = VGroup(
            Circle(rho1, color=WHITE),
            Circle(rho1, color=WHITE).rotate(pi / 2, Y, about_point=O),
            Circle(rho1, color=WHITE).rotate(-pressure_angle, X, about_point=O),
            Circle(rbw, color=RED)
            .rotate_about_origin(pi / 2 - gamma_pw, Y)
            .move_to(rho1 * cos(gamma_bw) * rotate(-gamma_pw, Y) * X),
            Circle(rbp, color=GREEN)
            .rotate_about_origin(pi / 2 + gamma_pp, Y)
            .move_to(rho1 * cos(gamma_bp) * rotate(gamma_pp, Y) * X)
        )

        dots = VGroup(
            Dot3D(O1, color=GREEN),
            Dot3D(O2, color=RED),
            Dot3D(P, color=WHITE),
            Dot3D(M1, color=GREEN),
            Dot3D(M2, color=RED),
        )

        tex = VGroup(
            MathTex("O_1", color=GREEN, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(O1, direction=-Z),
            MathTex("O_2", color=RED, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(O2, direction=0.3 * Z),
            MathTex("M_1", color=GREEN, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(M1, direction=Y),
            MathTex("M_2", color=RED, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(M2, direction=-Y),
            MathTex("P_0", color=WHITE, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(P, direction=0.3 * normalize(Z + Y)),
            MathTex("\\gamma_{b_1}", color=GREEN, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(
                rotate(-phip, O1) * rotate(gamma_pp, Y) * (rho1 * cos(gamma_bp * 0.5) * X + rbp * 0.5 * Z),
                direction=0.3 * Y,
            ),
            MathTex("\\gamma_{b_2}", color=RED, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(
                rotate(-phiw, O2) * rotate(-gamma_pw, Y) * (rho1 * cos(gamma_bw * 0.5) * X - rbw * 0.5 * Z),
                direction=0.3 * normalize(Z - Y)
            ),
            MathTex("\\gamma_{p_1}", color=GREEN, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(
                rotate(gamma_pp, Y) * (rho1 * cos(gamma_bp * 0.5) * X + rbp * 0.5 * Z),
                direction=-0.3 * Y,
            ),
            MathTex("\\gamma_{p_2}", color=RED, font_size=32)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(
                rotate(-gamma_pw, Y) * (rho1 * cos(gamma_bw * 0.5) * X - rbw * 0.5 * Z),
                direction=-0.3 * normalize(Z - Y)
            ),
        )

        n = 10

        # Arc between O1 and O2
        lines = arc(
            func=lambda t: rotate((gamma_pp + gamma_pw) * t, Y) * O2,
            n=n,
            color=WHITE,
        )

        # Arc between O1 and M1
        lines += arc(
            func=lambda t: rotate(gamma_bp * t, cross(O1, M1)) * O1,
            n=n,
            color=GREEN,
        )

        # Arc between O2 and M2
        lines += arc(
            func=lambda t: rotate(gamma_bw * t, cross(O2, M2)) * O2,
            n=n,
            color=RED,
        )

        # Arc for alpha
        arcs = arc_with_arrows(
            func=lambda t: rotate(-pressure_angle * t, X) * 0.5 * rho1 * Y,
            tex_content="\\alpha",
            direction=0.5 * Y,
            font_size=32,
            theta_rot=pi / 2,
            phi_rot=pi / 2,
        )

        # Arc between O1P and O1M1
        arcs += arc_with_arrows(
            func=lambda t: rotate(-phip * t, O1) * rotate(gamma_pp, Y) * (rho1 * cos(gamma_bp * 0.5) * X + rbp * 0.5 * Z),
            tex_content="\\phi_{p_1}",
            reverse=True,
            direction=0.7 * Z + 0.001 * Y,
            height=0.1,
            font_size=32,
            theta_rot=pi / 2,
            phi_rot=pi / 2,
        )

        # Arc between O2P and O2M2
        arcs += arc_with_arrows(
            func=lambda t: rotate(-phiw * t, O2) * rotate(-gamma_pw, Y) * (rho1 * cos(gamma_bw * 0.5) * X - rbw * 0.5 * Z),
            tex_content="\\phi_{p_2}",
            height=0.1,
            direction=-Z,
            font_size=32,
            theta_rot=pi / 2,
            phi_rot=pi / 2,
        )

        self.add(sphere, circles, dots, tex, *lines, *arcs)


class BevelAnimation(ThreeDScene):

    def construct(self):
        z = z_pinion = 10
        z_wheel = 20
        m = 0.5
        step = pi * m
        ka = 1
        kd = 1.25
        pressure_angle = pi / 9
        pitch_cone_angle = get_pitch_cone_angle(z_pinion, z_wheel)

        gamma_p = pitch_cone_angle  # for convenience
        gamma_b = asin(cos(pressure_angle) * sin(gamma_p))
        cos_b, sin_b = cos(gamma_b), sin(gamma_b)
        rp = z * step / (2 * pi)
        rho1 = rp / sin(gamma_p)
        rho0 = 2 * rho1 / 3
        k = sin(gamma_p) / z
        gamma_r = gamma_p - atan2(2 * kd * k, 1)
        gamma_f = gamma_p + atan2(2 * ka * k, 1)

        # self.set_camera_orientation(phi=110 * DEGREES, theta=-50 * DEGREES, zoom=0.6, focal_distance=1000)
        self.set_camera_orientation(phi=pi / 2, theta=pi / 2, zoom=0.6, focal_distance=1000)

        rb = rho1 * sin(gamma_b)

        rho1_circle = (
            VGroup(
                Circle(rho1, color=WHITE, fill_opacity=0.2)
                .rotate_about_origin(pi / 2 - gamma_b, Y),
                # Arrow3D(O, rho1 * rotate(-gamma_b, Y) *  Z),
                # Arrow3D(rho1 * rotate(-gamma_b, Y) *  Z, rho1 * rotate(-gamma_b, Y) *  Z + rb * Y),
                Dot3D(rho1 * rotate(-gamma_b, Y) * Z, radius=0.16),
            )
        )

        base_circle = (
            Circle(radius=rb, color=WHITE)
            .move_to(rho1 * cos(gamma_b) * Z)
        )
        cone = Cone(rb, rho1 * cos(gamma_b), -Z, fill_color=DARK_BLUE).set_opacity(0.05)
        sphere = Sphere(radius=rho1, fill_color=DARK_GRAY).set_opacity(0.01)

        shapes = [rho1_circle, base_circle, cone, sphere]

        dots = [
            Dot3D(O, radius=0.16),
            Dot3D(rho1 * cos(gamma_b) * Z, radius=0.16),
            Dot3D(rho1 * rotate(-gamma_b, Y) * Z, radius=0.16),
        ]

        lines = [
            DashedLine(-1.5 * rho1 * Z, 1.5 * rho1 * Z, dash_length=0.2),
            DashedLine(rho1 * cos(gamma_b) * Z + rb * X, rho1 * cos(gamma_b) * Z - rb * X, dash_length=0.2),
            DashedLine(O, rho1 * rotate(-gamma_b, Y) * Z, dash_length=0.2),
        ]

        group = VGroup(*shapes, *dots, *lines)

        t0 = pi / 6
        t = 2 * pi / 5

        self.add(group)

        self.wait()
        self.move_camera(phi=110 * DEGREES, theta=-50 * DEGREES, zoom=0.6, focal_distance=1000)
        self.wait()
        self.play(Rotate(rho1_circle, angle=-t0, axis=Z, about_point=O))


        # rho1_circle = rho1_circle.rotate(-t0 - t, Z, about_point=O)

        t0_line = VGroup(
            Line(rho1 * cos(gamma_b) * Z, rho1 * rotate(-t0, Z) * rotate(-gamma_b, Y) * Z),
            Dot3D(rho1 * rotate(-t0, Z) * rotate(-gamma_b, Y) * Z, radius=0.16)
        )

        t1_line = Line(rho1 * cos(gamma_b) * Z, rho1 * rotate(-t - t0, Z) * rotate(-gamma_b, Y) * Z)

        curve = ParametricFunction(
            lambda t: rho1 * rotate(-pi, Z) * spherical_involute(gamma_b, -t0, -t),
            t_range=[0, t],
            color=RED,
        )

        arc_rb = (
            Arc(radius=rb, angle=t, arc_center=rho1 * cos(gamma_b) * Z, color=ORANGE)
            .rotate_about_origin(pi - t - t0, Z)
        )
        arc_rho1 = (
            Arc(radius=rho1, angle=t0, arc_center=O, color=ORANGE)
            .rotate_about_origin(pi / 2 - gamma_b, Y)
            .rotate_about_origin(-t - t0, Z)
            .rotate_about_origin(pi, rotate(-t - t0, Z) * rotate(pi / 2 - gamma_b, Y) * Z)
        )
        sector = (
            Sector(rb, angle=t, fill_color=GREEN, fill_opacity=0.7, arc_center=rho1 * cos(gamma_b) * Z)
            .rotate_about_origin(pi - t - t0, Z)
        )

        theta_dot = Dot3D(rho1 * rotate(-pi, Z) * spherical_involute(gamma_b, -t0, -t), radius=0.16, color=WHITE)

        self.play(Create(t0_line), Create(t1_line))
        self.play(
            Rotate(rho1_circle, angle=-t, axis=Z, about_point=O, run_time=2.5, run_func=linear),
            Create(curve, run_time=2.5, run_func=linear),
            Create(sector, run_time=2.5, run_func=linear),
        )
        self.wait()
        self.play(
            Create(theta_dot, run_func=linear),
            Create(arc_rb, run_func=linear),
            Create(arc_rho1, run_func=linear),
        )
        self.wait()

class BevelStaticInvolute(ThreeDScene):

    def construct(self):
        z = z_pinion = 10
        z_wheel = 20
        m = 0.5
        step = pi * m
        ka = 1
        kd = 1.25
        pressure_angle = pi / 9
        pitch_cone_angle = get_pitch_cone_angle(z_pinion, z_wheel)

        gamma_p = pitch_cone_angle  # for convenience
        gamma_b = asin(cos(pressure_angle) * sin(gamma_p))
        cos_b, sin_b = cos(gamma_b), sin(gamma_b)
        rp = z * step / (2 * pi)
        rho1 = rp / sin(gamma_p)
        rho0 = 2 * rho1 / 3
        k = sin(gamma_p) / z
        gamma_r = gamma_p - atan2(2 * kd * k, 1)
        gamma_f = gamma_p + atan2(2 * ka * k, 1)

        self.set_camera_orientation(phi=110 * DEGREES, theta=-50 * DEGREES, zoom=0.6, focal_distance=1000)

        rb = rho1 * sin(gamma_b)

        gamma_arrow = Arrow3D(O, rho1 * rotate(-gamma_b, Y) *  Z)
        gamma_prime_arrow = Arrow3D(rho1 * rotate(-gamma_b, Y) *  Z, rho1 * rotate(-gamma_b, Y) *  Z + rb * Y)
        rho1_circle = (
            VGroup(
                Circle(rho1, color=WHITE, fill_opacity=0.2)
                .rotate_about_origin(pi / 2 - gamma_b, Y),
                gamma_arrow,
                gamma_prime_arrow,
                Dot3D(rho1 * rotate(-gamma_b, Y) * Z, radius=0.16),
            )
        )

        base_circle = (
            Circle(radius=rb, color=WHITE)
            .move_to(rho1 * cos(gamma_b) * Z)
        )
        cone = Cone(rb, rho1 * cos(gamma_b), -Z, fill_color=DARK_BLUE).set_opacity(0.05)
        sphere = Sphere(radius=rho1, fill_color=DARK_GRAY).set_opacity(0.01)

        shapes = [rho1_circle, base_circle, cone, sphere]

        dots = [
            Dot3D(O, radius=0.16),
            Dot3D(rho1 * cos(gamma_b) * Z, radius=0.16),
            Dot3D(rho1 * rotate(-gamma_b, Y) * Z, radius=0.16),
        ]

        lines = [
            DashedLine(-1.5 * rho1 * Z, 1.5 * rho1 * Z, dash_length=0.2),
            DashedLine(rho1 * cos(gamma_b) * Z + rb * X, rho1 * cos(gamma_b) * Z - rb * X, dash_length=0.2),
            DashedLine(O, rho1 * rotate(-gamma_b, Y) * Z, dash_length=0.2),
        ]

        group = VGroup(*shapes, *dots, *lines)

        t0 = pi / 6
        t = 2 * pi / 5

        rho1_circle = rho1_circle.rotate(-t0 - t, Z, about_point=O)

        t0_line = VGroup(
            Line(rho1 * cos(gamma_b) * Z, rho1 * rotate(-t0, Z) * rotate(-gamma_b, Y) * Z),
            Dot3D(rho1 * rotate(-t0, Z) * rotate(-gamma_b, Y) * Z, radius=0.16)
        )

        t1_line = Line(rho1 * cos(gamma_b) * Z, rho1 * rotate(-t - t0, Z) * rotate(-gamma_b, Y) * Z)

        curve = ParametricFunction(
            lambda t: rho1 * rotate(-pi, Z) * spherical_involute(gamma_b, -t0, -t),
            t_range=[0, t],
            color=RED,
        )

        arc_rb = (
            Arc(radius=rb, angle=t, arc_center=rho1 * cos(gamma_b) * Z, color=ORANGE)
            .rotate_about_origin(pi - t - t0, Z)
        )
        arc_rho1 = (
            Arc(radius=rho1, angle=t0, arc_center=O, color=ORANGE)
            .rotate_about_origin(pi / 2 - gamma_b, Y)
            .rotate_about_origin(-t - t0, Z)
            .rotate_about_origin(pi, rotate(-t - t0, Z) * rotate(pi / 2 - gamma_b, Y) * Z)
        )
        sector = (
            Sector(rb, angle=t, fill_color=GREEN, fill_opacity=0.7, arc_center=rho1 * cos(gamma_b) * Z)
            .rotate_about_origin(pi - t - t0, Z)
        )
        theta_dot = Dot3D(rho1 * rotate(-pi, Z) * spherical_involute(gamma_b, -t0, -t), radius=0.16, color=WHITE)

        rho1_arrow = Arrow3D(O, -rho1 * rotate(-t - t0, Z) * rotate(-gamma_b, Y) *  Z)
        rho1_group = VGroup(
            rho1_arrow,
            MathTex("\\rho_1")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(rho1_arrow)
        )

        t0_tex = (
            MathTex("t_0")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(rho1 * rotate(-t0, Z) * rotate(-gamma_b, Y) * Z, direction=Z)
        )

        t_tex = (
            MathTex("t")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(rho1 * rotate(-t - t0, Z) * rotate(-gamma_b, Y) * Z, direction=Z)
        )

        gamma_tex = (
            MathTex("\\overrightarrow{\\gamma(t)}")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(gamma_arrow)
        )

        gamma_prime = (
            MathTex("\\overrightarrow{\\gamma'(t)}")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(gamma_prime_arrow, direction=RIGHT)
        )

        rb_arrow = Arrow3D(rho1 * cos(gamma_b) * Z, rho1 * cos(gamma_b) * Z - rb * Y)
        rb_group = VGroup(
            rb_arrow,
            MathTex("r_b")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(rb_arrow)
        )

        theta_tex = (
            MathTex("\\theta(t)")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(theta_dot, direction=normalize(X + Y))
        )

        system = VGroup(
            Arrow3D(O, -0.4 * rho1 * X),
            Arrow3D(O, 0.4 * rho1 * Y),
            Arrow3D(O, 0.4 * rho1 * Z),
            MathTex("\\overrightarrow{x}")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(-0.4 * rho1 * X, direction=-X),
            MathTex("\\overrightarrow{y}")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(0.4 * rho1 * Y),
            MathTex("\\overrightarrow{z}")
            .rotate(110 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 50 * DEGREES, Z, about_point=O)
            .next_to(0.4 * rho1 * Z),
        )

        self.add(
            group,
            t0_line,
            t1_line,
            curve,
            sector,
            arc_rb,
            arc_rho1,
            theta_dot,
            rho1_group,
            t0_tex,
            t_tex,
            gamma_tex,
            gamma_prime,
            rb_group,
            theta_tex,
            system,
        )

def spherical_involuteof(pitch_cone_angle:float, t0:float, alpha:float, t:float) -> vec3:
    """
    Return the spherical interference function

    Parameters
    ----------
    t : float
        The angular position
    t0 : float
        The difference phase
    pitch_cone_angle : float
        The pitch cone angle
    alpha : float
        The height angle offset of the rack

    Returns
    -------
        a normalized `vec3`
    """
    cos_p, sin_p = cos(pitch_cone_angle), sin(pitch_cone_angle)
    return (
        cos(alpha) * spherical_involute(pitch_cone_angle, t0, t) 
        + sin(alpha) * vec3(-cos_p * cos(t + t0), -cos_p * sin(t + t0), sin_p)
    )


def derived_spherical_involute(cone_angle:float, t0:float):
    """
    Return the function of the derived spherical involute function.

    Parameters
    ----------
    cone_angle : float
        The cone angle
    t0 : float
        The phase difference
    """
    cos_g, sin_g = cos(cone_angle), sin(cone_angle)
    return lambda t: vec3(
        cos_g ** 2 * sin(t * sin_g) * cos(t + t0),
        cos_g ** 2 * sin(t * sin_g) * sin(t + t0),
        -cos_g * sin_g * sin(t * sin_g),
    )


def jacobian_spherical_involute(base_cona_angle:float, pitch_cone_angle:float, t01:float, t02:float, alpha:float) -> callable:
    """
    Return the function of the jacobian used for the newton method in `spherical_gearprofile`

    Parameters
    ----------
    base_cona_angle : float
        The base cone angle
    pitch_cone_angle : float
        The pitch cone angle
    t01 : float
        The phase of the spherical involute function
    t02 : float
        The phase of the spherical interference function
    alpha : float
        The height angle offset of the rack
    """
    dsi = derived_spherical_involute # for convenience
    derived_involute = dsi(base_cona_angle, t01)
    cos_p = cos(pitch_cone_angle)
    vec = lambda t: vec3(cos_p * sin(t), -cos_p * cos(t), 0)
    derived_interference = lambda t: dsi(pitch_cone_angle, t02)(t) * cos(alpha) + sin(alpha) * vec(t + t02)
    return lambda t1, t2: mat3(derived_involute(t1), -derived_interference(t2), vec3(0, 0, 1))


def spherical_rack_tools(z:float, pressure_angle:float=pi / 9, ka:float=1, kd:float=1.25):
    """
    Return a list of all information useful to generate a spherical rack.

    Parameters
    ----------
    z : float
        Number of tooth of the rack equal to `z_pinion / sin(pitch_cone_angle)`
        or `z_wheel / sin(shaft_angle - pitch_cone_angle)`
    pressure_angle : float
        The pressure angle of the gear
    ka : float
        The addendum coefficient
    kd : float
        The dedendum coefficient

    Returns
    -------
    tuple
        * the minimum abscissa for the function (fifth element)
        * the maximum abscissa for the function (fifth element)
        * the phase of a tooth
        * the phase of space
        * the function to generate a tooth
    """
    k = 1 / z
    gamma_p = 0.5 * pi
    gamma_b = asin(cos(pressure_angle))
    cos_b, sin_b = cos(gamma_b), sin(gamma_b)
    gamma_f = gamma_p + atan2(2 * ka * k, 1)
    gamma_r = gamma_p - atan2(2 * kd * k, 1)

    phi_p = acos(tan(gamma_b) / tan(gamma_p))
    theta_p = atan2(sin_b * tan(phi_p), 1) / sin_b - phi_p
    phase_diff = k * pi + 2 * theta_p

    t_min = acos(cos(gamma_r) / cos_b) / sin_b
    t_max = acos(cos(gamma_f) / cos_b) / sin_b

    involute = lambda t, t0: spherical_involute(gamma_b, t0, t)
    v = vec3(1, 1, 0)
    phase_empty = 2 * pi * k - anglebt(involute(t_min, 0) * v, involute(-t_min, phase_diff) * v)

    return [t_min, t_max, phase_diff, phase_empty, involute]

def spherical_gearprofile(
    z:int, 
    step: float,
    pitch_cone_angle:float, 
    pressure_angle:float = pi / 9, 
    ka:float = 1, 
    kd:float = 1.25, 
):
    """
    Generate 1-period tooth spherical profile for a bevel gear

    Parameters
    ----------
    z : int
        Number of tooth on the gear this profile is meant for
    step : float
        Step of a tooth
    pitch_cone_angle : float
        The pitch cone angle
    pressure_angle : float
        Pressure angle of the tooth
    ka : float
        Addendum coefficient
    kd : float
        Dedendum coefficient
    """
    # Initialization of parameters
    gamma_p = pitch_cone_angle # for convenience
    rp = z * step / (2 * pi)
    rho1 = rp / sin(gamma_p)
    gamma_b = asin(cos(pressure_angle) * sin(gamma_p))
    cos_b, sin_b = cos(gamma_b), sin(gamma_b)
    tooth_size = pi / z
    involute = lambda t, t0 : spherical_involute(gamma_b, t0, t)
    epsilon_p = acos(cos(gamma_p) / cos_b) / sin_b
    theta_p = anglebt(involute(0, 0) * vec3(1, 1, 0), involute(epsilon_p, 0) * vec3(1, 1, 0))
    phase_diff = tooth_size + 2 * theta_p
    phase_empty = phase_interference = 2 * pi / z - phase_diff
    # The following number `k` is useful to simplify some calculations
    # It's broadly speaking `1/z_rack` and `z_rack` is not an integer !
    k = sin(gamma_p) / z

    # Spherical involute part
    gamma_f = gamma_p + atan2(2 * ka * k, 1) # addendum cone angle
    gamma_r = gamma_p - atan2(2 * kd * k, 1) # dedendum cone angle
    t_min = 0
    t_max = acos(cos(gamma_f) / cos_b) / sin_b
    if gamma_r > gamma_b:
        v = vec3(1, 1, 0)
        t_min = acos(cos(gamma_r) / cos_b) / sin_b
        phase_empty = 2 * pi / z - anglebt(
            involute(t_min, 0) * v, involute(-t_min, phase_diff) * v
        )

    # Calculation of offsets due to geometry of spherical rack
    _, t_rack_max, phase1, _, rinvolute = spherical_rack_tools(1 / k, pressure_angle, ka, kd)
    interference = lambda t, t0 : spherical_involuteof(gamma_p, t0, alpha, t)
    alpha = atan2(2 * ka * k, 1)
    n1, n2 = rinvolute(t_rack_max, 0) * vec3(1, 1, 0), rinvolute(-t_rack_max, phase1) * vec3(1, 1, 0)
    beta = 0.5 * anglebt(n1, n2) * length(n1) / sin_b

    # Newton method to calculate the intersection between
    # the spherical involute and the spherical interference.
    # Objective function
    involuteat = lambda t2, t0 : spherical_involuteof(gamma_p, t0, alpha, t2)
    f = lambda t1, t2: involute(t1, 0) - involuteat(t2, -0.5 * phase_interference + beta)
    # Jacobian matrix
    J = jacobian_spherical_involute(gamma_b, gamma_p, 0, -0.5 * phase_interference + beta, alpha)

    # Compute the intersection values
    t1, t2, t3 = 0.5 * t_max, -0.5 * t_max, 0
    for i in range(8):
        t1, t2, t3 = vec3(t1, t2, t3) - inverse(J(t1, t2)) * f(t1, t2)

    # Build sides of a tooth
    interference1 = ParametricFunction(
        lambda t: rho1 * interference(t, -0.5 * phase_interference + beta),
        t_range=[t2, 0],
    )
    interference2 = ParametricFunction(
        lambda t: rho1 * interference(-t, phase_diff + 0.5 * phase_interference - beta),
        t_range=[t2, 0],
    )
    side1 = ParametricFunction(lambda t: rho1 * involute(t, 0), t_range=[t1, t_max])
    side2 = ParametricFunction(lambda t: rho1 * involute(-t, phase_diff), t_range=[t1, t_max])

    # Extreme points of sides to compute angle between them
    a = rho1 * interference(0, -0.5 * phase_interference + beta)
    b = rho1 * interference(0, phase_diff + 0.5 * phase_interference - beta)
    final_phase_empty = 2 * pi / z - anglebt(a * vec3(1, 1, 0), b * vec3(1, 1, 0))
    top = Line(rho1 * involute(t_max, 0), rho1 * involute(-t_max, phase_diff))
    bottom = Line(rotate(-final_phase_empty, vec3(0, 0, 1)) * a, a)

    aligned_phase = anglebt(
        (rho1 * involute(t_max, 0) + rho1 * involute(-t_max, phase_diff)) * 0.5 * vec3(1, 1, 0),
        X,
    )

    return VGroup(
        interference1,
        interference2,
        side1,
        side2,
        top,
        bottom,
        Dot(a, radius=0.03),
        Dot(b, radius=0.03),
        Dot(rho1 * involute(t1, 0), radius=0.03),
        Dot(rho1 * involute(t_max, 0), radius=0.03),
        Dot(rho1 * involute(-t1, phase_diff), radius=0.03),
        Dot(rho1 * involute(-t_max, phase_diff), radius=0.03),
        # Line(rotate(-final_phase_empty, vec3(0, 0, 1)) * a, rho1 * Z),
        # Line(b, rho1 * Z),
    ).rotate_about_origin(-aligned_phase, Z)

class BevelTransmission(ThreeDScene):

    def construct(self):
        z_pinion = 10
        z_wheel = 20
        m = 0.5
        step = pi * m
        pressure_angle = pi / 9
        pitch_cone_angle = get_pitch_cone_angle(z_pinion, z_wheel)

        gamma_p = pitch_cone_angle # for convenience
        rp = z_pinion * step / (2 * pi)
        rho1 = rp / sin(gamma_p)

        angle1tooth_pinion = 2 * pi / z_pinion
        # pinion_profile = spherical_gearprofile(z_pinion, step, pitch_cone_angle)
        pinion_profile = fast.revolution(spherical_gearprofile(z_pinion, step, pitch_cone_angle), angle1tooth_pinion, z_pinion)

        angle1tooth_wheel = 2 * pi / z_wheel
        # wheel_profile = spherical_gearprofile(z_wheel, step, pi / 2 - pitch_cone_angle)
        wheel_profile = fast.revolution(spherical_gearprofile(z_wheel, step, pi / 2 - pitch_cone_angle), angle1tooth_wheel, z_wheel)

        sphere = Sphere(radius=rho1, fill_color=DARK_GRAY).set_opacity(0.01)

        # Pinion
        gamma_pp = pitch_cone_angle
        gamma_bp = asin(cos(pressure_angle) * sin(gamma_pp))
        rpp = z_pinion * step / (2 * pi)
        rho1 = rho1p = rpp / sin(gamma_pp)
        rbp = rho1p * sin(gamma_bp)
        phip = acos(tan(gamma_bp) / tan(gamma_pp))

        # Wheel
        gamma_pw = pi * 0.5 - pitch_cone_angle
        gamma_bw = asin(cos(pressure_angle) * sin(gamma_pw))
        rpw = z_wheel * step / (2 * pi)
        rho1w = rpw / sin(gamma_pw)
        rbw = rho1w * sin(gamma_bw)
        phiw = acos(tan(gamma_bw) / tan(gamma_pw))

        circles = VGroup(
            # Circle(rho1, color=WHITE),
            # Circle(rho1, color=WHITE).rotate(pi / 2, Y, about_point=O),
            # Circle(rho1, color=WHITE).rotate(-pressure_angle, X, about_point=O),
            Circle(rbw, color=RED)
            .rotate_about_origin(pi / 2 - gamma_pw, Y)
            .move_to(rho1 * cos(gamma_bw) * rotate(-gamma_pw, Y) * X),
            DashedVMobject(Circle(rpw, color=RED))
            .rotate_about_origin(pi / 2 - gamma_pw, Y)
            .move_to(rho1 * cos(gamma_pw) * rotate(-gamma_pw, Y) * X),
            Circle(rbp, color=BLUE)
            .rotate_about_origin(pi / 2 + gamma_pp, Y)
            .move_to(rho1 * cos(gamma_bp) * rotate(gamma_pp, Y) * X),
            DashedVMobject(Circle(rpp, color=BLUE))
            .rotate_about_origin(pi / 2 + gamma_pp, Y)
            .move_to(rho1 * cos(gamma_pp) * rotate(gamma_pp, Y) * X),
        )

        lines = VGroup(
            DashedLine(rho1 * X, -0.5 * rho1 * Y + rho1 * X),
            Arrow3D(start=rho1 * X, end=-0.5 * rho1 * Y + rho1 * X, stroke_width=6)
            .rotate(-pressure_angle, X, about_point=rho1 * X),
        )

        O1 = rho1 * rotate(-gamma_pw, Y) * X
        O2 = rho1 * rotate(gamma_pp, Y) * X
        P = rho1 * X
        T1 = rotate(-phiw, O1) * rotate(-gamma_pw, Y) * (rho1 * cos(gamma_bw) * X - rbw * Z)
        T2 = rotate(-phip, O2) * rotate(gamma_pp, Y) * (rho1 * cos(gamma_bp) * X + rbp * Z)

        dots = VGroup(
            Dot3D(P, color=WHITE),
            Dot3D(O1, color=RED),
            Dot3D(T1, color=RED),
            Dot3D(O2, color=BLUE),
            Dot3D(T2, color=BLUE),
        )

        tex = VGroup(
            MathTex("O_1", color=RED)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(O1, direction=-Z),
            MathTex("O_2", color=BLUE)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(O2, direction=Z),
            MathTex("T_1", color=RED)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(T1, direction=Z),
            MathTex("T_2", color=BLUE)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(T2, direction=-Z),
            MathTex("P", color=WHITE)
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(P, direction=Z),
            MathTex("\\overrightarrow{F_{2 / 1}}")
            .rotate(pi / 2, X, about_point=O)
            .rotate(pi / 2, Z, about_point=O)
            .next_to(-rotate(-pressure_angle, X) * 0.5 * rho1 * Y, direction=normalize(Z + Y)),
        )


        # Arc for alpha
        arcs = arc_with_arrows(
            func=lambda t: -rotate(-pressure_angle * t, X) * 0.4 * rho1 * Y,
            tex_content="\\alpha",
            direction=-0.5 * Y,
            theta_rot=pi / 2,
            phi_rot=pi / 2,
        )

        self.set_camera_orientation(phi=90 * DEGREES, theta=00 * DEGREES, zoom=0.8, focal_distance=1000)
        # self.add(sphere, pinion_profile, wheel_profile)
        self.add(
            sphere,
            circles,
            lines,
            dots,
            tex,
            *arcs,
            pinion_profile.rotate_about_origin(pi + pi / (2 * z_pinion), Z).rotate_about_origin(pi / 2 + gamma_p, Y),
            wheel_profile.rotate_about_origin(pi / (2 * z_wheel), Z).rotate_about_origin(gamma_p, Y),
        )
