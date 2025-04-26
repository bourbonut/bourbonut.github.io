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
        z_pinion = 10
        z = z_wheel = 20
        m = 0.3
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
            Line(perp * 0.5, perp * 1.2),
            Line(
                rho1 * rotation(gamma_p) * X + perp * 0.5,
                rho1 * rotation(gamma_p) * X + perp * 1.2,
            ),
            DoubleArrow(
                perp * 0.85,
                rho1 * rotation(gamma_p) * X + perp * 0.85,
                stroke_width=3,
                max_tip_length_to_length_ratio=0.04,
                buff=0.025,
            ),
            MathTex("\\rho_1", font_size=32).move_to(
                0.5 * rho1 * rotation(gamma_p) * X + perp * 1.1
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
                    MathTex(tex, font_size=32, color=color).move_to(
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
            MathTex("k_a \\cdot m", font_size=32, color=BLUE)
            .move_to(rho1 * 1.1 * rotation(gamma_p) * X + ka * m * perp * 2.2),
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
            MathTex("k_d \\cdot m", font_size=32, color=BLUE)
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
            MathTex("r_p", font_size=32, color=YELLOW).move_to(0.96 * rp / tan(gamma_p) * X + rp * 0.5 * Y)
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
    tex_pos: vec3,
    n: int = 10,
    height: float = 0.15,
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
    tex_pos : vec3
        Tex position
    n : int
        Resolution of the arc
    height : float
        Height of arrows
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
    return [
        Line(func(i / n), func((i + 1) / n), color=color)
        for i in range(n)
    ] + [
        Line(func((i + 0.5) / n), func((i + 1.5) / n), color=color)
        for i in range(n - 1)
    ] + [
        Arrow3D(func(1 / n), func(0 / n), resolution=0, height=0.1, color=color),
        Arrow3D(func(9 / n), func(10 / n), resolution=0, height=height, color=color),
        MathTex(tex_content, color=color)
        .rotate(phi_rot, X, about_point=O)
        .rotate(theta_rot, Z, about_point=O)
        .move_to(tex_pos),
    ]

class SphericalRepr(ThreeDScene):

    def construct(self):
        z_pinion = 10
        z = z_wheel = 20
        m = 0.3
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

        self.set_camera_orientation(phi=75 * DEGREES, theta=-20 * DEGREES, zoom=0.7)

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
            .move_to(rho1 * rotate(gamma_b, Y) * (X - 0.03 * Z - 0.06 * Y)),
        )

        P = VGroup(
            Dot3D(positions["P"]),
            MathTex("P")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(rho1 * rotate(phi, Z) * (X + 0.05 * Z + 0.06 * Y)),
        )

        M = VGroup(
            Dot3D(positions["M"]),
            MathTex("M")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(rho1 * (X - 0.05 * Z - 0.06 * Y)),
        )

        Q = VGroup(
            Dot3D(positions["Q"]),
            MathTex("Q")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(rho1 * rotate(-epsilon, rotate(gamma_b, Y) * X) * (X + 0.06 * Y + 0.08 * Z)),
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
            color=GREEN,
        )

        # epsilon (OM, OQ)
        arrows = arc_with_arrows(
            func=lambda t: rotate(-epsilon * t, rotate(gamma_b, Y) * X) * rho1 * rotate(gamma_b * 0.5, Y) * X,
            tex_content="\\epsilon",
            tex_pos=rotate(-epsilon * 0.5, rotate(gamma_b, Y) * X) * rho1 * rotate(gamma_b * 0.6, Y) * X,
            n=n,
            height=0.15,
            color=BLUE,
        )

        # gamma_b (O_S O, O_S M) where O_S is the origin
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rho1 * rotate(gamma_b * t, Y) * X,
            tex_content="\\gamma_b",
            tex_pos=0.4 * rho1 * rotate(gamma_b * 0.5, Y) * X,
            n=n,
            height=0.1,
            color=GREEN,
        )

        # phi (O_S M, O_S P) where O_S is the origin
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rho1 * rotate(phi * t, Z) * X,
            tex_content="\\varphi",
            tex_pos=0.65 * rho1 * rotate(phi * 0.3, Z) * X,
            n=n,
            height=0.1,
            color=RED,
        )

        # epsilon (y2, y3)
        arrows += arc_with_arrows(
            func=lambda t: rotate(-epsilon * t, rotate(gamma_b, Y) * X) * 0.7 * 0.5 * rho1 * rotate(gamma_b, Y) * Z,
            tex_content="\\epsilon",
            tex_pos=rotate(-epsilon * 0.5, rotate(gamma_b, Y) * X) * 0.8 * 0.5 * rho1 * rotate(gamma_b, Y) * Z,
            n=n,
            height=0.15,
            color=BLUE,
        )

        # gamma_b (y0, y2)
        arrows += arc_with_arrows(
            func=lambda t: 0.7 * 0.5 * rho1 * rotate(gamma_b * t, Y) * Z,
            tex_content="\\gamma_b",
            tex_pos=0.8 * 0.5 * rho1 * rotate(gamma_b * 0.5, Y) * Z,
            n=n,
            height=0.1,
            color=GREEN,
        )

        # phi (x0, x1)
        arrows += arc_with_arrows(
            func=lambda t: 0.5 * rotate(phi * t, Z) * 1.1 * rho1 * Y,
            tex_content="\\varphi",
            tex_pos=0.65 * rotate(phi * 0.5, Z) * 1.1 * rho1 * Y,
            n=n,
            height=0.1,
            color=RED,
        )

        rho1_repr = VGroup(
            Arrow3D(O, rho1 * rotate(-pi / 4, Z) * X),
            MathTex("\\rho_1")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(rho1 * 0.5 * rotate(-pi / 4 * 1.3, Z) * X),
        )

        rb_repr = VGroup(
            Arrow3D(cos(gamma_b) * positions["O"], rho1 * rotate(epsilon, rotate(gamma_b, Y) * X) * X),
            MathTex("r_b")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(rho1 * rotate(epsilon, rotate(gamma_b, Y) * X) * rotate(gamma_b * 0.5, Y) * X),
        )

        base_text = (
            Text("Base circle", font_size=32)
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(rho1 * rotate(-gamma_b * 1.7, Z) * rotate(gamma_b * 1.2, Y) * X)
        )

        # Add curve
        ref = rho1 * rotate(-pi / 2 - phi, Z) * rotate(-pi / 2, X) * rotate(pi, Z) * spherical_involute(gamma_p, 0, 0)
        angle = anglebt(positions["M"], ref)
        rot = (
            rotate(pi / 2, positions["Q"]) *
            rotate(-epsilon, positions["O"]) *
            rotate(-pi / 2 - phi - angle, Z) *
            rotate(-pi / 2, X) * rotate(pi, Z)
        )
        curve = ParametricFunction(
            lambda t: rho1 * rot * spherical_involute(gamma_p, 0, t),
            t_range=[0, epsilon * 0.98], # 0.98 due to line width
            color=YELLOW
        )

        system0 = VGroup(
            Arrow3D(O, 1.1 * positions["P"]),
            Arrow3D(O, 1.1 * rotate(pi / 2, Z) * positions["P"]),
            Arrow3D(O, 0.5 * rho1 * Z),
            MathTex("\\overrightarrow{x_0}", font_size=80)
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(1.4 * rotate(pi / 2, Z) * rho1 * rotate(phi, Z) * (X + 0.05 * Y)),
            MathTex("\\overrightarrow{z_0}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(1.15 * positions["P"]),
            MathTex("\\overrightarrow{y_0}, \\overrightarrow{y_1}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(0.5 * rho1 * (Z + 0.7 * X)),
        )
        center = system0.get_center()
        system0 = system0.move_to(1.8 * Z + center)

        system1 = VGroup(
            Arrow3D(O, 1.1 * positions["M"]),
            Arrow3D(O, 1.1 * rotate(pi / 2, Z) * positions["M"]),
            MathTex("\\overrightarrow{x_1}, \\overrightarrow{x_2}", font_size=70)
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(1.35 * rotate(pi / 2, Z) * rho1 * (X + 0.04 * Y)),
            MathTex("\\overrightarrow{z_1}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(1.2 * positions["M"]),
        )
        center = system1.get_center()
        system1 = system1.move_to(1.8 * Z + center).set_color(RED)

        system2 = VGroup(
            Arrow3D(positions["O"], 1.1 * positions["O"]),
            Arrow3D(O, 0.5 * rotate(-pi / 2, Y) * positions["O"]),
            MathTex("\\overrightarrow{y_2}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(0.5 * rotate(-pi / 2, Y) * rho1 * cos(gamma_b) * rotate(gamma_b, Y) * (X + 0.15 * Y)),
            MathTex("\\overrightarrow{z_2}, \\overrightarrow{z_3}")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(1.2 * positions["O"]),
        )
        center = system2.get_center()
        system2 = system2.move_to(1.8 * Z + center).set_color(GREEN)

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
        )
        center = system3.get_center()
        system3 = system3.move_to(1.8 * Z + center).set_color(BLUE)

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
        self.add(group, system0, system1, system2, system3)
