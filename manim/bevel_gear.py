from manim import *
from glm import *
from utils import *
from utils.maths import O, X, Y, Z, rotation
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
    """
    return atan2(sin(shaft_angle), ((z_wheel / z_pinion) + cos(shaft_angle)))


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
        # arc = Arc(rho1, -pi / 2, angle=pi).rotate(pi / 2, Y, about_point=O).rotate(-pi / 6, Z, about_point=O)
        # self.add(arc)

        rb = rho1 * sin(gamma_b)
        base_circle = (
            Circle(radius=rb, color=WHITE)
            .rotate(pi / 2, Y, about_point=O)
            .move_to(rho1 * cos(gamma_b) * X)
            .rotate(gamma_b, Y, about_point=O)
        )

        positions = {
            "O": rho1 * cos(gamma_b) * rotate(gamma_b, Y) * X,
            "P": rho1 * rotate(phi, Z) * X,
            "M": rho1 * X,
            "Q": rho1 * rotate(-epsilon, rotate(gamma_b, Y) * X) * X,
        }

        center = VGroup(
            Dot3D(positions["O"]),
            MathTex("O")
            .rotate(75 * DEGREES, X, about_point=O)
            .rotate(pi / 2 - 20 * DEGREES, Z, about_point=O)
            .move_to(rho1 * cos(gamma_b) * rotate(gamma_b, Y) * (X - 0.07 * Z - 0.03 * Y)),
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
            Line3D(positions["O"], positions["Q"]),
            Line3D(positions["O"], positions["M"]),
            DashedLine(O, positions["O"], dash_length=0.2, color=GREEN),
            DashedLine(O, rotate(-epsilon, positions["O"]) * rotate(pi / 2, Z) * positions["M"], dash_length=0.2, color=BLUE),
        ]

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
            Arrow3D(
                rotate(-epsilon, positions["O"]) * rotate(pi / 2, Z) * positions["M"],
                1.1 * rotate(-epsilon, positions["O"]) * rotate(pi / 2, Z) * positions["M"]
            ),
        )
        center = system3.get_center()
        system3 = system3.move_to(1.8 * Z + center).set_color(BLUE)

        group = VGroup(
            *shapes,
            *points,
            *lines,
            rho1_circle,
            base_circle,
        ).move_to(0.2 * Z)
        self.add(group, system0, system1, system2, system3)
