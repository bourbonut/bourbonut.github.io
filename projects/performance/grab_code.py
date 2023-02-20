import requests
from operator import itemgetter
from itertools import starmap

html = """---
title: "Computation of the sum of prime numbers in different languages"
author: "Benjamin Bourbon"
date: "2023-01-20"
categories: [computing]
---

<script src="https://d3js.org/d3.v4.js"></script>

{}

{}

<div id="standard"></div>
<div id="logarithmic"></div>

<script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>
<script src="plot.js"></script>
<script src="selection.js"></script>
"""

raw_url = "https://raw.githubusercontent.com/bourbonut/prime_numbers_sum/main/programs/{}/primes.{}"


def grab_code(language, suffix):
    response = requests.get(raw_url.format(language, suffix))
    content = response.content.decode("utf-8")
    return f"```{language}\n{content}```"

def button(language):
    title = language.title() if language != "cpp" else "C++"
    return f"""<input id="{language}-button" type="button" value="{title}" style="color:white; text-decoration:none; border:none; border-bottom: white 1px solid;" onclick="select_language('{language}');"></input>"""

languages = [
    ("python", "py"),
    ("rust", "rs"),
    ("cpp", "cpp"),
    ("elixir", "ex"),
    ("java", "java"),
    ("lua", "lua"),
    ("javascript", "js"),
    ("scala", "scala"),
]

# print(grab_code("python", "py"))
code = "\n\n".join(starmap(grab_code, languages))
buttons = "\n".join(map(button, map(itemgetter(0), languages)))

content = html.format(buttons, code)
with open("index.qmd", "w") as file:
    file.write(content)

