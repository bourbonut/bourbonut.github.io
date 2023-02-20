var yellow = "#ffd36a";
var white = "white";
var current_selected = "python";

languages = {
    "python": 1,
    "rust": 2,
    "cpp": 3,
    "elixir": 4,
    "java": 5,
    "lua": 6,
    "javascript": 7,
    "scala": 8,
}

document.getElementById(current_selected + "-button").style.borderBottomColor = yellow;
for (language in languages){
    if (language != current_selected){
        document.getElementById("cb" + languages[language]).hidden = true;
    }
}


function select_language(language) {
    if (language != current_selected) {
        document.getElementById("cb" + languages[current_selected]).hidden = true;
        document.getElementById(current_selected + "-button").style.borderBottomColor = white;
        document.getElementById("cb" + languages[language]).hidden = false;
        document.getElementById(language + "-button").style.borderBottomColor = yellow;
        current_selected = language;
    }
}
