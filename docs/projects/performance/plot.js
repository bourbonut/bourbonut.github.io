// Set the dimensions and margins of the graph
var margin = { top: 10, right: 70, bottom: 60, left: 70 },
    width = 700 - margin.left - margin.right,
    height = 630 - margin.top - margin.bottom;

// Append the svg object to the body of the page
var svg_standard = d3.select("#standard")
    .append("svg")
    .attr("viewBox", "0 0 " + (width + margin.left + margin.right) + " " + (height + margin.top + margin.bottom))
    .attr("preserveAspectRatio", "xMinYMin meet")
    // .attr("width", width + margin.left + margin.right)
    // .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform",
        "translate(" + margin.left + "," + margin.top + ")");

// Append the svg object to the body of the page
var svg_logarithmic = d3.select("#logarithmic")
    .append("svg")
    .attr("viewBox", "0 0 " + (width + margin.left + margin.right) + " " + (height + margin.top + margin.bottom))
    .attr("preserveAspectRatio", "xMinYMin meet")
    // .attr("width", width + margin.left + margin.right)
    // .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform",
        "translate(" + margin.left + "," + margin.top + ")");

function zip(arrays) {
    var keys = [];
    for (key in arrays) { keys.push(key) };
    var arrays = keys.map(function(key) { return arrays[key] })
    return Array.apply(null, Array(arrays[0].length)).map(function(_, i) {
        return arrays.map(function(array) { return array[i] })
    });
}

function ceilInt(number){
    if (Number(number) === number && number % 1 === 0){
        var exposant = String(number).length - 1;
        return (Math.floor(number * 10 ** (-exposant)) + 1) * 10 ** exposant;
    } else {
        return number;
    }
}

function draw(svg, data, data_ready, keys, yfunc){
    // A color scale: one color for each group
    var randomcolors = d3.scaleOrdinal().domain(keys).range(d3.schemeSet2);

    var xmax = data["python"].x.map(x => +x).reduce((a, b) => Math.max(a, b));
    // Add X axis --> it is a date format
    var x = d3.scaleLinear().domain([0, xmax]).range([0, width]);
    svg.append("g")
        .attr("class", "xaxis")
        .attr("transform", "translate(0," + height + ")")
        .attr("stroke", "white")
        .call(d3.axisBottom(x))

    svg.append("text")
        .attr("class", "x label")
        .attr("text-anchor", "middle")
        .attr("x", width / 2)
        .attr("y", height + 30)
        .text("Power")
        .attr("fill", "white");

    var ymax = keys.map(key => yfunc(data[key])).flat(1).reduce((a, b) => Math.max(a, b));
    // Add Y axis
    var y = d3.scaleLinear().domain([0, ceilInt(ymax)]).range([height, 0]).nice();
    svg.append("g")
        .attr("class", "yaxis")
        .attr("stroke", "white")
        .call(d3.axisLeft(y).ticks(10, "s"));

    svg.selectAll("path.domain").attr("stroke", "white");
    svg.selectAll("g.tick").selectAll("line").attr("stroke", "white");
    svg.selectAll("g.tick").selectAll("text").attr("fill", "white").attr("stroke", "none");

    svg.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "middle")
        .attr("x", -height / 2)
        .attr("y", -2 * margin.left / 3)
        .attr("transform", "rotate(-90)")
        .text("Execution Time")
        .attr("fill", "white");

    // Add the lines
    var line = d3.line().x(function(d) { return x(+d.x) }).y(function(d) { return y(+yfunc(d)) })
    svg.selectAll("lines").data(data_ready).enter().append("path")
        .attr("class", function(d) { return d.name })
        .attr("d", function(d) { return line(d.values) })
        .attr("stroke", function(d) { return randomcolors(d.name) })
        .style("stroke-width", 4)
        .style("fill", "none")

    // Add the points
    svg
        // First we need to enter in a group
        .selectAll("dots")
        .data(data_ready)
        .enter()
        .append('g')
        .style("fill", function(d) { return randomcolors(d.name) })
        .attr("class", function(d) { return d.name })
        // Second we need to enter in the 'values' part of this group
        .selectAll("myPoints")
        .data(function(d) { return d.values })
        .enter()
        .append("circle")
        .attr("cx", function(d) { return x(d.x) })
        .attr("cy", function(d) { return y(yfunc(d)) })
        .attr("r", 5)
        .attr("stroke", "white")

    // Add a label at the end of each line
    svg
        .selectAll("myLabels")
        .data(data_ready)
        .enter()
        .append('g')
        .append("text")
        .attr("class", function(d) { return d.name })
        .datum(function(d) { return { name: d.name, value: d.values[d.values.length - 1] }; }) // keep only the last value of each time series
        .attr("transform", function(d) { return "translate(" + x(d.value.x) + "," + y(yfunc(d.value)) + ")"; }) // Put the text at the position of the last point
        .attr("x", 12) // shift the text a bit more right
        .text(function(d) { return d.name; })
        .style("fill", function(d) { return randomcolors(d.name) })
        .style("font-size", 15)

    var visible = keys.reduce((visible, key) => {
        visible[key] = true;
        return visible;
    }, {});
    
    function update(d){
        var name = d.name;
        visible[name] = !visible[name]; 
        var ymax = keys.filter(key => visible[key]).
                        map(key => yfunc(data[key])).
                        flat(1).
                        reduce((a, b) => Math.max(a, b));
        y.domain([0, ceilInt(ymax)]).nice();
        svg.selectAll("g.yaxis").
            transition().
            duration(1000).
            call(d3.axisLeft(y).ticks(10, "s"));
        svg.selectAll("g.tick").selectAll("line").attr("stroke", "white");

        var line = d3.line().x(function(d) { return x(d.x) }).y(function(d) { return y(yfunc(d)) })
        for (key in data){
            svg.selectAll("path."+ key)
                .transition()
                .duration(1000)
                .attr("d", function(d) { return line(d.values)})
                .style("opacity", visible[key] ? 1 : 0);

            svg.selectAll("text." + key)
                .transition()
                .duration(1000)
                .attr("transform", function(d) { return "translate(" + x(d.value.x) + "," + y(yfunc(d.value)) + ")"; })
                .attr("opacity", visible[key] ? 1 : 0);

            svg.selectAll("g." + key)
                .selectAll("circle")
                .transition()
                .duration(1000)
                .attr("cy", function(d){return y(yfunc(d))})
                .style("opacity", visible[key] ? 1 : 0);
        } 
    }

    var nbcolumns = 4;
    var middle = (width - 100 * nbcolumns) / 2;
    var rectsize = 20;

    // Add a legend (interactive)
    svg
        .selectAll("legend")
        .data(data_ready)
        .enter()
        .append('g')
        .append('rect')
        .attr('x', function(d, i) { return middle + (i % nbcolumns) * 100 })
        .attr('y', function(d, i) { return 16 + 30 * (Math.floor(i / nbcolumns))})
        .attr('width', rectsize)
        .attr('height', rectsize)
        .style("fill", function(d) {return randomcolors(d.name)})
        .on("click", update)

    svg
        .selectAll("squares")
        .data(data_ready)
        .enter()
        .append('g')
        .append("text")
        .attr('x', function(d, i) { return middle + rectsize + 5 + (i % nbcolumns) * 100 })
        .attr('y', function(d, i) { return 30 * (Math.floor(i / nbcolumns) + 1)})
        .text(function(d) { return d.name; })
        .style("fill", function(d) { return "white" })
        .style("font-size", 15)
        .on("click", update)
}

d3.json("https://raw.githubusercontent.com/bourbonut/prime_numbers_sum/main/plot/allresults.json", function(data) {
    var keys = [];
    for (key in data) { keys.push(key) };

    var data_ready = keys.map(function(key) {
        return {
            name: key,
            values: zip(data[key]).map(function(e) {
                return { x: +e[0], y: +e[1], log: +e[2] };
            })
        };
    })

    draw(svg_standard, data, data_ready, keys, function(d){return d.y});
    draw(svg_logarithmic, data, data_ready, keys, function(d){return d.log});
})
