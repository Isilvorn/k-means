// defining the size and padding for the scatterplot
var w = 900,
    h = 600,
    pad = 20,
    left_pad = 100;

// selecting the correct portion of the page to print the scatterplot 
var svg = d3.select("#scatterplot")
        .append("svg")
        .attr("width", w)
        .attr("height", h);

var plot_parms = parms(); 
var x = d3.scale.linear().domain([plot_parms[0], plot_parms[1]]).range([left_pad, w-pad]),
    y = d3.scale.linear().domain([plot_parms[2], plot_parms[3]]).range([h-pad*2, pad]);
 
var xAxis = d3.svg.axis().scale(x).orient("bottom")
        .ticks(10)
        .tickFormat(function (fpr, i) {
            return (fpr.toFixed(2));
        }),
    yAxis = d3.svg.axis().scale(y).orient("left")
        .ticks(10)
        .tickFormat(function (tpr, i) {
            return (tpr.toFixed(2));
        });
 
svg.append("g")
    .attr("class", "axis")
    .attr("transform", "translate(0, "+(h-pad)+")")
    .call(xAxis);
 
svg.append("g")
    .attr("class", "axis")
    .attr("transform", "translate("+(left_pad-pad)+")", 0)
    .call(yAxis);
 
svg.append("text")
    .attr("class", "loading")
    .text("Loading ...")
    .attr("x", function () { return w/2; })
    .attr("y", function () { return h/2-5; });
 
var r = d3.scale.linear()
          .domain([0, 1.0])
          .range([3, 8]);

function col(idx) {
    switch (idx) {
        case 1:  return "#0000FF"; break; // blue
        case 2:  return "#32CD32"; break; // lime green
        case 3:  return "#FF0000"; break; // red
        case 4:  return "#FFFF00"; break; // yellow
        case 5:  return "#8A2BE2"; break; // blue violet
        case 6:  return "#FFA500"; break; // orange
        case 7:  return "#008080"; break; // teal
        case 8:  return "#00FFFF"; break; // cyan
        case 9:  return "#D2B48C"; break; // tan
        case 10: return "#FF1493"; break; // deep pink
        case 11: return "#800000"; break; // maroon
        case 12: return "#C0C0C0"; break; // silver
        case 13: return "#BDB76B"; break; // dark khaki
        case 14: return "#FF4500"; break; // orange red
        case 15: return "#FFB6C1"; break; // light pink
        default: return "#000000"; break; // black
    }
}

var scatter_data = setdata();
svg.selectAll(".loading").remove();
 
svg.selectAll("circle")
    .data(scatter_data)
    .enter()
    .append("circle")
    .attr("class", "circle")
    .attr("cx", function (d)    { return x(d[0]);   })
    .attr("cy", function (d)    { return y(d[1]);   })
    .transition()
    .duration(800)
    .attr("r", function (d)     { return r(0.5);      })
    .attr("stroke", function(d) { return col(d[2]); })
    .attr("fill", function (d)  { return col(d[2]); });

svg.append("text")
    .attr("x", 130)
    .attr("y", 60)
    .text("K = " + plot_parms[4])
    .attr("font-size", "30px")
    .attr("font-family", "sans-serif")

svg.append("text")
    .attr("x", 130)
    .attr("y", 100)
    .text("SSE: " + plot_parms[5])
    .attr("font-size", "20px")
    .attr("font-family", "sans-serif")
