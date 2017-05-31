
var render_pair_substitutions = function(plot, substitutions) {
    var tip = tooltip(plot, "Pair substitution:");

    plot.hierarchy.eachBefore(function(node) {
        var node_substitutions = substitutions.filter(function(s) {
           return s.node === node.data.index;
        })
       .filter(function(s) {
            if (plot.options.sites === "all") { return true; }
            if (plot.options.sites === "none") { return false; }
            return plot.options.sites.indexOf(s.site) !== -1;
        });
        var subs = plot.svg.select("#node_" + node.data.index)
            .selectAll(".substitution")
            .data(node_substitutions)
            .enter()
            .append("g")
            .attrs({
                "class": "substitution",
                "transform": function(d) {
                    return "translate(" + plot.scale(d.time - node.data.length) + ",0)";
                    //return plot.scale(d.time - node.data.length);
                }
            })
            .styles({
                "fill-opacity": plot.options.circle_opacity
            })
            .call(tip);

        subs.append("circle")
            .attrs({
                "r": plot.options.circle_size,
                "cy": -plot.options.circle_size,
                "fill": function(d) {
                    if (d.hasOwnProperty("first_color")) return(d.first_color);
                    else return("#FC571F");
                }
            });

        subs.append("circle")
            .attrs({
                "r": plot.options.circle_size,
                "cy": plot.options.circle_size,
                "fill": function(d) {
                    if (d.hasOwnProperty("second_color")) return(d.second_color);
                    else return("#FC571F");
                }
            });
    });

};
