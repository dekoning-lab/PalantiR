var render_substitutions = function(plot, substitutions) {
    var tip = tooltip(plot, "Substitution:");

    plot.hierarchy.eachBefore(function(node) {
        var node_substitutions = substitutions
            .filter(function(s) {
                return s.node === node.data.index;
            })
            .filter(function(s) {
                if (plot.options.sites === "all") { return true; }
                if (plot.options.sites === "none") { return false; }
                return plot.options.sites.indexOf(s.site) !== -1;
            });

        plot.svg.select("#node_" + node.data.index)
            .selectAll(".substitution")
            .data(node_substitutions)
            .enter()
            .append("circle")
            .attrs({
                "class": "substitution",
                "r": plot.options.circle_size,
                "cx": function(d) {
                    return plot.scale(d.time - node.data.length);
                },
                "fill": function(d) {
                    if (d.hasOwnProperty("color")) return(d.color);
                    else return('#FC571F');
                }
            })
            .styles({
                "fill-opacity": plot.options.circle_opacity
            })
            .call(tip);
    });
};
