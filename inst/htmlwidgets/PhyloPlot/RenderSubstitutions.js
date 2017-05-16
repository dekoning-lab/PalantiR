var render_substitutions = function(plot, substitutions) {
    var tip = substitution_tooltip(plot);

    plot.hierarchy.eachBefore(function(node) {
        var node_substitutions = substitutions.filter(function(s) {
            return s.node === node.data.index;
        });

        d3.select("#node_" + node.data.index)
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
