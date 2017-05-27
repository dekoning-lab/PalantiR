var render_intervals = function(plot, intervals) {
    var tip = tooltip(plot, "Interval:");

    var color_scale = d3.scaleOrdinal(d3.schemeCategory10);

    plot.hierarchy.eachBefore(function(node) {
        var node_intervals = intervals
            .filter(function(i) {
                return i.node === node.data.index;
            });

            plot.svg.select("#node_" + node.data.index)
                .selectAll(".interval")
                .data(node_intervals)
                .enter()
                .append("line")
                .attrs({
                    "class": "interval",
                    y1: 0,
                    y2: 0,
                    x1: function(d) {
                        return plot.scale(d.from - node.data.length);
                    },
                    x2: function(d) {
                        return plot.scale(d.to - node.data.length);
                    },
                    "stroke": function(d) {
                        return color_scale(d.state);
                    },
                    "stroke-width": plot.options.line_width,
                    "stroke-opacity": plot.options.interval_opacity
                })
                .call(tip);
    });
};
