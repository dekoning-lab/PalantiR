var render_intervals = function(intervals, plot, options) {
    var interval_tip = SVGTip({
        content: function(d) {
            var description = ["Interval:"];
            Object.getOwnPropertyNames(d)
            .forEach(function(name, index, array) {
                if (name !== "color") {
                    var munged_name = name.split("_").join(" ");
                    description.push(munged_name + ": " + d[name]);
                }
            });
            return description;
        },
        parent: plot.svg
    });
    var color_scale = d3.scaleOrdinal(d3.schemeCategory10);
    var intervals = plot.nodes.selectAll(".interval")
        .data(intervals)
        .enter()
        .append("line")
        .filter(function(d){
            var node = d3.select(this.parentNode).datum();
            return node.index === d.node;
        })
        // Should we filter on sites?
        .attrs({
            "class": "interval",
            y1: 0,
            y2: 0,
            x1: function(d) {
                var node = d3.select(this.parentNode).datum();
                return plot.scale(d.from - node.length);
            },
            x2: function(d) {
                var node = d3.select(this.parentNode).datum();
                return plot.scale(d.to - node.length);
            }
        })
        .styles({
            "stroke": function(d) {
                return color_scale(d.state);
            },
            "stroke-width": options.line_width,
            "stroke-opacity": options.interval_opacity
        })
        .call(interval_tip);
};
