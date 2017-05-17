var SVGTip = function (options) {
    var get_default = function (x, d) {
        return x === undefined ? d : x;
    };

    var is_array = function (value) {
        return value &&
        typeof value === "object" &&
        typeof value.length === "number" &&
        typeof value.splice === "function" &&
        !(value.propertyIsEnumerable("length"));
    };

    var translate = function (coord) {
        return "translate(" + coord[0] + "," + coord[1] + ")";
    };

    var content = get_default(options.content, "");
    var style = get_default(options.style, "");
    var parent = get_default(options.parent, d3.select("svg"));
    var padding = 7; //px
    var corner_radius = 5; //px
    var line_height = 1.2; //em

    var tooltip = parent.append("g")
        .attrs({
            "class": "svg-tooltip",
            "text-align": "center",
            "font-size": "12px",
            "font-family": "Helvetica",
            "fill": "black",
            "opacity": "0.8",
            "pointer-events": "none"
        });
    var tooltip_style = {
        "fill": "black",
        "opacity": "0.8",
        "corner-radius": corner_radius
    };

    tooltip.container = tooltip.append("g");
    tooltip.box = tooltip.container.append("rect")
        .attrs({
            "class": "box",
            "rx": corner_radius,
            "ry": corner_radius
        })
        .attrs(tooltip_style);

    tooltip.text = tooltip.append("text")
        .attrs({
            "class": "text"
        })
        .attrs({
            "font-size": "12px",
            "font-family": "Helvetica",
            "fill": "#eee",
            "stroke": "none",
            "text-anchor": "middle",
            "alignment-baseline": "middle"
        });

    var show = function () {
        tooltip.style("display", "inline");
    };

    var hide = function () {
        tooltip.style("display", "none");
    };

    return function (selection) {
        selection
            .on("mouseover", show)
            .on("mouseout", hide)
            .on("mousemove", function (d, n) {
                var point = d3.mouse(parent.node());
                point.x = point[0];
                point.y = point[1];

                tooltip.text.selectAll("tspan").remove();

                var lines = ((typeof content === "function") ? content(d, n) : content);
                if (!is_array(lines)) lines = [lines];
                var s = ((typeof style === "function") ? style(d, n) : style);

                lines.forEach(function (line) {
                    tooltip.text.append("tspan").text(line)
                    .attrs({
                        x: 0,
                        dy: line_height + "em"
                    })
                    .attrs({
                        "text-anchor": "middle",
                        "alignment-baseline": "middle"
                    });
                    tooltip.text.attrs(s);
                });

                var bbox = tooltip.text.node().getBBox();

                var x_adj = (bbox.width / 2) + padding;
                var y_adj = - bbox.height - bbox.y - padding;

                console.log(bbox);
                console.log(point);

                if((point.x + bbox.width + padding) > parent.node().width.baseVal.value ||
                    (point.y - bbox.height - padding) < 0) {
                    x_adj = -x_adj;
                    y_adj = 0;
                }
                tooltip.attr("transform",
                    translate([
                        point.x + x_adj,
                        point.y + y_adj
                        ]));

                tooltip.box.attrs({
                    x: bbox.x - padding,
                    y: bbox.y - padding,
                    width: bbox.width + (padding * 2),
                    height: bbox.height + (padding * 2)
                });
            });
    };
};
