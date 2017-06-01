var right_angle_diagonal = function() {
    // Produce path data with right angles
    var projection = function(d) { return [d.y, d.x]; };

    var path = function(path_data) {
        return "M" + path_data[0] + " " + path_data[1] + " " + path_data[2];
    };

    function diagonal(diagonal_path, i) {
        var source = diagonal_path.source,
        target = diagonal_path.target,
        midpoint_x = (source.x + target.x) / 2,
        midpoint_y = (source.y + target.y) / 2,
        path_data = [source, {x: target.x, y: source.y}, target];
        path_data = path_data.map(projection);
        return path(path_data);
    }

    diagonal.projection = function(x) {
        if (!arguments.length) return projection;
        projection = x;
        return diagonal;
    };

    diagonal.path = function(x) {
        if (!arguments.length) return path;
        path = x;
        return diagonal;
    };

    return diagonal;
}();

var render_tree = function(plot, root) {
    var options = plot.options;

    // Initialize plot
    d3.select("#" + plot.container.id).selectAll("svg").remove();

    plot.svg = d3.select("#" + plot.container.id)
        .append("svg")
        .attrs({
            "width": options.width,
            "height": options.height,
            "class": "phylogram"
        });

    plot.vis =  plot.svg.append("g")
        .attr("transform", "translate(" + options.padding + "," + options.padding + ")");

    // Initialize data
    plot.hierarchy = d3.hierarchy(root);

    var layout = d3.cluster()
        .size([options.height - (3 * options.padding), options.width - (2 * options.padding)])
        .separation(function(a, b) { return 1; });

    layout(plot.hierarchy);

    var nodes = plot.hierarchy.descendants();
    var links = plot.hierarchy.links(nodes);
    var leaves = plot.hierarchy.leaves();

    // Assign helper attributes to nodes
    plot.hierarchy.eachBefore(function(node) {
        node.data.dist = (node.parent ? node.parent.data.dist : 0) + (node.data.length || 0);
        var box = string_bounding_box(node.data.label, "Helvetica", options.label_font_size, plot.svg);
        node.data.text_size = (node.data.label ? box.width : 0);
    });

    // Calculate scale, accounting for labels on leaves
    var max_dist = d3.max(attr_getter(nodes, "dist"));
    var max_text_size = options.show_labels ? d3.max(attr_getter(nodes, "text_size")) : 0;

    plot.scale = d3.scaleLinear()
        .domain([ 0, max_dist ])
        .range([ 0, options.width - (2 * options.padding) - max_text_size ]);

    // Reassign y location for each node
    plot.hierarchy.each(function(node) { node.y = plot.scale(node.data.dist); });

    // Create axis
    var inverse_scale = d3.scaleLinear()
        .domain(plot.scale.domain())
        .range(swap(plot.scale.range()));

    var axis = d3.axisBottom()
        .scale(inverse_scale)
        .ticks(options.axis_ticks)
        .tickSize(options.height - (3 * options.padding))
        .tickFormat(d3.format(".2f"));

    // Plot axis
    if (options.show_axis) {
        plot.vis.append("g")
            .attr("class", "axis")
            .call(axis)
            .selectAll("text")
            .styles({
                "font-family": "monospace",
                "font-size": options.axis_tick_font_size,
                "fill": "#888"
        });
        var axis_tick_style = {
                "fill": "none",
                "stroke": "#ddd",
                "shape-rendering": "crispEdges"
        };
        plot.vis.selectAll(".tick line").attrs(axis_tick_style);
        plot.vis.selectAll(".axis path").attrs(axis_tick_style);
    }

    // Plot tree

    // Branches
    var plot_links = plot.vis.selectAll(".link")
        .data(links)
        .enter()
        .append("path")
        .attrs({
            "class": "link",
            "stroke-width": options.line_width,
            "stroke": "#ccc",
            "fill": "none",
            "d": right_angle_diagonal
        })
        .call(SVGTip({
            content: function(d) {
                return [
                    "Branch:",
                    "index: " + d.target.data.index,
                    "length: " + d.target.data.length
                ];
            },
            parent: plot.svg
        }));

    // Nodes - will contain labels, substitutions, and intervals
    var plot_nodes = plot.vis.selectAll(".node")
        .data(nodes)
        .enter()
        .append("g")
        .attrs({
            "class": "node",
            "id": function(d) { return "node_" + d.data.index; },
            "transform": function(d) { return "translate(" + d.y + "," + d.x + ")"; },
            "font-family": "Helvetica, sans-serif"
        });

    // Leaf labels
    if (options.show_labels) {
        plot_nodes.append("text")
        .attrs({
            "dx": function(d) { return d.children ? -8 : 8; },
            "dy": 3
        })
        .styles({
            "font-size": options.label_font_size
        })
        .text(function(d) { return d.children ? null : d.data.label; });
    }

    // Downlod svg by pressing `S`
    document.onkeydown = function (event) {
        if (event.key === "s") download_svg("phylogram", "PalantiR_phylogeny.svg");
    };
};
