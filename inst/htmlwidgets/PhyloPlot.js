HTMLWidgets.widget({

    name: "PhyloPlot",

    type: "output",

    factory: function(el, width, height) {

        var data = {
            tree: null,
            substitutions: null,
            intervals: null
        };

        var plot = {
            tree: null,
            nodes: null,
            links: null,
            scale: null,
            vis: null,
            svg: null,
            options: null
        };

        var is_empty_object = function(object) {
            return Object.keys(object).length === 0 && object.constructor === Object;
        };

        var attr_getter = function(array, attr) {
            return array.map(function(n) { return n.data[attr]; });
        };

        var right_angle_diagonal = function() {
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
        };

        var translate = function(x, y) { return "translate(" + x + "," + y + ")"; };

        var swap = function(pair) {return [pair[1], pair[0]]};

        var verify_options = function(options) {
            var get_default = function (x, d) {
                return x === undefined ? d : x;
            };
            options.sites = get_default(options.sites, "all");

            // Promote to array if one element is passed
            if (options.sites.length === undefined) {
                options.sites = [options.sites];
            }

            options.show_substitutions = get_default(options.show_substitutions, true);
            options.show_labels = get_default(options.show_labels, true);
            options.plot_intervals = get_default(options.plot_intervals, false);
            options.interval_opacity = get_default(options.interval_opacity, 0.3);
            options.label_font_size = get_default(options.label_font_size, 10);
            options.show_axis = get_default(options.show_axis, true);
            options.axis_ticks = get_default(options.axis_ticks, 20);
            options.axis_tick_font_size = get_default(options.axis_tick_font_size, 10);
            options.padding = get_default(options.padding, 10);
            options.circle_size = get_default(options.circle_size, 4);
            options.circle_opacity = get_default(options.circle_opacity, 1);
            options.line_width = get_default(options.line_width, 4);
        };

        var render_tree = function(root, options, width, height) {
            var cluster_layout = d3.cluster()
                .size([height - (2 * options.padding), width - (2 * options.padding)])
                .separation(function(a, b) { return 1; });

            var root_node = d3.hierarchy(root);
            cluster_layout(root_node);
            var nodes = root_node.descendants();
            var links = root_node.links(nodes);
            var leaves = root_node.leaves();

            root_node.eachBefore(function(node) {
                node.data.dist = (node.parent ? node.parent.data.dist : 0) + (node.data.length || 0);
                node.data.text_size = (node.data.label ? (node.data.label.length * monospace_char_width(options.label_font_size)) : 0);
            });

            var max_dist = d3.max(attr_getter(nodes, "dist"));
            var max_text_size = options.show_labels ? d3.max(attr_getter(nodes, "text_size")) : 0;

            plot.scale = d3.scaleLinear()
                .domain([ 0, max_dist ])
                .range([ 0, width - (2 * options.padding) - max_text_size ]);

            root_node.each(function(node) {
                node.y = plot.scale(node.data.dist);
            });


            var inverse_scale = d3.scaleLinear()
                .domain(plot.scale.domain())
                .range(swap(plot.scale.range()));

            var axis = d3.axisBottom()
                .scale(inverse_scale)
                .ticks(options.axis_ticks)
                .tickSize(height - (3 * options.padding))
                .tickFormat(d3.format("2r"));

            var h = options.padding;
            var spacing = (height - (options.padding * 4)) / leaves.length;

            var diagonal = right_angle_diagonal();

            document.onkeydown = function (event) {
                if (event.key === "s") download_svg("phylogram", "PalantiR_phylogeny.svg");
            };

            d3.select("#" + el.id).selectAll("svg").remove();

            plot.svg = d3.select("#" + el.id)
                .append("svg")
                .attrs({
                    "width": width,
                    "height": height,
                    "id": "phylogram"
                })
                .append("g")
                .attr("transform", translate(options.padding, options.padding));

            plot.vis = plot.svg.append("g");

            if (options.show_axis) {
                plot.vis.append("g")
                .attr("class", "axis")
                .call(axis)
                .selectAll("text")
                .styles({
                    "font-family": "monospace",
                    "font-size": options.axis_tick_font_size,
                    "fill": "#ddd"
                });
                plot.vis.selectAll(".axis")
                .styles({
                    "fill": "none",
                    "stroke": "#ddd",
                    "shape-rendering": "crispEdges"
                });
            }

            plot.links = plot.vis.selectAll(".link")
            .data(links)
            .enter()
            .append("path")
            .attrs({
                "class": "link",
                "stroke-width": options.line_width,
                "stroke": "#ccc",
                "fill": "none",
                "d": diagonal
            })
            .call(SVGTip({
                content: function(d) {
                    return ["Branch:",
                    "index: " + d.target.data.index,
                    "length: " + d.target.data.length];
                },
                parent: plot.svg
            }));

            plot.nodes = plot.vis.selectAll(".node")
            .data(nodes)
            .enter()
            .append("g")
            .attrs({
                "class": "node",
                "id": function(d) { return "node_" + d.data.index; },
                "transform": function(d) { return translate(d.y, d.x); },
                "font-family": "monospace"
            });

            if (options.show_labels) {
                plot.nodes.append("text")
                .attrs({
                    "dx": function(d) { return d.children ? -8 : 8; },
                    "dy": 3
                })
                .styles({
                    "font-size": options.label_font_size

                })
                .text(function(d) { return d.children ? null : d.data.label; });
            }
            return root_node;
        };

        var render_substitutions = function(root_node, substitutions, plot, options) {
            var subs_tip = SVGTip({
                content: function(d) {
                    var site = (options.n_sites === 1 ? 1 : d.site);
                    var description = ["Substitution:"];
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

            root_node.eachBefore(function(node) {
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
                        "r": options.circle_size,
                        "cx": function(d) {
                            return plot.scale(d.time - node.data.length);
                        },
                        "fill": function(d) {
                            if (d.hasOwnProperty("color")) return(d.color);
                            else return('#FC571F');
                        }
                    })
                    .styles({
                        "fill-opacity": options.circle_opacity
                    })
                    .call(subs_tip);
            });
        };

        var render_pair_substitutions = function(substitutions, plot, options) {
            var subs_tip = SVGTip({
                content: function(d) {
                    var site = (options.n_sites === 1 ? 1 : d.site);
                    var description = ["Pair Substitution:"];
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

            var subs = plot.nodes.selectAll(".substitution")
                .data(substitutions)
                .enter()
                .append("g")
                .filter(function(d) {
                    var node = d3.select(this.parentNode).datum();
                    return node.index === d.node;
                })
                .filter(function(d) {
                    if (options.sites === "all") { return true; }
                    if (options.sites === "none") { return false; }
                    return options.sites.indexOf(d.site) !== -1;
                })
                .attrs({
                    "class": "substitution",
                    "transform": function(d) {
                        var node = d3.select(this.parentNode).datum();
                        return translate(plot.scale(d.time - node.length), 0);
                    }
                })
                .styles({
                    "fill-opacity": options.circle_opacity

                })
                .call(subs_tip);

            subs.append("circle")
            .attrs({
                "r": options.circle_size,
                "cy": -options.circle_size,
                "fill": function(d) {
                    if (d.hasOwnProperty("first_color")) return(d.first_color);
                    else return("#FC571F");
                }
            });

            subs.append("circle")
            .attrs({
                "r": options.circle_size,
                "cy": options.circle_size,
                "fill": function(d) {
                    if (d.hasOwnProperty("first_color")) return(d.first_color);
                    else return("#FC571F");
                }
            });
        };

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

        var render = function(data, plot, options, width, height) {
            var root_node = render_tree(data.tree, plot.options, width, height);

            if (!is_empty_object(data.intervals) && plot.options.plot_intervals) {
                var intervals = render_intervals(root_node, data.intervals, plot, plot.options);
            }
            if (!is_empty_object(data.substitutions)) {
                if(plot.options.sim_type == "codon_pair") {
                    var subs = render_pair_substitutions(data.substitutions, plot, plot.options);
                } else {
                    var subs = render_substitutions(root_node, data.substitutions, plot, plot.options);
                }

            }
        };

        return {
            renderValue: function(object) {
                data.tree = object.tree;
                data.substitutions = object.substitutions;
                data.intervals = object.intervals;
                //assign_parents(data.tree);

                verify_options(object.options);
                plot.options = object.options;

                render(data, plot, plot.options, width, height);
            },
            resize: function(width, height) {
                render(data, plot, plot.options, width, height);
            },
            plot: plot
        };
    }
});
