
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
