var substitution_tooltip = function(plot) {
    return SVGTip({
        content: function(d) {
            var site = (plot.options.n_sites === 1 ? 1 : d.site);
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
};
