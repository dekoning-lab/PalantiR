var tooltip = function(plot, name) {
    return SVGTip({
        content: function(d) {
            // FIX (2026-08-20, MIN11): dropped a `site` variable computed from
            // plot.options.n_sites -- an option the R side has never emitted, so
            // the read was always undefined -- and never used afterwards.
            var description = [name];
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
