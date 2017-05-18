HTMLWidgets.widget({

    name: "PhyloPlot",

    type: "output",

    factory: function(el, width, height) {

        var plot = {
            container: null,    // Element containing the figure
            hierarchy: null,    // Tree node hierarchy
            scale: null,        // d3 scale object
            options: null,      // Optional flags
            svg: null,          // SVG element container (include tooltips)
            vis: null           // SVG group for visible plot elements
        };

        var data;

        var render = function(data, width, height) {
            plot.container = el;
            plot.options = data.options;
            plot.options.width = width;
            plot.options.height = height;

            // promote to array
            if (plot.options.sites.length === undefined) {
                plot.options.sites = [plot.options.sites];
            }

            render_tree(plot, data.tree);

            if (!is_empty_object(data.intervals) && plot.options.plot_intervals && data.intervals) {
                render_intervals(plot, data.intervals);
            }
            if (!is_empty_object(data.substitutions)) {
                if(plot.options.sim_type == "codon_pair") {
                    render_pair_substitutions(plot, data.substitutions);
                } else {
                    render_substitutions(plot, data.substitutions);
                }

            }
        };

        return {
            renderValue: function(object) {
                data = object;
                render(data, width, height);
            },
            resize: function(width, height) {
                render(data, width, height);
            },
            plot: plot
        };
    }
});
