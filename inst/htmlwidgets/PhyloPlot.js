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

            if (data.intervals && !is_empty_object(data.intervals) && plot.options.plot_intervals) {
                render_intervals(plot, data.intervals);
            }
            if (!is_empty_object(data.substitutions)) {
                if(plot.options.sim_type == "codon_pair") {
                    render_pair_substitutions(plot, data.substitutions);
                } else {
                    render_substitutions(plot, data.substitutions);
                }
            }

            // FIX (2026-08-20, MIN12): a click anywhere on the plot used to
            // download the whole svg whenever the widget was not in the RStudio
            // viewer pane -- undocumented, unguarded, and surprising in a
            // browser or a Shiny app where a click means something else. It is
            // now opt-in via the download_on_click option (PhyloPlot(),
            // default FALSE); the viewer-pane exclusion is kept because the
            // download does not work there.
            if(plot.options.download_on_click &&
               window.location.href.indexOf("viewer_pane=1") <= -1) {
                plot.svg.on("click", function() {
                    download_svg("#" + el.id + ">svg.phylogram", "PalantiR_phylogeny.svg");
                });
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
