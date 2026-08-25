HTMLWidgets.widget({

    name: 'AlignmentPlot',

    type: 'output',

    factory: function(el, width, height) {

        var data = {};
        var plot = {};

        var render = function(data, plot, width, height) {
            // Historical widgets receive the sequence array directly. Site-
            // class-aware widgets wrap that array together with one metadata
            // record per alignment column.
            var sequences = Array.isArray(data) ? data : data.sequences;
            var site_classes = Array.isArray(data) ? null : data.site_classes;
            var has_site_classes = site_classes &&
                site_classes.length === sequences[0].sequence.length;

            var sites = Array.apply(null, Array(sequences[0].sequence.length)).map(function (_, i) {
                if (has_site_classes) return site_classes[i];
                return i;
            });

            var table_styles = {
                "font-family": "Courier, monospace",
                "text-align": "center",
                "border-collapse": "collapse"
            };

            // FIX (2026-08-20, MIN10): render() appended a fresh container
            // without removing the previous one, so every Shiny re-render (and
            // every resize, which calls render() again) stacked another copy of
            // the alignment table below the last. RenderTree.js already clears
            // its svg the same way.
            d3.select("#" + el.id).selectAll("div").remove();

            plot.container = d3.select("#" + el.id)
                .append("div")
                .styles({
                    "border": "1px dashed grey",
                    "position": "relative",
                    "overflow": "hidden",
                    "width": width + "px",
                    "height": height + "px"
                });

            // Alignment
            plot.alignment = {};
            plot.alignment.div = plot.container
                .append("div")
                .styles({
                    "position": "absolute",
                    "overflow-x": "scroll",
                    "overflow-y": "scroll"
                });

            plot.alignment.table = plot.alignment.div
                .append("table")
                .styles(table_styles);

            var cells = plot.alignment.table
                .selectAll("tr")
                .data(sequences)
                .enter()
                .append("tr")
                .selectAll("td")
                .data(function(d) { return d.sequence; })
                .enter()
                .append("td")
                .styles({
                    "border": "1px solid black",
                    "padding": "3px",
                    "background-color": function(d) {
                        return d.color;
                    }
                })
                .text(function(d) {
                    return d.state;
                });

            var table_width = plot.alignment.table.node().offsetWidth;
            var table_height = plot.alignment.table.node().offsetHeight;

            // Taxa
            plot.taxa = {};
            plot.taxa.div = plot.container
                .append("div")
                .styles({
                    "position": "absolute",
                    "overflow": "hidden"
                });

            plot.taxa.table = plot.taxa.div
                .append("table")
                .styles({
                    "table-layout": "fixed",
                    "height": table_height + "px"
                })
                .styles(table_styles);

            plot.taxa.table.selectAll("tr")
                .data(sequences)
                .enter()
                .append("tr")
                .append("td")
                .styles({
                    "border": "1px solid black",
                    "font-size": "10px",
                })
                .text(function(d) {
                    return d.taxon;
                });

            // Sites
            plot.site = {};
            plot.site.div = plot.container
                .append("div")
                .styles({
                    "position": "absolute",
                    "overflow": "hidden"
                });

            plot.site.table = plot.site.div
                .append("table")
                .styles(table_styles)
                .styles({
                    "table-layout": "fixed",
                    "width": table_width + "px"
                });

            plot.site.table
                .append("tr")
                .selectAll("td")
                .data(sites)
                .enter()
                .append("td")
                .styles({
                    "border": "1px solid black",
                    "font-size": "10px",
                    "padding-bottom": has_site_classes ? "55px" : "15px",
                    "transform": "rotate(-90deg)  translate(0px, 5px)",
                })
                .attr("title", function(d) {
                    if (!has_site_classes) return null;
                    return "Site " + d.site + "; site class " + d.site_class +
                        " (class index " + d.class_index + ")";
                })
                .text(function(d) {
                    if (!has_site_classes) return d;
                    return d.site + " \u00b7 " + d.site_class;
                });

            // reposition
            var taxa_width = plot.taxa.table.node().offsetWidth;
            var site_height = plot.site.table.node().offsetHeight;

            var scrollbar_width = get_scrollbar_width();

            plot.alignment.div
                .styles({
                    "left": taxa_width + "px",
                    "top": site_height + "px",
                    "max-height": height - site_height + scrollbar_width + "px",
                    "max-width": width - taxa_width  + scrollbar_width + "px",
                });

            plot.site.div
                .styles({
                    "left": taxa_width + "px",
                    "max-width": width - taxa_width + "px",
                });

            plot.taxa.div
                .styles({
                    "top": site_height + "px",
                    "max-height": height - site_height + "px",
                });

            plot.alignment.div.on("scroll", function() {
                plot.taxa.div.node().scrollTop = this.scrollTop;
                plot.site.div.node().scrollLeft = this.scrollLeft;
            });
        };

        return {
            renderValue: function(object) {
                data = object;
                render(object, plot, width, height);
            },
            resize: function(width, height) {
                var taxa_width = plot.taxa.table.node().offsetWidth;
                var site_height = plot.site.table.node().offsetHeight;
                var scrollbar_width = get_scrollbar_width();
                plot.container
                    .styles({
                        "width": width + "px",
                        "height": height + "px"
                    });
                plot.alignment.div
                    .styles({
                        "max-height": height - site_height + scrollbar_width + "px",
                        "max-width": width - taxa_width + scrollbar_width + "px",
                    });
                plot.site.div
                    .styles({
                        "max-width": width - taxa_width + "px",
                    });
                plot.taxa.div
                    .styles({
                        "max-height": height - site_height + "px",
                    });
            }
        };
    }
});
