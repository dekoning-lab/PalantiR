HTMLWidgets.widget({

    name: 'AlignmentPlot',

    type: 'output',

    factory: function(el, width, height) {

        var data = {};
        var plot = {};

        var render = function(data, plot, width, height) {
            data = data;

            var sites = Array.apply(null, Array(data[0].sequence.length)).map(function (_, i) { return i; });

            var table_styles = {
                "font-family": "Courier, monospace",
                "text-align": "center",
                "border-collapse": "collapse"
            };

            plot.container = d3.select("#" + el.id)
                .append("div")
                .styles({
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
                .data(data)
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

            var cell_width = cells.node().offsetWidth;
            var cell_height = cells.node().offsetHeight;
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
                .data(data)
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
                    "padding-bottom": "15px",
                    "transform": "rotate(-90deg)  translate(0px, 5px)",
                })
                .text(function(d) {
                    return d;
                });

            // reposition
            var taxa_width = plot.taxa.table.node().offsetWidth;
            var site_height = plot.site.table.node().offsetHeight;

            plot.alignment.div
                .styles({
                    "left": taxa_width + "px",
                    "top": site_height + "px",
                    "max-height": height - site_height + "px",
                    "max-width": width - taxa_width + "px",
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

            // Sync scrolling
            var syncing_taxa_scroll = false;
            var syncing_alignment_scroll = false;
            var syncing_site_scroll = false;

            //plot.taxa.div.on("scroll", function() {
            //    if(!syncing_taxa_scroll) {
            //        syncing_alignment_scroll = true;
            //        plot.alignment.div.node().scrollTop = this.scrollTop;
            //    }
            //    syncing_taxa_scroll = false;
            //});

            //plot.site.div.on("scroll", function() {
            //    if(!syncing_site_scroll) {
            //        syncing_alignment_scroll = true;
            //        plot.alignment.div.node().scrollLeft = this.scrollLeft;
            //    }
            //    syncing_site_scroll = false;
            //});

            plot.alignment.div.on("scroll", function() {
                if(!syncing_alignment_scroll) {
                    syncing_taxa_scroll = true;
                    plot.taxa.div.node().scrollTop = this.scrollTop;
                    plot.site.div.node().scrollLeft = this.scrollLeft;
                }
                syncing_alignment_scroll = false;
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
                plot.container
                    .styles({
                        "width": width + "px",
                        "height": height + "px"
                    });
                plot.alignment.div
                    .styles({
                        "max-height": height - site_height + "px",
                        "max-width": width - taxa_width + "px",
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
