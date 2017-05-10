HTMLWidgets.widget({

  name: 'AlignmentPlot',

  type: 'output',

  factory: function(el, width, height) {

    var data = {
      alignment: null,
      taxa: null,
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

    var render = function(data, plot, width, height) {
        d3.select("#" + el.id).selectAll("div").remove();

        plot.container = d3.select("#" + el.id)
            .append("div")
            .style({
                "display": "flex"
            });

        plot.taxa = plot.container
            .append("div")
            .style({
                "overflow-y": "scroll",
                "flex-basis": "content",
                "text-align": "right",
                "max-widht": width + "px",
                "max-height": height + "px"});

        plot.taxa.table = plot.taxa
            .append('table')
            .style({
                "border-collapse": "collapse"
            })
            .append('tbody');

        plot.taxa.rows = plot.taxa.table.selectAll("tr")
            .data(data)
            .enter()
            .append("tr");

        plot.taxa.cells = plot.taxa.rows.selectAll("td")
            .data(function(d) {
                return d.taxon;
            })
            .enter()
            .append("td")
            .style({
                "border": "1px solid white",
            })
            .text(function(d, i) {
                return d;
            });

        plot.alignment = plot.container
            .append("div")
            .style({
                "text-align": "center",
                "overflow-x": "scroll",
                "overflow-y": "scroll",
                "flex": 1,
                "max-widht": width + "px",
                "max-height": height + "px"});

        plot.alignment.table = plot.alignment
            .append('table')
            .style({
                "border-collapse": "collapse"
            })
            .append('tbody');

        plot.alignment.rows = plot.alignment.table.selectAll("tr")
            .data(data)
            .enter()
            .append("tr");

        plot.alignment.cells = plot.alignment.rows.selectAll("td")
            .data(function(d) {
                return d.sequence;
            })
            .enter()
            .append("td")
            .style({
                "border": "1px solid black"
            })
            .text(function(d, i) {
                return d;
            });

        // Sync scrolling

        var taxa_div  = plot.taxa.node();
        var syncing_taxa_scroll = false;
        var alignment_div  = plot.alignment.node();
        var syncing_alignment_scroll = false;

        taxa_div.onscroll = function() {
            if(!syncing_taxa_scroll) {
                syncing_alignment_scroll = true;
                alignment_div.scrollTop = this.scrollTop;
            }
            syncing_taxa_scroll = false;
        };

        alignment_div.onscroll = function() {
            if(!syncing_alignment_scroll) {
                syncing_taxa_scroll = true;
                taxa_div.scrollTop = this.scrollTop;
            }
            syncing_alignment_scroll = false;
        };
    };

    return {
      renderValue: function(object) {
        render(object, plot, width, height);
      },
      resize: function(width, height) {
        render(data, plot, width, height);
      },
      plot: plot
    };
  }
});
