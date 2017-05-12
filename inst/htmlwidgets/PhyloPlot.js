HTMLWidgets.widget({

  name: 'PhyloPlot',

  type: 'output',

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

    var ifelse = function(expr, yes, no) {
      if (expr) {
        return yes;
      } else {
        return no;
      }
    };

    var max_attr = function(array, attr) {
      var attributes = array.map(function(n) { return n[attr]; });
      return d3.max(attributes);
    };

    var assign_parents = function(node, parent) {
      if (parent !== undefined) {
        node.parent = parent;
      }
      if (node.children) {
        for (var i = node.children.length - 1; i >= 0; i--) {
          assign_parents(node.children[i], node);
        }
      }
    };

    var traverse = function(node, callback) {
      callback(node);
      if (node.children) {
        for (var i = node.children.length - 1; i >= 0; i--) {
          traverse(node.children[i], callback);
        }
      }
    };

    var traversal = function(node, filter) {
      var descendants = [node];
      function walk(n) {
        if (n.children) {
           for (var i = n.children.length - 1; i >= 0; i--) {
            descendants.push(n.children[i]);
            walk(n.children[i]);
          }
        }
      }
      walk(node);
      return descendants;
    };

    var is_leaf = function(node) { return !node.hasOwnProperty('children'); };

    var right_angle_diagonal = function() {
      var projection = function(d) { return [d.y, d.x]; };

      var path = function(path_data) {
        return 'M' + path_data[0] + ' ' + path_data[1] + ' ' + path_data[2];
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

    var translate = function(x, y) { return 'translate(' + x + ',' + y + ')'; };

    var swap = function(pair) {return [pair[1], pair[0]]};

    var verify_options = function(options) {
      var get_default = function (x, d) {
        return x === undefined ? d : x;
      };
      options.sites = get_default(options.sites, 'all');

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

    var download_svg = function(element_id, file_name) {
    	var svg = document.getElementById(element_id);
    	var serializer = new XMLSerializer();
    	var source = serializer.serializeToString(svg);
    	source = '<?xml version="1.1" standalone="no"?>\r\n' + source;
    	var pom = document.createElement("a");
    	pom.setAttribute("href","data:image/svg+xml;charset=utf-8," + encodeURIComponent(source));
    	pom.setAttribute("download", file_name);
    	if (document.createEvent) {
    		var event = document.createEvent("MouseEvents");
    		event.initEvent("click",true,true);
    		pom.dispatchEvent(event);
    	} else {
    		pom.click();
    	}
    };

    var render_tree = function(root, options, width, height) {
      var layout = d3.layout.cluster()
          .size([height - (2 * options.padding), width - (2 * options.padding)]);

      var nodes = layout(root);
      var links = layout.links(nodes);
      var leaves = traversal(root).filter(is_leaf);

      var all_nodes = traversal(root);

      traverse(root, function(node) {
        // Assign distance from root
        node.dist = (node.parent ? node.parent.dist : 0) + (node.length || 0);
        // Assign visible label size
        node.text_size = (node.label ? (node.label.length * monospace_char_width(options.label_font_size)) : 0);
      });

      var max_dist = max_attr(nodes, 'dist');
      var max_text_size = options.show_labels ? max_attr(nodes, 'text_size') : 0;

      plot.scale = d3.scale.linear()
                  .domain([ 0, max_dist ])
                  .range([ 0, width - (2 * options.padding) - max_text_size ]);

      var inverse_scale = d3.scale.linear()
                              .domain(plot.scale.domain())
                              .range(swap(plot.scale.range()));

      var axis = d3.svg.axis()
                    .scale(inverse_scale)
                    .orient('bottom')
                    .ticks(options.axis_ticks)
                    .tickSize(height - (3 * options.padding))
                    .tickFormat(d3.format('2r'));

      var h = options.padding;
      var spacing = (height - (options.padding * 4)) / leaves.length;
      for (var i = leaves.length - 1; i >= 0; i--) {
        leaves[i].x = h;
        h += spacing;
      }

      for (var k = all_nodes.length - 1; k >= 0; k--) {
        var n = all_nodes[k];
        n.y = plot.scale(n.dist);
        if (n.children) {
          var a = n.children[0];
          var b = n.children[n.children.length-1];
          n.x = (n.children[0].x + n.children[n.children.length-1].x ) / 2;
        }
      }

      var diagonal = right_angle_diagonal();

      document.onkeydown = function (event) {
        if (event.key === 's') download_svg("phylogram", "PalantiR_phylogeny.svg");
      };

      d3.select("#" + el.id).selectAll("svg").remove();

      plot.svg = d3.select('#' + el.id)
          .append('svg')
          .attr('width', width)
          .attr('height', height)
          .attr('id', 'phylogram')
          .append('g')
          .attr('transform', translate(options.padding, options.padding));

      plot.vis = plot.svg.append('g');

      if (options.show_axis) {
        plot.vis.append('g')
            .attr('class', 'axis')
            .call(axis)
            .selectAll('text')
            .style({
              'font-family': 'monospace',
              'font-size': options.axis_tick_font_size,
              'fill': '#ddd'
            });
        plot.vis.selectAll('.axis')
            .style({
              'fill': 'none',
              'stroke': '#ddd',
              'shape-rendering': 'crispEdges'
            });
      }

      plot.links = plot.vis.selectAll('.link')
          .data(links)
          .enter()
          .append('path')
          .attr('class', 'link')
          .attr({
            'stroke-width': options.line_width,
            'stroke': '#ccc',
            'fill': 'none'
          })
          .attr('d', diagonal)
          .call(SVGTip({
            content: function(d) {
              return ["Branch:",
                      "index: " + d.target.index,
                      "length: " + d.target.length];
            },
            parent: plot.svg
          }));

      plot.nodes = plot.vis.selectAll('.node')
          .data(nodes)
          .enter()
          .append('g')
          .attr('class', 'node')
          .attr('transform', function(d) { return translate(d.y, d.x); })
          .attr({
            'font-family': 'monospace'
          });

      if (options.show_labels) {
        plot.nodes.append('text')
            .style('font-size', options.label_font_size)
            .attr('dx', function(d) { return d.children ? -8 : 8; })
            .attr('dy', 3)
            .text(function(d) { return d.children ? null : d.label; });
      }
    };

    var render_substitutions = function(substitutions, plot, options) {
      var subs_tip = SVGTip({
        content: function(d) {
          var site = ifelse(options.n_sites === 1, 1, d.site);
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

      var subs = plot.nodes.selectAll('.substitution')
          .data(substitutions)
          .enter()
          .append('circle')
          .filter(function(d) {
            var node = d3.select(this.parentNode).datum();
            return node.index === d.node;
          })
          .filter(function(d) {
            if (options.sites === 'all') { return true; }
            if (options.sites === 'none') { return false; }
            return options.sites.indexOf(d.site) !== -1;
          })
          .attr('class', 'substitution')
          .attr('r', options.circle_size)
          .attr('cx', function(d) {
            var node = d3.select(this.parentNode).datum();
            return plot.scale(d.time - node.length);
          })
          .attr('fill', function(d) {
            if (d.hasOwnProperty("color")) return(d.color);
            else return('#FC571F');
          })
          .style({'fill-opacity': options.circle_opacity})
          .call(subs_tip);
    };

    var render_pair_substitutions = function(substitutions, plot, options) {
      var subs_tip = SVGTip({
        content: function(d) {
          var site = ifelse(options.n_sites === 1, 1, d.site);
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

      var subs = plot.nodes.selectAll('.substitution')
          .data(substitutions)
          .enter()
          .append('g')
          .filter(function(d) {
            var node = d3.select(this.parentNode).datum();
            return node.index === d.node;
          })
          .filter(function(d) {
            if (options.sites === 'all') { return true; }
            if (options.sites === 'none') { return false; }
            return options.sites.indexOf(d.site) !== -1;
          })
          .attr('class', 'substitution')
          .attr('transform', function(d) {
              var node = d3.select(this.parentNode).datum();
              return translate(plot.scale(d.time - node.length), 0);
          })
          .style({'fill-opacity': options.circle_opacity})
          .call(subs_tip);

      subs.append('circle')
          .attr('r', options.circle_size)
          .attr('cy', function(d) {
              return -options.circle_size;
          })
          .attr('fill', function(d) {
            if (d.hasOwnProperty("first_color")) return(d.first_color);
            else return('#FC571F');
          });

      subs.append('circle')
          .attr('r', options.circle_size)
          .attr('cy', function(d) {
              return options.circle_size;
          })
          .attr('fill', function(d) {
            if (d.hasOwnProperty("second_color")) return(d.second_color);
            else return('#FC571F');
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
      var color_scale = d3.scale.category10();
      var intervals = plot.nodes.selectAll('.interval')
          .data(intervals)
          .enter()
          .append('line')
          .filter(function(d){
            var node = d3.select(this.parentNode).datum();
            return node.index === d.node;
          })
          // Should we filter on sites?
          .attr({
            'class': 'interval',
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
          .style({
            'stroke': function(d) {
              return color_scale(d.state);
            },
            'stroke-width': options.line_width,
            'stroke-opacity': options.interval_opacity
          })
          .call(interval_tip);
    };

    var render = function(data, plot, options, width, height) {
      render_tree(data.tree, plot.options, width, height);

      if (!is_empty_object(data.intervals) && plot.options.plot_intervals) {
        var intervals = render_intervals(data.intervals, plot, plot.options);
      }
      if (!is_empty_object(data.substitutions)) {
          if(plot.options.sim_type == "codon_pair") {
              var subs = render_pair_substitutions(data.substitutions, plot, plot.options);
          } else {
              var subs = render_substitutions(data.substitutions, plot, plot.options);
          }

      }
    };

    return {
      renderValue: function(object) {
        data.tree = object.tree;
        data.substitutions = object.substitutions;
        data.intervals = object.intervals;
        assign_parents(data.tree);

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
