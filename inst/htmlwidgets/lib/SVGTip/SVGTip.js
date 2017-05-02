var SVGTip = function (options) {
  var getDefault = function (x, d) {
    return x === undefined ? d : x;
  };

  var isArray = function (value) {
    return value &&
      typeof value === 'object' &&
      typeof value.length === 'number' &&
      typeof value.splice === 'function' &&
      !(value.propertyIsEnumerable('length'));
  };

  var tip = function (x, y) {
    return [
      [x, y],
      [x + 3, y - 4],
      [x - 3, y - 4]
    ].map(function (x) {
      return x.join(',');
    }).join(' ');
  };

  var translate = function (coord) {
    return 'translate(' + coord[0] + ',' + coord[1] + ')';
  };

  var content = getDefault(options.content, '');
  var style = getDefault(options.style, '');
  var parent = getDefault(options.parent, d3.select('svg'));
  var padding = 3; //px
  var cornerRadius = 5; //px
  var lineHeight = 1.2; //em

  var tooltip = parent.append('g')
    .attr('class', 'svg-tooltip')
    .style({
      'text-align': 'center',
      'font': '12px sans-serif',
      'fill': 'black',
      'opacity': '0.8',
      'pointer-events': 'none'
    });
    var tooltipStyle = {
      'fill': 'black',
      'opacity': '0.8',
      'corner-radius': '2'
    };

  tooltip.container = tooltip.append('g');
  tooltip.tip = tooltip.container.append('polygon')
    .attr('class', 'tip')
    .style(tooltipStyle);
  tooltip.box = tooltip.container.append('rect')
    .attr({
      class: 'box',
      rx: cornerRadius,
      ry: cornerRadius
    })
    .style(tooltipStyle);

  tooltip.text = tooltip.append('text')
    .attr('class', 'text')
    .style({
      'fill': 'white',
      'stroke': 'none',
      'text-anchor': 'middle',
      'alignment-baseline': 'middle',
      'text-anchor': 'middle'
    });

  var show = function () {
    tooltip.style('display', 'inline');
  };

  var hide = function () {
    tooltip.style('display', 'none');
  };

  return function (selection) {
    selection.on('mouseover', show)
      .on('mouseout', hide)
      .on('mousemove', function (d, n) {
        var point = d3.mouse(parent.node());
        point.x = point[0];
        point.y = point[1];

        tooltip.text.selectAll('tspan').remove();

        var lines = ((typeof content === 'function') ? content(d, n) : content);
        if (!isArray(lines)) lines = [lines];
        var s = ((typeof style === 'function') ? style(d, n) : style);

        lines.forEach(function (line) {
          tooltip.text.append('tspan').text(line).attr({
            x: 0,
            dy: lineHeight + 'em'
          })
          .style({
            'text-anchor': 'middle',
            'alignment-baseline': 'middle'
          });
          tooltip.text.style(s);
        });

        var bbox = tooltip.text.node().getBBox();

        tooltip.attr('transform',
          translate([
            point.x,
            point.y - (bbox.height + bbox.y + padding + 4)
          ]));

        tooltip.tip.attr({
          points: tip(
            bbox.x + (bbox.width / 2),
            bbox.y + bbox.height + padding + 4)
        });
        tooltip.box.attr({
          x: bbox.x - padding,
          y: bbox.y - padding,
          width: bbox.width + (padding * 2),
          height: bbox.height + (padding * 2)
        });
      });
  };
};
