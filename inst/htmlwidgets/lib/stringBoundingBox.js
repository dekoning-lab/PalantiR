var string_bounding_box = function(text, font_family, font_size, parent_selection) {
    var t = parent_selection.append("text")
        .attrs({
            "x": -1000,
            "y": -1000,
            "font-size": font_size,
            "font_family": font_family
        })
        .text(text);

    var box = t.node().getBBox();
    t.remove();

    return {width: box.width, height: box.height};
};
