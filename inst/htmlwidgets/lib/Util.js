var is_empty_object = function(object) {
    return Object.keys(object).length === 0 && object.constructor === Object;
};

var attr_getter = function(array, attr) {
    return array.map(function(n) { return n.data[attr]; });
};

var swap = function(pair) {return [pair[1], pair[0]]};
