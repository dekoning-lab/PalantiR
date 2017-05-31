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
