//http://html5-demos.appspot.com/static/a.download.html
function downlaod_alignment(element) {
	window.URL = window.webkitURL || window.URL;
	var output = document.querySelector('alignmentOutput');
	var prevLink = output.querySelector('a');
	if (prevLink) {
		window.URL.revokeObjectURL(prevLink.href);
		output.innerHTML = '';
	}
	var svg_export = al.exportSvg({
		whitespace: true,
		excludeAttributes: ["title"]
	});
	var bb = new Blob([svg_export], {type: "image/svg+xml"});

	var a = document.createElement('a');
	a.download = "alignment.svg";
	a.href = window.URL.createObjectURL(bb);
	a.textContent = 'Download ready';

	a.dataset.downloadurl = ["image/svg+xml", a.download, a.href].join(':');
	a.draggable = true; // Don't really need, but good practice.
	a.classList.add('dragout');

	output.appendChild(a);

	a.onclick = function(e) {
		if ('disabled' in this.dataset) {
			return false;
		}

		this.textContent = 'Downloaded';

		// Need a small delay for the revokeObjectURL to work properly.
		setTimeout(function() {
			window.URL.revokeObjectURL(this.href);
		}, 1500);
	};
}
