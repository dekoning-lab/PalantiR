var monospace_char_width = function(font_size) {
  var dummy = document.getElementById('dummy');
  if (dummy === null) {
    dummy = document.createElement('span');
    dummy.textContent = "a";
    dummy.id = "dummy";
    dummy.style.position = 'absolute';
    dummy.style.top = '-100px';
    dummy.style.fontSize = String(font_size) + 'px';
    dummy.style.fontFamily = 'monospace';
    document.body.appendChild(dummy);
  }
  return dummy.offsetWidth;
};
