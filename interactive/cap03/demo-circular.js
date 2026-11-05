/* Simulación «Del momento de cada episodio al vector». Datos sintéticos (interactive/cap03/circular.js).
 * No lee nada del proyecto ni de data/. Script clásico, sin red. */
(function () {
  'use strict';
  var C = window.BRCap03Circular, raiz = document.getElementById('demo-circular');
  if (!C || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };
  var DEF = { media: 11, dispersion: 2 };
  var st = Object.assign({}, DEF);

  function num(v, d) { return v.toFixed(d).replace('-', '−').replace('.', ','); }

  function escribir() {
    $('k-media').value = st.media; $('k-media-val').textContent = num(st.media, 2) + ' h tras L5';
    $('k-disp').value = st.dispersion; $('k-disp-val').textContent = num(st.dispersion, 1) + ' h';
  }

  var W = 760;
  function medir() { var w = $('k-plot').getBoundingClientRect().width; W = Math.max(300, Math.round(w || 760)); }

  function dibujar(e, m) {
    var apilado = W < 600, H = apilado ? 640 : 320, s = '';
    var r = apilado ? Math.min(120, W / 2 - 50) : 105, cx = apilado ? W / 2 : 50 + r, cy = apilado ? 160 : 150;
    // Reloj: 0 h (punto medio de L5) arriba, horas en el sentido de las agujas.
    function P(h, rad) { var f = 2 * Math.PI * h / 24; return [cx + rad * Math.sin(f), cy - rad * Math.cos(f)]; }
    s += '<circle cx="' + cx + '" cy="' + cy + '" r="' + r + '" fill="none" class="eje"/>';
    for (var h = 0; h < 24; h += 3) {
      var a = P(h, r), b = P(h, r + 6), t = P(h, r + 20);
      s += '<line class="eje" x1="' + a[0] + '" y1="' + a[1] + '" x2="' + b[0] + '" y2="' + b[1] + '"/>';
      s += '<text x="' + t[0] + '" y="' + (t[1] + 4) + '" text-anchor="middle" font-size="12">' + (h === 0 ? '0 h (L5)' : h + ' h') + '</text>';
    }
    e.forEach(function (p) { var q = P(p.hora, r - 8); s += '<circle class="ep" cx="' + q[0].toFixed(1) + '" cy="' + q[1].toFixed(1) + '" r="' + (1.5 + 1.5 * p.peso).toFixed(1) + '"/>'; });
    var v = P(m.theta, r * m.R);
    s += '<line class="vector" x1="' + cx + '" y1="' + cy + '" x2="' + v[0].toFixed(1) + '" y2="' + v[1].toFixed(1) + '"/><circle class="punto" cx="' + v[0].toFixed(1) + '" cy="' + v[1].toFixed(1) + '" r="5"/>';
    var lin = P(C.mediaLineal(e), r * 0.55);
    s += '<line class="curva-1" x1="' + cx + '" y1="' + cy + '" x2="' + lin[0].toFixed(1) + '" y2="' + lin[1].toFixed(1) + '"/>';
    s += '<text x="' + cx + '" y="' + (cy + r + 44) + '" text-anchor="middle" font-size="12">' + (apilado ? 'Hora de cada episodio respecto al punto medio de L5' : 'Hora tras el punto medio de L5') + '</text>';
    // Plano del modelo: x = R·cos θ, y = R·sin θ.
    var lado = apilado ? Math.min(W - 80, 240) : 210, px = apilado ? (W - lado) / 2 : W - lado - 16, py = apilado ? 370 : 40;
    function X(c) { return px + (c + 1) / 2 * lado; } function Y(c) { return py + (1 - c) / 2 * lado; }
    s += '<rect x="' + px + '" y="' + py + '" width="' + lado + '" height="' + lado + '" fill="none" class="eje"/>';
    s += '<line class="rejilla" x1="' + X(-1) + '" x2="' + X(1) + '" y1="' + Y(0) + '" y2="' + Y(0) + '"/><line class="rejilla" x1="' + X(0) + '" x2="' + X(0) + '" y1="' + Y(-1) + '" y2="' + Y(1) + '"/>';
    s += '<circle class="guia" cx="' + X(0) + '" cy="' + Y(0) + '" r="' + (lado / 2) + '"/>';
    [-1, 0, 1].forEach(function (c) {
      s += '<text x="' + X(c) + '" y="' + (py + lado + 16) + '" text-anchor="middle" font-size="12">' + num(c, 0) + '</text>';
      s += '<text x="' + (px - 6) + '" y="' + (Y(c) + 4) + '" text-anchor="end" font-size="12">' + num(c, 0) + '</text>';
    });
    s += '<text x="' + (px + lado / 2) + '" y="' + (py + lado + 34) + '" text-anchor="middle" font-size="12">R·cos θ</text>';
    s += '<text transform="translate(' + (px - 28) + ' ' + (py + lado / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">R·sin θ</text>';
    s += '<circle class="punto" cx="' + X(m.C).toFixed(1) + '" cy="' + Y(m.S).toFixed(1) + '" r="6"/>';
    $('k-plot').setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    $('k-plot').innerHTML = s;
  }

  function lectura(e, m) {
    var lin = C.mediaLineal(e), err = C.distanciaHoras(lin, m.theta), p = [];
    p.push('<p><strong>Media circular:</strong> θ = ' + num(m.theta, 2) + ' h tras L5, con R = ' + num(m.R, 2) + ' (1 = todos los episodios a la misma hora; 0 = repartidos por todo el día).</p>');
    p.push('<p><strong>Coordenadas para el modelo:</strong> R·cos θ = ' + num(m.C, 2) + ' y R·sin θ = ' + num(m.S, 2) + '. ' + (m.R < 0.2 ? 'Con R tan bajo, el punto queda cerca del origen: el niño apenas tiene hora preferida y el modelo lo trata casi como «sin <em>timing</em>», sin inventarle un θ.' : 'Es el punto del plano que entra en el modelo.') + '</p>');
    p.push('<p><strong>Media lineal de las horas:</strong> ' + num(lin, 2) + ' h tras L5' + (err > 1 ? ', a <strong>' + num(err, 1) + ' h</strong> de la circular en el reloj: promediar horas como números falla cuando los episodios rodean las 0 h.' : ', parecida a la circular: aquí los episodios no rodean las 0 h.') + '</p>');
    $('k-lectura').innerHTML = p.join('');
  }

  function render() { medir(); escribir(); var e = C.episodios(st.media, st.dispersion), m = C.mediaCircular(e); dibujar(e, m); lectura(e, m); }

  $('k-media').addEventListener('input', function () { st.media = Number($('k-media').value); render(); });
  $('k-disp').addEventListener('input', function () { st.dispersion = Number($('k-disp').value); render(); });
  $('k-tarde').addEventListener('click', function () { st = { media: 11, dispersion: 2 }; render(); });
  $('k-dispersa').addEventListener('click', function () { st = { media: 11, dispersion: 8 }; render(); });
  $('k-cero').addEventListener('click', function () { st = { media: 0, dispersion: 1.5 }; render(); });
  var tm = null; window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  render();
})();
