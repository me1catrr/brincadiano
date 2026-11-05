/* Simulación «¿Qué frecuencias sobreviven?» (interactive/cap04/frecuencias.js): filtro idealizado y aliasing.
 * No lee nada del proyecto ni de data/. Script clásico, sin red. */
(function () {
  'use strict';
  var F = window.BRCap04Frecuencias, raiz = document.getElementById('demo-frecuencias');
  if (!F || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };
  var st = { f: 2 };
  var F_MIN = 0.1, F_MAX = 25;

  function lmap(v) { return F_MIN * Math.pow(F_MAX / F_MIN, v / 1000); }
  function linv(f) { return Math.round(1000 * Math.log(f / F_MIN) / Math.log(F_MAX / F_MIN)); }
  function num(v, d) { return v.toFixed(d).replace('-', '−').replace('.', ','); }
  function fHz(f) { return num(f, f < 1 ? 2 : 1) + ' Hz'; }

  function escribir() { $('f-f').value = linv(st.f); $('f-f-val').textContent = fHz(st.f); }

  var W = 760, Lm = 50, Rm = 12;
  function medir() { var w = $('f-plot').getBoundingClientRect().width; W = Math.max(300, Math.round(w || 760)); Lm = W < 480 ? 40 : 50; Rm = W < 480 ? 8 : 12; }

  function dibujar() {
    var H = 380, T1 = 14, h1 = 140, T2 = 210, h2 = 130, s = '', i;
    // Respuesta del filtro (eje de frecuencia logarítmico)
    function X(f) { return Lm + Math.log(f / F_MIN) / Math.log(F_MAX / F_MIN) * (W - Lm - Rm); }
    function Y1(g) { return T1 + (1 - g) * h1; }
    s += '<rect class="ventana-dia" x="' + X(15) + '" y="' + T1 + '" width="' + (X(F_MAX) - X(15)) + '" height="' + h1 + '"/>';
    s += '<line class="guia" x1="' + X(15) + '" x2="' + X(15) + '" y1="' + T1 + '" y2="' + (T1 + h1) + '"/><text x="' + (X(15) - 4) + '" y="' + (T1 + 14) + '" text-anchor="end" font-size="11" class="tramo">Nyquist, 15 Hz</text>';
    [0, 0.5, 1].forEach(function (g) { s += '<line class="rejilla" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y1(g) + '" y2="' + Y1(g) + '"/><text x="' + (Lm - 6) + '" y="' + (Y1(g) + 4) + '" text-anchor="end" font-size="12">' + num(g, 1) + '</text>'; });
    [0.1, 0.3, 1, 3, 10, 25].forEach(function (f) { s += '<text x="' + X(f) + '" y="' + (T1 + h1 + 16) + '" text-anchor="middle" font-size="12">' + num(f, f < 1 ? 1 : 0) + '</text>'; });
    s += '<text x="' + ((Lm + W - Rm) / 2) + '" y="' + (T1 + h1 + 32) + '" text-anchor="middle" font-size="12">Frecuencia (Hz)</text>';
    s += '<text transform="translate(12 ' + (T1 + h1 / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Ganancia</text>';
    var d = '';
    for (i = 0; i <= 300; i++) { var f = F_MIN * Math.pow(F_MAX / F_MIN, i / 300); d += (i ? 'L' : 'M') + X(f).toFixed(1) + ',' + Y1(F.ganancia(f)).toFixed(1); }
    s += '<path class="curva-2" d="' + d + '"/>';
    s += '<line class="guia" x1="' + X(st.f) + '" x2="' + X(st.f) + '" y1="' + T1 + '" y2="' + (T1 + h1) + '"/><circle class="punto" cx="' + X(st.f) + '" cy="' + Y1(F.ganancia(st.f)) + '" r="6"/>';
    s += '<rect x="' + Lm + '" y="' + T1 + '" width="' + (W - Lm - Rm) + '" height="' + h1 + '" fill="none" class="eje"/>';
    // Señal en el tiempo: original, tras el filtro (solo amplitud) y muestras a 30 Hz
    var dur = st.f < 1 ? 2 / st.f : 1, fa = F.aparente(st.f), g = F.ganancia(st.f);
    function X2(t) { return Lm + t / dur * (W - Lm - Rm); }
    function Y2(v) { return T2 + (1 - (v + 1.1) / 2.2) * h2; }
    // Las muestras de sin(2π·f·t) a 30 Hz coinciden con las de sin(2π·δ·t), con δ = f − 30·redondeo(f/30).
    function muestraAlias(t) { var k = Math.round(st.f / F.FS), sg = st.f - F.FS * k >= 0 ? 1 : -1; return sg * Math.sin(2 * Math.PI * fa * t); }
    var orig = '', filt = '', ali = '', n = 600;
    for (i = 0; i <= n; i++) {
      var t = dur * i / n;
      orig += (i ? 'L' : 'M') + X2(t).toFixed(1) + ',' + Y2(Math.sin(2 * Math.PI * st.f * t)).toFixed(1);
      filt += (i ? 'L' : 'M') + X2(t).toFixed(1) + ',' + Y2(g * Math.sin(2 * Math.PI * st.f * t)).toFixed(1);
      if (st.f > 15) ali += (i ? 'L' : 'M') + X2(t).toFixed(1) + ',' + Y2(muestraAlias(t)).toFixed(1);
    }
    s += '<line class="rejilla" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y2(0) + '" y2="' + Y2(0) + '"/>';
    s += '<path class="curva-0" d="' + orig + '"/><path class="curva-2" d="' + filt + '"/>';
    if (st.f > 15) s += '<path class="curva-1" d="' + ali + '"/>';
    var nm = Math.min(Math.floor(dur * F.FS), 400);
    if (dur <= 2) for (i = 0; i <= nm; i++) { var tm = i / F.FS; s += '<circle class="ep" cx="' + X2(tm).toFixed(1) + '" cy="' + Y2(Math.sin(2 * Math.PI * st.f * tm)).toFixed(1) + '" r="2.6"/>'; }
    s += '<rect x="' + Lm + '" y="' + T2 + '" width="' + (W - Lm - Rm) + '" height="' + h2 + '" fill="none" class="eje"/>';
    s += '<text x="' + ((Lm + W - Rm) / 2) + '" y="' + (H - 6) + '" text-anchor="middle" font-size="12">Tiempo (' + num(dur, dur < 10 ? 1 : 0) + ' s)</text>';
    $('f-plot').setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    $('f-plot').innerHTML = s;
  }

  function lectura() {
    var g = F.ganancia(st.f), p = [];
    p.push('<p><strong>Filtro de los <em>counts</em> (idealizado):</strong> a ' + fHz(st.f) + ' deja pasar el ' + num(100 * g, g < 0.1 ? 1 : 0) + ' % de la amplitud (' + num(F.dB(g), 1) + ' dB). A 10 Hz, donde aún llega el movimiento según Fridolfsson et al. (2019), pasa el ' + num(100 * F.ganancia(10), 1) + ' %.</p>');
    if (st.f > 15) p.push('<p><strong>Por encima de Nyquist:</strong> muestreada a 30 Hz sin filtro previo, una oscilación de ' + fHz(st.f) + ' deja las mismas muestras que una de ' + fHz(F.aparente(st.f)) + ' (línea discontinua): es el <em>aliasing</em>.</p>');
    else p.push('<p><strong>Muestreo a 30 Hz:</strong> ' + fHz(st.f) + ' está por debajo de 15 Hz; las muestras (puntos) bastan para reconstruirla.</p>');
    $('f-lectura').innerHTML = p.join('');
  }

  function render() { medir(); escribir(); dibujar(); lectura(); }

  $('f-f').addEventListener('input', function () { st.f = lmap(Number($('f-f').value)); render(); });
  [['f-1', 1], ['f-3', 3], ['f-10', 10], ['f-20', 20]].forEach(function (b) { $(b[0]).addEventListener('click', function () { st.f = b[1]; render(); }); });
  var tm = null; window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  render();
})();
