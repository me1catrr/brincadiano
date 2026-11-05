/* Simulación «Del segundo al minuto». Datos sintéticos (interactive/cap02/epocas.js): una hora generada con semilla fija.
 * No lee nada del proyecto ni de data/. Script clásico, sin red. */
(function () {
  'use strict';
  var E = window.BRCap02Epocas, raiz = document.getElementById('demo-epocas');
  if (!E || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };
  var hora = E.generarHora(), x = hora.enmo;
  var DEF = { epoca: 60, umbral: E.UMBRAL_VIGOROSA };
  var st = Object.assign({}, DEF);
  var U_MIN = 100, U_MAX = 900, Y_MAX = 1400;

  function num(v, d) { return v.toFixed(d).replace('-', '−').replace('.', ','); }
  function mediana(a) { var s = a.slice().sort(function (p, q) { return p - q; }), n = s.length; return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2; }

  function escribir() {
    st.umbral = Math.min(U_MAX, Math.max(U_MIN, st.umbral));
    $('e-umbral').value = Math.round(st.umbral);
    $('e-umbral-val').textContent = num(st.umbral, st.umbral % 1 ? 1 : 0) + ' mg';
    Array.prototype.forEach.call(raiz.querySelectorAll('input[name="e-epoca"]'), function (r) { r.checked = Number(r.value) === st.epoca; });
  }

  var W = 760, Lm = 58, Rm = 14;
  function medir() { var w = $('e-plot').getBoundingClientRect().width; W = Math.max(300, Math.round(w || 760)); Lm = W < 480 ? 44 : 58; Rm = W < 480 ? 8 : 14; }

  function dibujar() {
    var H = 300, B = 38, T = 22, h = H - B - T, s = '', i;
    function X(seg) { return Lm + seg / 3600 * (W - Lm - Rm); }
    function Y(mg) { return T + (1 - Math.min(mg, Y_MAX) / Y_MAX) * h; }
    [0, 200, 400, 600, 800, 1000, 1200, 1400].forEach(function (v) {
      s += '<line class="rejilla" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y(v) + '" y2="' + Y(v) + '"/>';
      if (W >= 480 || v % 400 === 0) s += '<text x="' + (Lm - 6) + '" y="' + (Y(v) + 4) + '" text-anchor="end" font-size="12">' + v + '</text>';
    });
    for (var m = 0; m <= 60; m += 10) s += '<text x="' + X(m * 60) + '" y="' + (H - B + 16) + '" text-anchor="middle" font-size="12">' + m + '</text>';
    s += '<text x="' + ((Lm + W - Rm) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" font-size="12">Minuto de la hora sintética</text>';
    s += '<text transform="translate(13 ' + (T + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">ENMO (mg)</text>';
    E.TRAMOS.forEach(function (tr) {
      if (X(tr[1]) - X(tr[0]) < 44) return;   // tramo demasiado estrecho para su rótulo (carrera; en móvil, los de 5 min)
      s += '<text x="' + ((X(tr[0]) + X(tr[1])) / 2) + '" y="' + (T - 7) + '" text-anchor="middle" font-size="11" class="tramo">' + tr[2] + '</text>';
    });
    var med = E.agregar(x, st.epoca);
    for (i = 0; i < med.length; i++) if (med[i] >= st.umbral) s += '<rect class="marca" x="' + X(i * st.epoca) + '" y="' + T + '" width="' + Math.max(1, X((i + 1) * st.epoca) - X(i * st.epoca)) + '" height="' + h + '"/>';
    var d = '';
    for (i = 0; i < x.length; i++) d += (i ? 'L' : 'M') + X(i + 0.5).toFixed(1) + ',' + Y(x[i]).toFixed(1);
    s += '<path class="curva-0" d="' + d + '"/>';
    if (st.epoca > 1) {
      d = '';
      for (i = 0; i < med.length; i++) d += (i ? 'L' : 'M') + X(i * st.epoca).toFixed(1) + ',' + Y(med[i]).toFixed(1) + 'L' + X((i + 1) * st.epoca).toFixed(1) + ',' + Y(med[i]).toFixed(1);
      s += '<path class="curva-2" d="' + d + '"/>';
    }
    s += '<line class="curva-1" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y(st.umbral) + '" y2="' + Y(st.umbral) + '"/>';
    s += '<rect x="' + Lm + '" y="' + T + '" width="' + (W - Lm - Rm) + '" height="' + h + '" fill="none" class="eje"/>';
    $('e-plot').setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    $('e-plot').innerHTML = s;
  }

  function lectura() {
    var u = st.umbral, min = E.minutosSobre(x, st.epoca, u), p = [];
    var todos = E.EPOCAS.map(function (e) { return e + ' s → ' + num(E.minutosSobre(x, e, u), 2) + ' min'; }).join('; ');
    p.push('<p><strong>' + num(min, 2) + ' min</strong> de esta hora tienen un <em>epoch</em> de ' + st.epoca + ' s cuya media llega a ' + num(u, u % 1 ? 1 : 0) + ' mg.</p>');
    p.push('<p>Con el mismo umbral y los mismos movimientos: ' + todos + '.</p>');
    var dur = hora.rafagas.map(function (r) { return r.duracion; }), cortas = dur.filter(function (v) { return v < 15; }).length;
    p.push('<p>La hora tiene ' + dur.length + ' ráfagas sintéticas (mediana ' + num(mediana(dur), 0) + ' s; ' + cortas + ' duran menos de 15 s) y una carrera de 2 min. La media más alta con <em>epochs</em> de ' + st.epoca + ' s es ' + num(E.maximo(E.agregar(x, st.epoca)), 0) + ' mg; el segundo más alto, ' + num(E.maximo(x), 0) + ' mg.</p>');
    $('e-lectura').innerHTML = p.join('');
  }

  function render() { medir(); escribir(); dibujar(); lectura(); }

  Array.prototype.forEach.call(raiz.querySelectorAll('input[name="e-epoca"]'), function (r) {
    r.addEventListener('change', function () { st.epoca = Number(r.value); render(); });
  });
  $('e-umbral').addEventListener('input', function () { st.umbral = Number($('e-umbral').value); render(); });
  $('e-vig').addEventListener('click', function () { st.umbral = E.UMBRAL_VIGOROSA; render(); });
  $('e-mod').addEventListener('click', function () { st.umbral = E.UMBRAL_MODERADA; render(); });
  $('e-reset').addEventListener('click', function () { st = Object.assign({}, DEF); render(); });
  var tm = null; window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  render();
})();
