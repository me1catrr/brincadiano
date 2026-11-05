/* Simulación «ENMO, MAD y un error de calibración». Datos sintéticos (interactive/cap04/metricas.js).
 * No lee nada del proyecto ni de data/. Script clásico, sin red. */
(function () {
  'use strict';
  var M = window.BRCap04Metricas, raiz = document.getElementById('demo-metricas');
  if (!M || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };
  var st = { offset: 0 };

  function num(v, d) { return v.toFixed(d).replace('-', '−').replace('.', ','); }

  function escribir() { $('m-off').value = st.offset; $('m-off-val').textContent = num(st.offset, 0) + ' mg'; }

  var W = 760, Lm = 58, Rm = 12;
  function medir() { var w = $('m-plot').getBoundingClientRect().width; W = Math.max(300, Math.round(w || 760)); Lm = W < 480 ? 44 : 58; Rm = W < 480 ? 8 : 12; }

  function dibujar(a, ep) {
    var H = 360, T1 = 22, h1 = 110, T2 = 190, h2 = 130, s = '', i;
    function X(t) { return Lm + t / 60 * (W - Lm - Rm); }
    function Y1(g) { return T1 + (1 - (Math.min(Math.max(g, 0.4), 1.8) - 0.4) / 1.4) * h1; }
    var maxv = Math.max.apply(null, ep.map(function (e) { return Math.max(e.enmo, e.mad, e.normaMenos1); }));
    var pasoY = maxv <= 250 ? 50 : 100, top = Math.max(200, Math.ceil(maxv / pasoY) * pasoY);
    function Y2(mg) { return T2 + (1 - Math.max(0, mg) / top) * h2; }
    M.TRAMOS.forEach(function (tr) { s += '<text x="' + ((X(tr[0]) + X(tr[1])) / 2) + '" y="' + (T1 - 8) + '" text-anchor="middle" font-size="11" class="tramo">' + tr[2] + '</text>'; });
    [0.5, 1, 1.5].forEach(function (g) { s += '<line class="rejilla" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y1(g) + '" y2="' + Y1(g) + '"/><text x="' + (Lm - 6) + '" y="' + (Y1(g) + 4) + '" text-anchor="end" font-size="12">' + num(g, 1) + ' g</text>'; });
    var d = '', paso = W < 480 ? 3 : 1;
    for (i = 0; i < a.length; i += paso) d += (i ? 'L' : 'M') + X(i / M.FS).toFixed(1) + ',' + Y1(a[i]).toFixed(1);
    s += '<path class="curva-0" d="' + d + '"/><line class="guia" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y1(1) + '" y2="' + Y1(1) + '"/>';
    s += '<rect x="' + Lm + '" y="' + T1 + '" width="' + (W - Lm - Rm) + '" height="' + h1 + '" fill="none" class="eje"/>';
    s += '<text transform="translate(12 ' + (T1 + h1 / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">‖a‖</text>';
    for (var v = 0; v <= top; v += pasoY) s += '<line class="rejilla" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y2(v) + '" y2="' + Y2(v) + '"/><text x="' + (Lm - 6) + '" y="' + (Y2(v) + 4) + '" text-anchor="end" font-size="12">' + num(v, 0) + '</text>';
    var anchoEp = (X(M.EPOCH) - X(0)), b = anchoEp * 0.26;
    ep.forEach(function (e) {
      var x0 = X(e.t0) + anchoEp * 0.1;
      [['enmo', 'b-enmo'], ['normaMenos1', 'b-norma'], ['mad', 'b-mad']].forEach(function (c, k) {
        var val = Math.max(0, e[c[0]]);
        s += '<rect class="' + c[1] + '" x="' + (x0 + k * b).toFixed(1) + '" y="' + Y2(val).toFixed(1) + '" width="' + (b - 1).toFixed(1) + '" height="' + (Y2(0) - Y2(val)).toFixed(1) + '"/>';
      });
    });
    s += '<rect x="' + Lm + '" y="' + T2 + '" width="' + (W - Lm - Rm) + '" height="' + h2 + '" fill="none" class="eje"/>';
    s += '<text transform="translate(12 ' + (T2 + h2 / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">mg por epoch de 5 s</text>';
    for (var t = 0; t <= 60; t += 10) s += '<text x="' + X(t) + '" y="' + (T2 + h2 + 16) + '" text-anchor="middle" font-size="12">' + t + '</text>';
    s += '<text x="' + ((Lm + W - Rm) / 2) + '" y="' + (H - 6) + '" text-anchor="middle" font-size="12">Tiempo (s)</text>';
    $('m-plot').setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    $('m-plot').innerHTML = s;
  }

  function lectura(ep) {
    var f = function (c, k) { return num(M.mediaTramo(ep, c, M.TRAMOS[k][0], M.TRAMOS[k][1]), 1); }, p = [];
    p.push('<p><strong>De pie y quieto:</strong> ENMO ' + f('enmo', 0) + ' mg; MAD ' + f('mad', 0) + ' mg. <strong>Tumbado:</strong> ENMO ' + f('enmo', 2) + ' mg; MAD ' + f('mad', 2) + ' mg.</p>');
    p.push('<p><strong>Andando:</strong> ENMO ' + f('enmo', 1) + ' mg, pero la norma media menos 1 solo ' + f('normaMenos1', 1) + ' mg: el truncado de cada muestra es lo que da a ENMO su valor.</p>');
    if (st.offset > 0) p.push('<p>Con ' + num(st.offset, 0) + ' mg de error en el eje vertical, el niño de pie y quieto «hace» ' + f('enmo', 0) + ' mg de ENMO: actividad ficticia. Tumbado, el mismo error apenas cuenta (la gravedad cae sobre otro eje). MAD casi no cambia.</p>');
    else p.push('<p>Sin error de calibración, ENMO en reposo no es exactamente 0: el truncado se queda con la mitad positiva del ruido y de la cuantización (1/256 g ≈ 3,9 mg).</p>');
    $('m-lectura').innerHTML = p.join('');
  }

  function render() { medir(); escribir(); var a = M.norma(M.generarMinuto(st.offset / 1000)), ep = M.metricas(a); dibujar(a, ep); lectura(ep); }

  $('m-off').addEventListener('input', function () { st.offset = Number($('m-off').value); render(); });
  $('m-0').addEventListener('click', function () { st.offset = 0; render(); });
  $('m-20').addEventListener('click', function () { st.offset = 20; render(); });
  $('m-50').addEventListener('click', function () { st.offset = 50; render(); });
  var tm = null; window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  render();
})();
