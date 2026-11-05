/* Simulación «Misma composición, distinto momento». Datos sintéticos (interactive/cap03/luz24.js).
 * No lee nada del proyecto ni de data/. Script clásico, sin red. */
(function () {
  'use strict';
  var L = window.BRCap03Luz24, raiz = document.getElementById('demo-luz24');
  if (!L || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };
  var DEF = { inicio: 10 * 60, duracion: 120 };
  var st = Object.assign({}, DEF);
  var INI_MIN = 7 * 60, FIN = L.SUPUESTOS.acostarse;

  function num(v, d) { return v.toFixed(d).replace('-', '−').replace('.', ','); }
  function hhmm(min) { var h = Math.floor(min / 60), m = Math.round(min - 60 * h); if (m === 60) { h++; m = 0; } return (h < 10 ? '0' : '') + h + ':' + (m < 10 ? '0' : '') + m; }

  function escribir() {
    st.duracion = Math.min(st.duracion, FIN - INI_MIN);
    st.inicio = Math.max(INI_MIN, Math.min(st.inicio, FIN - st.duracion));
    $('l-ini').max = FIN - st.duracion; $('l-ini').value = st.inicio; $('l-ini-val').textContent = hhmm(st.inicio) + '–' + hhmm(st.inicio + st.duracion);
    $('l-dur').value = st.duracion; $('l-dur-val').textContent = st.duracion + ' min';
  }

  var W = 760, Lm = 64, Rm = 14;
  function medir() { var w = $('l-plot').getBoundingClientRect().width; W = Math.max(300, Math.round(w || 760)); Lm = W < 480 ? 50 : 64; Rm = W < 480 ? 8 : 14; }

  function dibujar(x, c) {
    var estrecho = W < 520, H = estrecho ? 300 : 270, T = 14, h = 130, s = '', i;
    function X(min) { return Lm + min / 1440 * (W - Lm - Rm); }
    function Y(lx) { return T + h - Math.log10(1 + lx) / Math.log10(1 + 3000) * h; }
    s += '<rect class="ventana-dia" x="' + X(360) + '" y="' + T + '" width="' + (X(1080) - X(360)) + '" height="' + h + '"/>';
    s += '<text x="' + ((X(360) + X(1080)) / 2) + '" y="' + (T + 12) + '" text-anchor="middle" font-size="11" class="tramo">ventana «día» (06–18 h)</text>';
    [0, 1, 100, 3000].forEach(function (v) {
      s += '<line class="rejilla" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y(v) + '" y2="' + Y(v) + '"/><text x="' + (Lm - 6) + '" y="' + (Y(v) + 4) + '" text-anchor="end" font-size="12">' + v + ' lx</text>';
    });
    for (var hh = 0; hh <= 24; hh += (W < 480 ? 6 : 3)) s += '<text x="' + X(hh * 60) + '" y="' + (T + h + 16) + '" text-anchor="middle" font-size="12">' + hh + '</text>';
    var d = '';
    for (i = 0; i < 1440; i++) d += (i ? 'L' : 'M') + X(i).toFixed(1) + ',' + Y(x[i]).toFixed(1) + 'L' + X(i + 1).toFixed(1) + ',' + Y(x[i]).toFixed(1);
    s += '<path class="curva-2" d="' + d + '"/>';
    s += '<rect x="' + Lm + '" y="' + T + '" width="' + (W - Lm - Rm) + '" height="' + h + '" fill="none" class="eje"/>';
    // Barra de composición (minutos en los cuatro niveles de Ding et al.)
    var by = T + h + 44, bh = 26, partes = [['diaAlta', 'día ≥ 100 lx', 'c1'], ['diaBaja', 'día < 100 lx', 'c2'], ['nocheAlta', 'noche ≥ 1 lx', 'c3'], ['nocheBaja', 'noche < 1 lx', 'c4']], xa = Lm;
    s += '<text x="' + (estrecho ? 8 : Lm) + '" y="' + (by - 6) + '" font-size="12">' + (estrecho ? 'Composición de las 24 h' : 'Composición de las 24 h (lo que entra en el modelo de Ding et al.)') + '</text>';
    partes.forEach(function (p, k) {
      var w = c[p[0]] / 1440 * (W - Lm - Rm);
      s += '<rect class="comp ' + p[2] + '" x="' + xa.toFixed(1) + '" y="' + by + '" width="' + Math.max(0, w).toFixed(1) + '" height="' + bh + '"/>';
      var lx = estrecho ? 8 : Lm + (k % 2) * (W - Lm - Rm) / 2, ly = by + bh + 18 + (estrecho ? k : Math.floor(k / 2)) * 18;
      s += '<rect class="comp ' + p[2] + '" x="' + lx + '" y="' + (ly - 10) + '" width="12" height="12"/><text x="' + (lx + 18) + '" y="' + ly + '" font-size="12">' + p[1] + ': ' + c[p[0]] + ' min</text>';
      xa += w;
    });
    $('l-plot').setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    $('l-plot').innerHTML = s;
  }

  function lectura(x, c) {
    var ref = L.composicion(L.dia(DEF)), igual = ['diaAlta', 'diaBaja', 'nocheAlta', 'nocheBaja'].every(function (k) { return c[k] === ref[k]; });
    var tras18 = Math.max(0, st.inicio + st.duracion - 1080) - Math.max(0, st.inicio - 1080), p = [];
    p.push('<p><strong>Hora media de la luz ≥ 100 lx:</strong> ' + hhmm(L.horaMedia(x, 100) * 60) + '. <strong>Composición:</strong> ' + c.diaAlta + ' / ' + c.diaBaja + ' / ' + c.nocheAlta + ' / ' + c.nocheBaja + ' min (día ≥ 100, día < 100, noche ≥ 1, noche < 1).</p>');
    if (igual) p.push('<p>La composición es <strong>idéntica</strong> a la de la luz exterior a media mañana (10:00–12:00): para el modelo composicional, los dos días son el mismo, aunque la luz llegue a otra hora.</p>');
    if (tras18 > 0) p.push('<p><strong>' + Math.min(tras18, st.duracion) + ' min de luz exterior caen después de las 18:00</strong>: cuentan como «noche ≥ 1 lx», igual que la lámpara del salón. La ventana de reloj fija no sabe si ese niño ya está cenando o aún en el parque.</p>');
    else if (!igual) p.push('<p>La composición cambia solo porque cambia la duración de la luz exterior, no su momento.</p>');
    $('l-lectura').innerHTML = p.join('');
  }

  function render() { medir(); escribir(); var x = L.dia(st), c = L.composicion(x); dibujar(x, c); lectura(x, c); }

  $('l-ini').addEventListener('input', function () { st.inicio = Number($('l-ini').value); render(); });
  $('l-dur').addEventListener('input', function () { st.duracion = Number($('l-dur').value); render(); });
  $('l-man').addEventListener('click', function () { st = { inicio: 600, duracion: 120 }; render(); });
  $('l-tar').addEventListener('click', function () { st = { inicio: 960, duracion: 120 }; render(); });
  $('l-cruce').addEventListener('click', function () { st = { inicio: 1020, duracion: 120 }; render(); });
  var tm = null; window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  render();
})();
