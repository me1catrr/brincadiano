/* Núcleo numérico de la simulación «Misma composición, distinto momento». Script clásico (file:// y Node), sin dependencias.
 *
 * Un día SINTÉTICO de luz, minuto a minuto (1440 valores, en lux): a oscuras durmiendo (0 lx), luz de interior de 07:00 a
 * 21:30 y un bloque de luz brillante (exterior) que se puede mover. Se resume como Ding et al. (2025): minutos en cuatro
 * niveles con ventanas de reloj fijas, luz diurna (06–18 h) ≥ 100 lx y < 100 lx, y luz nocturna (18–06 h) ≥ 1 lx y < 1 lx.
 * Y, además, con una medida de momento: la hora media de los minutos con luz ≥ 100 lx.
 * No es un registro: los niveles de luz son supuestos de la simulación (SUPUESTOS). */
(function (root) {
  'use strict';
  var SUPUESTOS = { despertar: 7 * 60, acostarse: 21 * 60 + 30, interior: 60, exterior: 3000 };

  // brillo: { inicio (min desde las 00:00), duracion (min) }.
  function dia(brillo) {
    var S = SUPUESTOS, x = new Float64Array(1440);
    for (var m = 0; m < 1440; m++) {
      if (m < S.despertar || m >= S.acostarse) { x[m] = 0; continue; }
      x[m] = S.interior;
      if (m >= brillo.inicio && m < brillo.inicio + brillo.duracion) x[m] = S.exterior;
    }
    return x;
  }

  function esDia(m) { return m >= 360 && m < 1080; }   // 06:00–18:00

  function composicion(x) {
    var c = { diaAlta: 0, diaBaja: 0, nocheAlta: 0, nocheBaja: 0 };
    for (var m = 0; m < 1440; m++) {
      if (esDia(m)) { if (x[m] >= 100) c.diaAlta++; else c.diaBaja++; }
      else { if (x[m] >= 1) c.nocheAlta++; else c.nocheBaja++; }
    }
    return c;
  }

  // Hora media (h) de los minutos con luz ≥ umbral; NaN si no hay ninguno.
  function horaMedia(x, umbral) {
    var s = 0, n = 0;
    for (var m = 0; m < 1440; m++) if (x[m] >= umbral) { s += m + 0.5; n++; }
    return n ? s / n / 60 : NaN;
  }

  var api = { SUPUESTOS: SUPUESTOS, dia: dia, composicion: composicion, horaMedia: horaMedia };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.BRCap03Luz24 = api;
})(typeof window !== 'undefined' ? window : this);
