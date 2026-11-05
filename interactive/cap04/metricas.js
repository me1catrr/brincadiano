/* Núcleo numérico de la simulación «ENMO, MAD y un error de calibración». Script clásico (file:// y Node), sin dependencias.
 *
 * Un MINUTO SINTÉTICO de aceleración triaxial a 30 Hz (en g): 20 s de pie y quieto (la gravedad cae sobre el eje vertical, y),
 * 20 s andando y 20 s tumbado (la gravedad cae sobre z). Ruido gaussiano en cada eje y cuantización a 1/256 g, como el
 * wGT3X-BT. Se le puede sumar un error de desplazamiento (offset) al eje y, como el de un dispositivo mal calibrado.
 * Por epoch de 5 s se calculan:
 *   ENMO      = media de max(‖a‖ − 1, 0)        (truncado en cada muestra, antes de promediar)
 *   ‖a‖ − 1   = media de ‖a‖, menos 1           (lo que NO es ENMO)
 *   MAD       = media de |‖a‖ − media de ‖a‖|   (desviación absoluta media de la norma en el epoch)
 * No es un registro: las señales son sinusoides y ruido inventados (ver SUPUESTOS). */
(function (root) {
  'use strict';
  var FS = 30, SEMILLA = 2026, EPOCH = 5;
  var SUPUESTOS = { ruido: 0.003, lsb: 1 / 256, paso: 2 };   // g, g, Hz (frecuencia de paso al andar)
  var TRAMOS = [[0, 20, 'de pie, quieto'], [20, 40, 'andando'], [40, 60, 'tumbado']];

  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // offset: error de desplazamiento del eje y, en g.
  function generarMinuto(offset, semilla) {
    var r = mulberry32(semilla === undefined ? SEMILLA : semilla), n = 60 * FS, S = SUPUESTOS;
    var x = new Float64Array(n), y = new Float64Array(n), z = new Float64Array(n);
    function gauss() { var u = 1 - r(), v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
    function q(v) { return Math.round(v / S.lsb) * S.lsb; }
    for (var i = 0; i < n; i++) {
      var t = i / FS, w = 2 * Math.PI * S.paso * t, gx = 0, gy = 0, gz = 0;
      if (t < 20) gy = 1;
      else if (t < 40) { gx = 0.15 * Math.sin(w); gy = 1 + 0.35 * Math.sin(w) + 0.12 * Math.sin(2 * w + 0.5); gz = 0.10 * Math.sin(w / 2); }
      else gz = 1;
      x[i] = q(gx + S.ruido * gauss()); y[i] = q(gy + offset + S.ruido * gauss()); z[i] = q(gz + S.ruido * gauss());
    }
    return { x: x, y: y, z: z };
  }

  function norma(s) {
    var n = s.x.length, a = new Float64Array(n);
    for (var i = 0; i < n; i++) a[i] = Math.sqrt(s.x[i] * s.x[i] + s.y[i] * s.y[i] + s.z[i] * s.z[i]);
    return a;
  }

  // Métricas por epoch de 5 s, en mg.
  function metricas(a) {
    var m = EPOCH * FS, k = Math.floor(a.length / m), out = [];
    for (var e = 0; e < k; e++) {
      var sumN = 0, sumE = 0, j;
      for (j = 0; j < m; j++) { var v = a[e * m + j]; sumN += v; sumE += Math.max(v - 1, 0); }
      var media = sumN / m, mad = 0;
      for (j = 0; j < m; j++) mad += Math.abs(a[e * m + j] - media);
      out.push({ t0: e * EPOCH, enmo: 1000 * sumE / m, normaMenos1: 1000 * (media - 1), mad: 1000 * mad / m });
    }
    return out;
  }

  // Media de una métrica sobre los epochs de un tramo [t0, t1) en s.
  function mediaTramo(ep, campo, t0, t1) {
    var s = 0, n = 0; ep.forEach(function (e) { if (e.t0 >= t0 && e.t0 < t1) { s += e[campo]; n++; } }); return s / n;
  }

  var api = { FS: FS, EPOCH: EPOCH, SEMILLA: SEMILLA, SUPUESTOS: SUPUESTOS, TRAMOS: TRAMOS, generarMinuto: generarMinuto,
              norma: norma, metricas: metricas, mediaTramo: mediaTramo };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.BRCap04Metricas = api;
})(typeof window !== 'undefined' ? window : this);
