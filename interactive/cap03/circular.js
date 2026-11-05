/* Núcleo numérico de la simulación «Del momento de cada episodio al vector». Script clásico (file:// y Node), sin dependencias.
 *
 * Cada episodio vigoroso tiene una hora h relativa al punto medio de L5 (0–24 h) y un peso w (su dosis). Su ángulo es
 * φ = 2π·h/24. La media circular ponderada es el vector (C, S) = Σ w·(cos φ, sin φ) / Σ w; su módulo es R (0–1, concentración)
 * y su ángulo, θ = atan2(S, C) pasado a horas (0–24 h). El modelo de la hipótesis usa las coordenadas R·cos θ = C y
 * R·sin θ = S, que son números corrientes: 0 h y 24 h dan lo mismo, y un niño sin preferencia horaria (R ≈ 0) queda cerca
 * del origen en lugar de recibir un θ arbitrario.
 * Los episodios son SINTÉTICOS (normal envuelta con semilla fija): no son episodios de ningún registro. */
(function (root) {
  'use strict';
  var SEMILLA = 2026, N = 40;

  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function mod24(h) { return ((h % 24) + 24) % 24; }

  // N episodios: hora = media + dispersión·z (z normal estándar, la misma para cada semilla), envuelta en [0, 24); peso 1–3.
  function episodios(media, dispersion, semilla) {
    var r = mulberry32(semilla === undefined ? SEMILLA : semilla), e = [];
    for (var i = 0; i < N; i++) {
      var u1 = 1 - r(), u2 = r(), z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
      e.push({ hora: mod24(media + dispersion * z), peso: 1 + 2 * r() });
    }
    return e;
  }

  function mediaCircular(e) {
    var C = 0, S = 0, W = 0;
    e.forEach(function (p) { var f = 2 * Math.PI * p.hora / 24; C += p.peso * Math.cos(f); S += p.peso * Math.sin(f); W += p.peso; });
    C /= W; S /= W;
    var R = Math.sqrt(C * C + S * S);
    return { C: C, S: S, R: R, theta: mod24(Math.atan2(S, C) * 24 / (2 * Math.PI)) };
  }

  function mediaLineal(e) {
    var s = 0, W = 0; e.forEach(function (p) { s += p.peso * p.hora; W += p.peso; }); return s / W;
  }

  // Distancia en el reloj entre dos horas (0–12 h).
  function distanciaHoras(a, b) { var d = Math.abs(mod24(a) - mod24(b)); return Math.min(d, 24 - d); }

  var api = { SEMILLA: SEMILLA, N: N, episodios: episodios, mediaCircular: mediaCircular, mediaLineal: mediaLineal,
              distanciaHoras: distanciaHoras, mod24: mod24 };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.BRCap03Circular = api;
})(typeof window !== 'undefined' ? window : this);
