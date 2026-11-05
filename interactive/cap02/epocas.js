/* Núcleo numérico de la simulación «Del segundo al minuto». Script clásico (file:// y Node), sin dependencias.
 *
 * Genera UNA HORA SINTÉTICA de ENMO (mg) segundo a segundo con una semilla fija, la promedia en epochs de 1, 5, 15 o 60 s
 * y cuenta el tiempo de los epochs cuya media llega al umbral. Así se ve por qué el protocolo (duración del epoch y punto
 * de corte) cambia los minutos de actividad vigorosa con los mismos movimientos.
 * No es un registro: los tramos, las ráfagas y sus amplitudes son supuestos de la simulación (ver SUPUESTOS). No trabaja
 * con counts de ActiLife: los umbrales son de ENMO, como en Brincadiano. */
(function (root) {
  'use strict';
  var SEMILLA = 2026;
  // Tramos de la hora (segundos): reposo sentado, juego con ráfagas, carrera mantenida de 2 min.
  var TRAMOS = [[0, 600, 'sentado'], [600, 1500, 'juego'], [1500, 1620, 'carrera'], [1620, 2400, 'sentado'],
                [2400, 3300, 'juego'], [3300, 3600, 'sentado']];
  var SUPUESTOS = {
    sentado: [8, 23],            // mg, uniforme
    juegoBase: [60, 240],        // mg, paseo aleatorio acotado
    rafagaProb: 1 / 40,          // probabilidad por segundo de que empiece una ráfaga durante el juego
    rafagaMedia: 4,              // s; duración 1 + exponencial redondeada, como máximo 20 s
    rafagaMax: 20,
    rafagaAmp: [550, 1200],      // mg
    carrera: [620, 840]          // mg
  };
  var EPOCAS = [1, 5, 15, 60];
  var UMBRAL_VIGOROSA = 464.6;   // mg: Hildebrand, cadera, 7–11 años (cap. 10 del informe)
  var UMBRAL_MODERADA = 142.6;

  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function generarHora(semilla) {
    var r = mulberry32(semilla === undefined ? SEMILLA : semilla), S = SUPUESTOS;
    var x = new Float64Array(3600), rafagas = [], base = 120, resto = 0, amp = 0;
    TRAMOS.forEach(function (tr) {
      for (var s = tr[0]; s < tr[1]; s++) {
        if (tr[2] === 'sentado') { x[s] = S.sentado[0] + (S.sentado[1] - S.sentado[0]) * r(); resto = 0; continue; }
        if (tr[2] === 'carrera') { x[s] = S.carrera[0] + (S.carrera[1] - S.carrera[0]) * r(); resto = 0; continue; }
        base = Math.min(S.juegoBase[1], Math.max(S.juegoBase[0], base + (r() - 0.5) * 30));
        if (resto === 0 && r() < S.rafagaProb) {
          resto = Math.min(S.rafagaMax, 1 + Math.round(-Math.log(1 - r()) * (S.rafagaMedia - 1)));
          amp = S.rafagaAmp[0] + (S.rafagaAmp[1] - S.rafagaAmp[0]) * r();
          rafagas.push({ inicio: s, duracion: Math.min(resto, tr[1] - s) });
        }
        if (resto > 0) { x[s] = amp * (0.85 + 0.3 * r()); resto--; } else x[s] = base + 20 * (r() - 0.5);
      }
    });
    return { enmo: x, rafagas: rafagas };
  }

  // Media de cada epoch de `epoca` segundos (3600 es múltiplo de 1, 5, 15 y 60).
  function agregar(x, epoca) {
    var n = Math.floor(x.length / epoca), m = new Float64Array(n);
    for (var i = 0; i < n; i++) { var s = 0; for (var j = 0; j < epoca; j++) s += x[i * epoca + j]; m[i] = s / epoca; }
    return m;
  }

  // Minutos de los epochs cuya media llega al umbral.
  function minutosSobre(x, epoca, umbral) {
    var m = agregar(x, epoca), k = 0;
    for (var i = 0; i < m.length; i++) if (m[i] >= umbral) k++;
    return k * epoca / 60;
  }

  function maximo(a) { var M = -Infinity; for (var i = 0; i < a.length; i++) if (a[i] > M) M = a[i]; return M; }

  var api = { SEMILLA: SEMILLA, TRAMOS: TRAMOS, SUPUESTOS: SUPUESTOS, EPOCAS: EPOCAS, UMBRAL_VIGOROSA: UMBRAL_VIGOROSA,
              UMBRAL_MODERADA: UMBRAL_MODERADA, generarHora: generarHora, agregar: agregar, minutosSobre: minutosSobre,
              maximo: maximo };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.BRCap02Epocas = api;
})(typeof window !== 'undefined' ? window : this);
