/* Núcleo numérico de la simulación «¿Qué frecuencias sobreviven?». Script clásico (file:// y Node), sin dependencias.
 *
 * (1) Filtro de los counts, IDEALIZADO: un paso alto de Butterworth de orden 2 a 0,29 Hz en cascada con un paso bajo de
 *     orden 2 a 1,63 Hz, las frecuencias de corte aproximadas que da el informe. Ganancia:
 *         |H(f)| = (f/fl)² / √(1 + (f/fl)⁴) · 1 / √(1 + (f/fh)⁴)
 *     NO es el filtro de ActiGraph, que tiene otra forma (Neishabouri et al., 2022); solo muestra el orden de magnitud.
 * (2) Muestreo a fs = 30 Hz sin filtro previo: una sinusoide de frecuencia f aparece con la frecuencia
 *         f_aparente = |f − fs · redondeo(f / fs)|,
 *     igual a f por debajo de Nyquist (15 Hz) y «disfrazada» de más lenta por encima. */
(function (root) {
  'use strict';
  var FL = 0.29, FH = 1.63, FS = 30;

  function ganancia(f) {
    var a = Math.pow(f / FL, 2), b = Math.pow(f / FH, 4);
    return a / Math.sqrt(1 + a * a) / Math.sqrt(1 + b);
  }
  function dB(g) { return 20 * Math.log10(g); }
  function aparente(f, fs) { fs = fs || FS; return Math.abs(f - fs * Math.round(f / fs)); }

  var api = { FL: FL, FH: FH, FS: FS, ganancia: ganancia, dB: dB, aparente: aparente };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.BRCap04Frecuencias = api;
})(typeof window !== 'undefined' ? window : this);
