/* 粒子示意：把烧杯里"看不见的粒子"画出来。
   完全由 systems.js 的物种数据驱动 —— 新增系统不需要改这里。 */

(function (root) {
  'use strict';

  // 从吸收系数反推粒子的固有色：液体显色越深，粒子越蓝
  function particleColor(sp) {
    if (sp.phase === 'solid') {
      var s = sp.scat;
      return [Math.round(s[0] * 235), Math.round(s[1] * 235), Math.round(s[2] * 235)];
    }
    var A = sp.absorb;
    if (A && A[0] + A[1] + A[2] > 0) {
      var t0 = Math.exp(-0.9 * A[0]), t1 = Math.exp(-0.9 * A[1]), t2 = Math.exp(-0.9 * A[2]);
      return [Math.round(t0 * 250), Math.round(t1 * 250), Math.round(t2 * 250)];
    }
    return [172, 172, 178];          // 无色离子：淡灰
  }

  // 一个带高光的球，画出体积感
  function ball(g, x, y, r, col) {
    var grd = g.createRadialGradient(x - r * 0.35, y - r * 0.4, r * 0.1, x, y, r);
    grd.addColorStop(0, 'rgba(255,255,255,0.85)');
    grd.addColorStop(0.35, 'rgba(' + col[0] + ',' + col[1] + ',' + col[2] + ',1)');
    grd.addColorStop(1, 'rgba(' + Math.round(col[0] * 0.55) + ',' +
                        Math.round(col[1] * 0.55) + ',' + Math.round(col[2] * 0.55) + ',1)');
    g.fillStyle = grd;
    g.beginPath(); g.arc(x, y, r, 0, 6.284); g.fill();
    g.strokeStyle = 'rgba(28,26,23,0.28)';
    g.lineWidth = 0.7;
    g.stroke();
  }

  // 挑出值得展示的物种：金属离子、碱、主要固相、络离子
  function pick(sys) {
    var out = [], seen = {};
    function push(k) {
      if (!k || seen[k]) return;
      var sp = sys.byKey[k];
      if (sp) { seen[k] = 1; out.push(sp); }
    }
    push(sys.metal);
    push(sys.base);
    for (var i = 0; i < sys.solids.length; i++) push(sys.solids[i].key);
    for (i = 0; i < sys.solutes.length; i++) {
      if (sys.solutes[i].flag) push(sys.solutes[i].key);
    }
    return out.slice(0, 4);
  }

  function draw(canvas, sys, t) {
    var g = canvas.getContext('2d');
    var dpr = Math.min(2, window.devicePixelRatio || 1);
    var w = canvas.clientWidth || 460, h = canvas.clientHeight || 132;
    if (canvas.width !== Math.round(w * dpr)) {
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
    }
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);

    var list = pick(sys);
    if (!list.length) return;
    var n = list.length;
    var cw = w / n;
    if (cw < 56) { h = Math.max(h, 132); }   // 槽太窄时留白由 CSS 兜底
    var font = '"Songti SC","SimSun","Noto Serif SC",Georgia,serif';

    for (var i = 0; i < n; i++) {
      var sp = list[i];
      var col = particleColor(sp);
      var cx = cw * (i + 0.5);
      var cy = h * 0.42;
      var pad = 5;
      var boxW = cw - pad * 2 - 4;

      // 小槽
      g.strokeStyle = 'rgba(28,26,23,0.16)';
      g.lineWidth = 1;
      g.beginPath();
      g.rect(cw * i + pad + 0.5, 4.5, boxW, h * 0.66);
      g.stroke();

      var count = sp.phase === 'solid' ? 5 : 4;
      var rr = Math.min(8.5, boxW * 0.16);
      var spanX = Math.max(6, boxW * 0.30 - rr);
      for (var k = 0; k < count; k++) {
        // 慢慢浮动，固相沉在下面
        var ph = t * (0.35 + 0.12 * k) + k * 1.7 + i * 2.1;
        var settle = sp.phase === 'solid' ? 1 : 0;
        var px = cx + Math.sin(ph) * spanX + Math.cos(ph * 0.7) * 5;
        var py = cy + Math.cos(ph * 0.9) * (h * 0.11) + settle * (h * 0.11)
                 + (k % 2) * (settle ? h * 0.06 : 0);
        ball(g, px, py, rr * (0.82 + 0.18 * ((k * 7) % 3) / 2), col);
      }

      g.textAlign = 'center';
      g.textBaseline = 'middle';
      g.font = '11.5px ' + font;
      g.fillStyle = 'rgba(28,26,23,0.85)';
      g.fillText(sp.name, cx, h * 0.82);
      g.font = '9.5px "Cascadia Mono",Consolas,monospace';
      g.fillStyle = 'rgba(28,26,23,0.45)';
      var tag = sp.phase === 'solid' ? '固相' : '溶液';
      g.fillText(tag, cx, h * 0.94);
    }
  }

  root.MoleculeView = { draw: draw };
})(window);