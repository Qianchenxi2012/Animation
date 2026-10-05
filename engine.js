/* 反应引擎：由 systems.js 的数据驱动。
   场量模型 —— 每一格记录各物种浓度（mol·L-1），每帧依次做
   平流 → 扩散 → 局部反应 → 湍流松弛 → 沉降。
   token 台账：每个物种声明自己结合了多少 token，反应只搬运物种，
   守恒式由 holds 自动给出，不需要任何特判。 */

(function (root) {
  'use strict';

  var NX = 56, NY = 92, N = NX * NY;
  var ROWS_PER_ML = 0.62;
  var CV = 1 / (ROWS_PER_ML * NX);
  var MAXV = 132;
  var POUR_RATE = 17;

  var T_ROOM = 20, T_BOIL = 100;
  var HEAT_RATE = 9, COOL_RATE = 4;

  var BE_W = 176, BE_H = 290, WALL = 7;
  var ML = 56, MR = 126, TOP = 100, BURNER_H = 46;
  var BG = [233, 229, 219];
  var FONT = '"Songti SC","SimSun","Noto Serif SC",Georgia,serif';
  var FONT2 = '"Cascadia Mono",Consolas,monospace';

  // exp 查表
  var LUT_N = 2048, LUT_MAX = 14, LUT_K = (LUT_N - 1) / LUT_MAX;
  var LUT_E = new Float32Array(LUT_N);
  for (var li = 0; li < LUT_N; li++) LUT_E[li] = Math.exp(-li / LUT_K);
  function lutE(v) {
    if (v <= 0) return 1;
    if (v >= LUT_MAX) return 0;
    return LUT_E[(v * LUT_K) | 0];
  }

  // 絮状物的平滑值噪声
  var noise = new Float32Array(N);
  (function () {
    var s = 20240117;
    function rnd() { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; }
    function octave(step, amp) {
      var gx = Math.ceil(NX / step) + 2, gy = Math.ceil(NY / step) + 2;
      var lat = new Float32Array(gx * gy), i;
      for (i = 0; i < lat.length; i++) lat[i] = rnd() - 0.5;
      for (var y = 0; y < NY; y++) {
        var fy = y / step, y0 = fy | 0, ty = fy - y0;
        ty = ty * ty * (3 - 2 * ty);
        var r0 = y0 * gx, r1 = (y0 + 1) * gx;
        for (var x = 0; x < NX; x++) {
          var fx = x / step, x0 = fx | 0, tx = fx - x0;
          tx = tx * tx * (3 - 2 * tx);
          var a = lat[r0 + x0], b = lat[r0 + x0 + 1];
          var c = lat[r1 + x0], e = lat[r1 + x0 + 1];
          var top = a + (b - a) * tx, bot = c + (e - c) * tx;
          noise[y * NX + x] += (top + (bot - top) * ty) * amp;
        }
      }
    }
    octave(9, 0.7);
    octave(4.5, 0.3);
  })();

  /* ============================ 引擎 ============================ */

  function createEngine(canvas) {
    var ctx = canvas.getContext('2d');
    var off = document.createElement('canvas');
    off.width = NX; off.height = NY;
    var octx = off.getContext('2d');
    var img = octx.createImageData(NX, NY);

    var tmp = new Float32Array(N), dv = new Float32Array(N);
    var velX = new Float32Array(N), velY = new Float32Array(N);
    var wetRow = new Uint8Array(NY + 2);
    var wbuf = new Float32Array(128);

    var F = {}, S = null;          // 物种场 / 当前系统
    var BX = 0, BY = 0, SC = 1, dpr = 1;

    var V = 0, LV = 0;
    var temp = T_ROOM, heating = false, boil = 0;
    var swirl = 0, simT = 0, jetT = 0, jetPos = NX / 2;
    var labelsOn = 0;

    var pour = null;               // {key, remain, total, pos, conc}
    var phase = 'pour', phaseKey = null, wait = 2.0, queue = [];

    /* ---------- 载入系统 ---------- */

    function load(sys, preset) {
      S = sys;
      F = {};
      S.solutes = []; S.solids = [];
      S.byKey = {}; S.reagentByKey = {};
      S.tokens = {};
      sys.reagents.forEach(function (rd) { S.reagentByKey[rd.key] = rd; });
      for (var i = 0; i < sys.species.length; i++) {
        var sp = sys.species[i];
        var f = new Float32Array(N);
        F[sp.key] = f;
        S.byKey[sp.key] = sp;
        sp.field = f;
        if (sp.phase === 'solid') S.solids.push(sp); else S.solutes.push(sp);
        if (sp.holds) for (var tk in sp.holds) S.tokens[tk] = true;
      }
      S.tokenList = Object.keys(S.tokens);
      S.orderKey = sys.orders[0].key;
      compileReactions();
      restart(preset);
    }

    function compileReactions() {
      S.rx = (S.reactions || []).map(function (R) {
        var rf = [], rn = [], pf = [], pn = [], k;
        for (var a in R.r) { rf.push(F[a]); rn.push(R.r[a]); }
        for (var b in R.p) { pf.push(F[b]); pn.push(R.p[b]); }
        k = (typeof R.k === 'function') ? R.k : null;
        return { rf: rf, rn: rn, pf: pf, pn: pn, nk: rn.length, np: pn.length,
                 kNum: typeof R.k === 'number' ? R.k : 0, kFn: k,
                 gate: R.gate || null, t0: R.t0, thresh: R.thresh,
                 grp: R.grp || null, when: R.when, scarce: R.scarce,
                 kdt: 0, off: false };
      });
    }

        // 计划投料量（实验开始前就已知），用来判断产物类别
    function planned() {
      var out = {};
      for (var t = 0; t < S.tokenList.length; t++) out[S.tokenList[t]] = 0;
      for (var r = 0; r < S.reagents.length; r++) {
        var rd = S.reagents[r];
        var n = rd.conc * rd.vol / 1000;
        for (var k in rd.injects) {
          var h = S.byKey[k].holds;
          if (!h) continue;
          for (var tk in h) out[tk] += n * h[tk];
        }
      }
      return out;
    }

    function tokenOf(spKey) {
      var h = S.byKey[spKey] && S.byKey[spKey].holds;
      return h ? Object.keys(h)[0] : spKey;
    }

    function updateGates() {
      var pl = planned();
      // 总碱够 2:1 就不会有碱式硫酸铜
      var b = pl[tokenOf(S.base)] || 0, m = pl[tokenOf(S.metal)] || 0;
      S.planBasic = b < S.ratio.stoich * m;
      for (var i = 0; i < S.rx.length; i++) {
        S.rx[i].off = (S.rx[i].grp === 'basic') && !S.planBasic;
      }
    }

    function restart(preset) {
      for (var k in F) F[k].fill(0);
      V = 0; temp = T_ROOM; heating = false; boil = 0;
      swirl = 0; simT = 0; labelsOn = 0; pour = null;
      if (preset) {
        for (var rk in preset.set) {
          var rd = S.reagentByKey[rk];
          rd.conc = preset.set[rk][0];
          rd.vol = preset.set[rk][1];
        }
      }
      updateGates();
      var seq = currentOrder().seq;
      queue = seq.slice(1);
      phase = 'pour'; phaseKey = seq[0];
      startPour(seq[0]);
    }

    function currentOrder() {
      var o = null;
      for (var i = 0; i < S.orders.length; i++) if (S.orders[i].key === S.orderKey) o = S.orders[i];
      return o || S.orders[0];
    }

    function startPour(rkey) {
      var rd = S.reagentByKey[rkey];
      var idx = S.reagents.indexOf(rd);
      var pos = NX / 2 - 9 + (idx % 2) * 17;
      pour = { key: rkey, remain: rd.vol, total: rd.vol, pos: pos, conc: rd.conc, rd: rd };
      S.orderKeyUsed = rkey;
    }

    function reagentByKey(k) { return S.reagentByKey[k]; }

    /* ---------- 几何 ---------- */

    function wet(y) { return y >= 0 && y < NY && (NY - 1 - y) < LV; }
    function markWet() {
      var top = NY - Math.ceil(LV);
      for (var y = 0; y < NY; y++) wetRow[y] = (y >= top) ? 1 : 0;
    }
    function levelPx() { return BY + (NY - Math.min(NY, V * ROWS_PER_ML)) / NY * BE_H * SC; }

    function gauss(x, s) { return Math.exp(-(x * x) / (2 * s * s)); }

    /* ---------- 注入 ---------- */

    function inject(fx, fy, mol, field) {
      var x0 = Math.max(0, Math.floor(fx - 3.2)), x1 = Math.min(NX - 1, Math.ceil(fx + 3.2));
      var y0 = Math.max(0, Math.floor(fy - 2.6)), y1 = Math.min(NY - 1, Math.ceil(fy + 2.6));
      var s = 0, k = 0;
      for (var y = y0; y <= y1; y++) {
        if (!wet(y)) continue;
        for (var x = x0; x <= x1; x++) {
          var w = gauss(x + 0.5 - fx, 2.0) * gauss(y + 0.5 - fy, 1.7);
          wbuf[k++] = w; s += w;
        }
      }
      if (s <= 0) return;
      var inv = mol * 1000 / CV / s;
      k = 0;
      for (var y2 = y0; y2 <= y1; y2++) {
        if (!wet(y2)) continue;
        for (var x2 = x0; x2 <= x1; x2++) field[y2 * NX + x2] += wbuf[k++] * inv;
      }
    }

    /* ---------- 输运 ---------- */

    function buildVel() {
      var depth = Math.max(3, LV);
      var roll = (6.5 + 30 * Math.exp(-simT / 2.2)) * Math.min(1.6, depth / 9);
      for (var y = 0; y < NY; y++) {
        var dRow = NY - 1 - y;
        var fy = 1 - dRow / depth; if (fy < 0) fy = 0; if (fy > 1) fy = 1;
        for (var x = 0; x < NX; x++) {
          var i = y * NX + x;
          var fx = (x + 0.5) / NX;
          var vx = roll * Math.sin(Math.PI * fx) * Math.cos(Math.PI * fy);
          var vy = -roll * Math.cos(Math.PI * fx) * Math.sin(Math.PI * fy);
          if (jetT > 0) vy -= jetT * 46 * gauss(x + 0.5 - jetPos, 3.0) * Math.exp(-dRow / 30);
          if (boil > 0) {
            var rw = boil * 9 * Math.sin(Math.PI * fx) * Math.sin(Math.PI * fy);
            vx += rw * 0.7; vy -= rw;
          }
          if (swirl > 1e-3) {
            var cx = NX / 2, cy = NY - depth / 2;
            var dx = (x + 0.5 - cx) / depth, dy = (y + 0.5 - cy) / depth;
            var damp = Math.sin(Math.PI * fx) * Math.sin(Math.PI * Math.max(0.05, fy));
            vx += swirl * 26 * (-dy) * damp;
            vy += swirl * 26 * (dx) * damp;
          }
          velX[i] = vx; velY[i] = vy;
        }
      }
    }

    function sample(f, x, y) {
      var x0 = Math.floor(x), y0 = Math.floor(y);
      if (x0 < 0) x0 = 0; else if (x0 > NX - 2) x0 = NX - 2;
      if (y0 < 0) y0 = 0; else if (y0 > NY - 2) y0 = NY - 2;
      var tx = x - x0, ty = y - y0;
      var a = f[y0 * NX + x0], b = f[y0 * NX + x0 + 1];
      var c = f[(y0 + 1) * NX + x0], d = f[(y0 + 1) * NX + x0 + 1];
      return (a + (b - a) * tx) + ((c + (d - c) * tx) - (a + (b - a) * tx)) * ty;
    }

    function advect(dt, f, vx, vy) {
      var yTop = NY - Math.ceil(LV);
      if (yTop < 0) yTop = 0;
      var lo = yTop - 0.5;
      for (var y = 0; y < NY; y++) {
        for (var x = 0; x < NX; x++) {
          var i = y * NX + x;
          var sx = x + 0.5 - vx[i] * dt, sy = y + 0.5 - vy[i] * dt;
          if (sx < 0.5) sx = 0.5; else if (sx > NX - 0.5) sx = NX - 0.5;
          if (sy < lo) sy = lo; else if (sy > NY - 0.5) sy = NY - 0.5;
          tmp[i] = sample(f, sx - 0.5, sy - 0.5);
        }
      }
      f.set(tmp);
    }

    function sumWet(f) {
      var s = 0;
      for (var y = 0; y < NY; y++) {
        if (!wetRow[y]) continue;
        var o = y * NX;
        for (var x = 0; x < NX; x++) s += f[o + x];
      }
      return s;
    }

    // 半拉格朗日平流，取样不越过液面，再做有界总量修正
    function advectKeep(dt, f, vx, vy) {
      var s = sumWet(f);
      if (s <= 1e-12) {
        for (var y0 = 0; y0 < NY; y0++) if (!wetRow[y0]) {
          var o0 = y0 * NX;
          for (var q = 0; q < NX; q++) f[o0 + q] = 0;
        }
        return;
      }
      advect(dt, f, vx, vy);
      var s2 = sumWet(f);
      var kk = s2 > 1e-12 ? s / s2 : 0;
      for (var y = 0; y < NY; y++) {
        var o = y * NX;
        if (!wetRow[y]) {
          for (var x0 = 0; x0 < NX; x0++) f[o + x0] = 0;
          continue;
        }
        for (var x = 0; x < NX; x++) f[o + x] *= kk;
      }
    }

    function diffuse(dt, f, D) {
      var kk = Math.min(D, 0.2 / dt) * dt;
      for (var y = 0; y < NY; y++) {
        if (!wetRow[y]) { var o = y * NX; for (var z = 0; z < NX; z++) tmp[o + z] = 0; continue; }
        var upWet = wetRow[y - 1], dnWet = wetRow[y + 1];
        var row = y * NX, up = row - NX, dn = row + NX;
        for (var x = 0; x < NX; x++) {
          var i = row + x, c = f[i];
          var l = (x > 0) ? f[i - 1] : c;
          var r = (x < NX - 1) ? f[i + 1] : c;
          var u = upWet ? f[up + x] : c;
          var d2 = dnWet ? f[dn + x] : c;
          tmp[i] = c + kk * (l + r + u + d2 - 4 * c);
        }
      }
      f.set(tmp);
    }

    function relax(f, rate, dt) {
      var s = 0, n = 0;
      for (var y = 0; y < NY; y++) {
        if (!wetRow[y]) continue;
        var o = y * NX;
        for (var x = 0; x < NX; x++) { s += f[o + x]; n++; }
      }
      if (!n) return;
      var m = s / n, kk = 1 - Math.exp(-rate * dt);
      for (var y2 = 0; y2 < NY; y2++) {
        if (!wetRow[y2]) continue;
        var o2 = y2 * NX;
        for (var x2 = 0; x2 < NX; x2++) f[o2 + x2] += (m - f[o2 + x2]) * kk;
      }
    }

    /* ---------- 沉降 ---------- */

    function packOne(f, density) {
      var x, y, i;
      var kh = 0.5 * (1 - Math.exp(-(7 + swirl * 26) * dtOf));
      dv.fill(0);
      for (y = 0; y < NY; y++) {
        if (!wetRow[y]) continue;
        var row = y * NX;
        for (x = 0; x < NX - 1; x++) {
          var a = f[row + x], b = f[row + x + 1];
          if (a <= 0 && b <= 0) continue;
          var m = (a - b) * kh * 0.5;
          dv[row + x] -= m; dv[row + x + 1] += m;
        }
      }
      for (i = 0; i < N; i++) if (dv[i] !== 0) { var t = f[i] + dv[i]; f[i] = t > 0 ? t : 0; }

      var fr = 8.5 * (1 - swirl * 0.9) * dtOf;
      var up = (swirl * 11 + boil * 5) * dtOf;
      for (y = 0; y < NY; y++) {
        if (!wetRow[y]) continue;
        var below = wetRow[y + 1], above = wetRow[y - 1];
        var base = y * NX;
        for (x = 0; x < NX; x++) {
          i = base + x;
          if (above && up > 0) {
            var u = f[i] * up;
            if (u > 0) { f[i] -= u; f[i - NX] += u; }
          }
          if (f[i] <= 0) continue;
          if (below) {
            var mv = f[i] * fr;
            var room = density - f[i + NX];
            if (room <= 0) continue;
            if (mv > room) mv = room;
            if (mv > 0) { f[i] -= mv; f[i + NX] += mv; }
            continue;
          }
          var mv2 = f[i] * fr;
          if (mv2 <= 0) continue;
          var nl = (x > 0 ? 1 : 0) + (x < NX - 1 ? 1 : 0);
          if (!nl) continue;
          var share = mv2 / nl;
          f[i] -= mv2;
          if (x > 0) f[i - 1] += share;
          if (x < NX - 1) f[i + 1] += share;
          f[i] += mv2 - share * nl;
        }
      }
    }

    var dtOf = 1 / 60;

    function pack() {
      for (var i = 0; i < S.solids.length; i++) {
        packOne(S.solids[i].field, S.solids[i].density || 2.6);
      }
    }

    /* ---------- 反应 ---------- */

    function gateFactor(R, i) {
      var g = R.gate;
      if (!g) return 1;
      if (g === 'metal_rich') {
        if (!S.planBasic) return 0;
        var m = F[S.metal], b = F[S.base];
        return S.ratio.stoich * b[i] < m[i] ? 1 : 0;
      }
      if (g === 'metal_poor') {
        // 总碱充足时只看整体配比，不用局部比例卡住沉淀反应
        if (!S.planBasic) return 1;
        var m2 = F[S.metal], b2 = F[S.base];
        return S.ratio.stoich * b2[i] >= m2[i] ? 1 : 0;
      }
      if (g === 'hot') {
        var ov = (temp - R.t0) / 45;
        return ov > 0 ? ov * ov : 0;
      }
      if (g === 'conc') return F[S.base][i] > R.thresh ? 1 : 0;
      return 1;
    }

    function react(dt) {
      var R, a;
      for (a = 0; a < S.rx.length; a++) {
        R = S.rx[a];
        if (R.grp === 'basic' && !S.planBasic) { R.kdt = 0; continue; }
        R.kdt = R.kFn ? R.kFn(temp) : 1 - Math.exp(-R.kNum * dt);
        if (R.when === 'base_scarce' && !S.base) { R.kdt = 0; R.off = true; }
        else if (R.grp === 'basic' && !S.planBasic) { R.kdt = 0; R.off = true; }
        else R.off = false;
      }
      var baseF = S.base ? F[S.base] : null;
      var scarce = R.when === 'base_scarce' ? 0.4 : 0;   // OH 很稀时门限
      for (var y = 0; y < NY; y++) {
        if (!wetRow[y]) continue;
        for (var x = 0; x < NX; x++) {
          var i = y * NX + x;
          for (a = 0; a < S.rx.length; a++) {
            R = S.rx[a];
            if (R.off) continue;
            if (R.when === 'base_scarce' && baseF[i] > (R.scarce || 0.4)) continue;
            var f = gateFactor(R, i);
            if (f <= 0) continue;
            var ext = Infinity;
            for (var b = 0; b < R.nk; b++) {
              var q = R.rf[b][i] / R.rn[b];
              if (q < ext) ext = q;
            }
            if (!(ext < Infinity)) continue;
            ext *= R.kdt * f;
            if (ext <= 0) continue;
            for (b = 0; b < R.nk; b++) R.rf[b][i] -= R.rn[b] * ext;
            for (b = 0; b < R.np; b++) R.pf[b][i] += R.pn[b] * ext;
          }
        }
      }
    }

    /* ---------- 统计 ---------- */

    var tally = { mol: {}, L: 0 };

    function doTally() {
      var i, s;
      var acc = {};
      for (i = 0; i < S.species.length; i++) acc[S.species[i].key] = 0;
      for (var y = 0; y < NY; y++) {
        if (!wetRow[y]) continue;
        for (var x = 0; x < NX; x++) {
          var j = y * NX + x;
          for (i = 0; i < S.species.length; i++) {
            var f = S.species[i].field;
            var v = f[j];
            if (v > 0) acc[S.species[i].key] += v;
          }
        }
      }
      tally.mol = acc;
      var fct = CV / 1000;
      for (i = 0; i < S.species.length; i++) S.species[i].mol = acc[S.species[i].key] * fct;
      tally.L = V / 1000;
    }

    /* ---------- 单步 ---------- */

    function step(dt) {
      dtOf = dt;
      if (heating) temp = Math.min(T_BOIL, temp + HEAT_RATE * dt);
      else temp = Math.max(T_ROOM, temp - COOL_RATE * dt);
      var ov = (temp - 96) / 6;
      boil = ov > 0 ? Math.min(1, ov) : 0;

      var r = pour, mol = 0, inj = null;
      if (r) {
        var room = MAXV - V;
        if (room <= 0) { pour = null; r = null; }
        else {
          var add = Math.min(r.remain, POUR_RATE * dt, room);
          r.remain -= add; V += add;
          mol = add * r.conc / 1000;
          inj = r;
          if (r.remain <= 1e-6) pour = null;
        }
      }
      LV = Math.min(NY, V * ROWS_PER_ML);
      markWet();
      if (inj && mol > 0) {
        var row = Math.min(NY - 0.6, V * ROWS_PER_ML);
        var ins = S.reagentByKey[inj.key].injects;
        for (var k in ins) inject(inj.pos, NY - row, mol * ins[k], F[k]);
        jetT = 0.55; jetPos = inj.pos;
      }
      if (LV < 0.3) { jetT = 0; swirl *= Math.pow(0.22, dt); doTally(); return; }

      simT += dt;
      buildVel();
      var q;
      for (q = 0; q < S.solutes.length; q++) advectKeep(dt, S.solutes[q].field, velX, velY);

      var base = S.byKey[S.base];
      var mix = (pour || jetT > 0) ? 2.2 : 0.9;
      if (base) {
        diffuse(dt / 2, base.field, 11); diffuse(dt / 2, base.field, 11);
        diffuse(dt / 2, F[S.metal], 7); diffuse(dt / 2, F[S.metal], 7);
      }
      for (q = 0; q < S.solutes.length; q++) {
        var sp = S.solutes[q];
        if (sp === base) continue;
        var D = (sp.key === S.metal) ? 7 : 8;
        diffuse(dt / 2, sp.field, D); diffuse(dt / 2, sp.field, D);
      }
      react(dt);
      for (q = 0; q < S.solutes.length; q++) relax(S.solutes[q].field, mix, dt);
      pack();
      doTally();
      jetT = Math.max(0, jetT - dt * 2.2);
      swirl *= Math.pow(0.22, dt);
    }

    /* ---------- 着色 ---------- */

    function buildImage() {
      var D = img.data;
      var sols = S.solutes, slds = S.solids;
      var ns = sols.length, nd = slds.length;
      for (var i = 0, p = 0; i < N; i++, p += 4) {
        var ar = 0, ag = 0, ab = 0;
        for (var q = 0; q < ns; q++) {
          var f = sols[q].field, v = f[i];
          if (v > 1e-5) {
            var A = sols[q].absorb;
            ar += 1.7 * A[0] * v; ag += 1.7 * A[1] * v; ab += 1.7 * A[2] * v;
          }
        }
        var tr = lutE(ar), tg = lutE(ag), tbb = lutE(ab);
        var tot = 0, sr = 0, sg = 0, sb = 0;
        for (q = 0; q < nd; q++) {
          var sp = slds[q], pz = sp.field[i];
          if (pz <= 1e-4) continue;
          var dn = sp.density || 2.6;
          var tt = pz > 0.3 ? Math.min(1, (pz - 0.3) / (dn - 0.3)) : 0;
          sr += (sp.scat[0] + (sp.scatD[0] - sp.scat[0]) * tt) * pz;
          sg += (sp.scat[1] + (sp.scatD[1] - sp.scat[1]) * tt) * pz;
          sb += (sp.scat[2] + (sp.scatD[2] - sp.scat[2]) * tt) * pz;
          tot += pz;
        }
        if (tot > 0) {
          var it = 1 / tot;
          var nz = 1 + noise[i] * 0.26 * Math.min(1, tot * 2.5);
          var R = sr * it * nz * 255, G = sg * it * nz * 255, B = sb * it * nz * 255;
          var al = 1 - lutE(tot * 1.7), w = 1 - al;
          D[p] = (R * al + tr * BG[0] * w) | 0;
          D[p + 1] = (G * al + tg * BG[1] * w) | 0;
          D[p + 2] = (B * al + tbb * BG[2] * w) | 0;
        } else {
          D[p] = (tr * BG[0]) | 0;
          D[p + 1] = (tg * BG[1]) | 0;
          D[p + 2] = (tbb * BG[2]) | 0;
        }
        D[p + 3] = 255;
      }
      octx.putImageData(img, 0, 0);
    }

    /* ---------- 绘制 ---------- */

    function pathInterior(k) {
      var x = BX + k * SC, y = BY + k * SC;
      var w = BE_W * SC - 2 * k * SC, h = BE_H * SC - 2 * k * SC;
      var rr = Math.max(3, (15 - k * 1.4) * SC);
      var fl = Math.max(2.5, (5 - k * 0.6) * SC);
      ctx.beginPath();
      ctx.moveTo(x - fl, y);
      ctx.lineTo(x, y + fl);
      ctx.lineTo(x, y + h - rr);
      ctx.quadraticCurveTo(x, y + h, x + rr, y + h);
      ctx.lineTo(x + w - rr, y + h);
      ctx.quadraticCurveTo(x + w, y + h, x + w, y + h - rr);
      ctx.lineTo(x + w, y + fl);
      ctx.lineTo(x + w + fl, y);
      ctx.closePath();
    }

    function pathLiquid() {
      var y = levelPx(), m = 2.1 * SC;
      var x0 = BX, x1 = BX + BE_W * SC;
      ctx.beginPath();
      ctx.moveTo(x0, y - m);
      ctx.quadraticCurveTo((x0 + x1) / 2, y + m * 0.55, x1, y - m);
      ctx.lineTo(x1, BY + BE_H * SC);
      ctx.lineTo(x0, BY + BE_H * SC);
      ctx.closePath();
    }

    function cellX(i) { return BX + ((i % NX) + 0.5) / NX * BE_W * SC; }

    function drawBench() {
      ctx.fillStyle = '#e9e5db';
      ctx.fillRect(0, 0, cvW(), cvH());
      var yb = BY + (BE_H + 2 * WALL) * SC;
      ctx.strokeStyle = 'rgba(28,26,23,0.10)';
      ctx.lineWidth = 1 * SC;
      ctx.beginPath();
      ctx.moveTo(0, yb); ctx.lineTo(cvW(), yb); ctx.stroke();
    }
    function cvW() { return canvas.width / (dpr * SC); }
    function cvH() { return canvas.height / (dpr * SC); }

    function drawShadow() {
      var cy = BY + (BE_H + WALL * 1.6) * SC;
      var g = ctx.createRadialGradient(cellX(27), cy, 2, cellX(27), cy, 96 * SC);
      g.addColorStop(0, 'rgba(28,26,23,0.20)');
      g.addColorStop(1, 'rgba(28,26,23,0)');
      ctx.fillStyle = g;
      ctx.save();
      ctx.translate(cellX(27), cy); ctx.scale(1, 0.13); ctx.translate(-cellX(27), -cy);
      ctx.beginPath(); ctx.arc(cellX(27), cy, 96 * SC, 0, 6.284); ctx.fill();
      ctx.restore();
    }

    function drawBurette() {
      var rk = pour ? pour.key : phaseKey;
      var rd = S.reagentByKey[rk];
      if (!rd) return;
      var pos = pour ? pour.pos : (NX / 2 - 9 + (S.reagents.indexOf(rd) % 2) * 17);
      var tx = cellX(Math.round(pos));
      var tipY = BY - 30 * SC;
      var w = 7.5 * SC;
      var topY = -10 * SC;
      var h = tipY - topY - 18 * SC;
      var fill = pour ? (1 - pour.remain / pour.total) : 1;
      if (pour && fill > 0) {
        ctx.save();
        ctx.beginPath(); ctx.rect(tx - w + 1.2 * SC, topY, w * 2 - 2.4 * SC, h); ctx.clip();
        ctx.fillStyle = 'rgba(190,206,214,0.5)';
        ctx.fillRect(tx - w, topY, w * 2, h * (1 - fill));
        ctx.restore();
      }
      var g = ctx.createLinearGradient(tx - w, 0, tx + w, 0);
      g.addColorStop(0, 'rgba(120,132,138,0.42)');
      g.addColorStop(0.22, 'rgba(255,255,255,0.42)');
      g.addColorStop(0.5, 'rgba(236,240,241,0.20)');
      g.addColorStop(0.82, 'rgba(255,255,255,0.30)');
      g.addColorStop(1, 'rgba(120,132,138,0.45)');
      ctx.fillStyle = g;
      ctx.fillRect(tx - w, topY, w * 2, h);
      ctx.strokeStyle = 'rgba(28,26,23,0.42)'; ctx.lineWidth = 1 * SC;
      ctx.strokeRect(tx - w, topY, w * 2, h);
      ctx.strokeStyle = 'rgba(28,26,23,0.30)';
      for (var m = 0; m <= 8; m++) {
        var yy = topY + 12 * SC + m * 8 * SC;
        if (yy > topY + h - 8 * SC) break;
        var len = (m % 4 === 0 ? 6.5 : 3.2) * SC;
        ctx.beginPath(); ctx.moveTo(tx + w - len, yy); ctx.lineTo(tx + w, yy); ctx.stroke();
      }
      ctx.fillStyle = '#3a352f';
      ctx.fillRect(tx - w - 2.5 * SC, tipY - 14 * SC, w * 2 + 5 * SC, 5 * SC);
      ctx.fillStyle = '#57504a';
      ctx.beginPath(); ctx.arc(tx, tipY - 11.5 * SC, 2.6 * SC, 0, 6.284); ctx.fill();
      ctx.beginPath();
      ctx.moveTo(tx - 3.4 * SC, tipY - 8 * SC);
      ctx.lineTo(tx - 1.3 * SC, tipY + 1.5 * SC);
      ctx.lineTo(tx + 1.3 * SC, tipY + 1.5 * SC);
      ctx.lineTo(tx + 3.4 * SC, tipY - 8 * SC);
      ctx.closePath();
      ctx.fillStyle = 'rgba(200,214,220,0.55)'; ctx.fill();
      ctx.strokeStyle = 'rgba(28,26,23,0.38)'; ctx.stroke();

      var lbl = rd.formula;
      ctx.textBaseline = 'middle';
      ctx.font = (11.5 * SC) + 'px ' + FONT;
      var f = fitText(lbl, tx + 14 * SC, 'left');
      var lw = f.w, lx = f.x;
      if (lx + lw + 4 * SC > tx) lx = fitText(lbl, tx - 14 * SC, 'right').x;
      ctx.textAlign = 'left';
      ctx.fillStyle = '#4d4840';
      ctx.fillText(lbl, lx, topY + h * 0.34);
      ctx.strokeStyle = 'rgba(28,26,23,0.35)'; ctx.lineWidth = 1 * SC;
      ctx.beginPath();
      if (lx > tx) { ctx.moveTo(tx + 8 * SC, topY + h * 0.34); ctx.lineTo(lx - 3 * SC, topY + h * 0.34); }
      else { ctx.moveTo(tx - 8 * SC, topY + h * 0.34); ctx.lineTo(lx + lw + 3 * SC, topY + h * 0.34); }
      ctx.stroke();
    }

    function drawStream() {
      if (!pour) return;
      var x = cellX(Math.round(pour.pos));
      var y0 = BY - 27 * SC;
      var y1 = levelPx();
      var len = Math.max(4, y1 - y0);
      var w = (1.5 + Math.min(2.4, len * 0.011)) * SC;
      var cr = 210, cg = 210, cb = 210;
      var ins = S.reagentByKey[pour.key].injects;
      for (var k in ins) {
        var sp = S.byKey[k];
        if (!sp.absorb) continue;
        cr = Math.exp(-1.7 * sp.absorb[0] * pour.conc) * 255;
        cg = Math.exp(-1.7 * sp.absorb[1] * pour.conc) * 255;
        cb = Math.exp(-1.7 * sp.absorb[2] * pour.conc) * 255;
        break;
      }
      var g = ctx.createLinearGradient(x - w, 0, x + w, 0);
      g.addColorStop(0, 'rgba(' + (cr | 0) + ',' + (cg | 0) + ',' + (cb | 0) + ',0.25)');
      g.addColorStop(0.38, 'rgba(' + (cr | 0) + ',' + (cg | 0) + ',' + (cb | 0) + ',0.62)');
      g.addColorStop(0.58, 'rgba(255,255,255,0.62)');
      g.addColorStop(1, 'rgba(' + (cr | 0) + ',' + (cg | 0) + ',' + (cb | 0) + ',0.3)');
      ctx.fillStyle = g;
      ctx.beginPath();
      var n = 16;
      for (var i = 0; i <= n; i++) {
        var f2 = i / n;
        var wob = Math.sin(f2 * 8 - simT * 12) * 0.45 * SC;
        ctx.lineTo(x - w * (0.8 + 0.5 * f2) + wob, y0 + len * f2);
      }
      for (i = n; i >= 0; i--) {
        var f3 = i / n;
        var wob2 = Math.sin(f3 * 8 - simT * 12 + 0.9) * 0.5 * SC;
        ctx.lineTo(x + w * (0.8 + 0.6 * f3) + wob2, y0 + len * f3);
      }
      ctx.closePath(); ctx.fill();
      ctx.strokeStyle = 'rgba(255,255,255,0.32)';
      ctx.lineWidth = 1 * SC;
      ctx.beginPath();
      ctx.ellipse(x, y1 - 0.5 * SC, w * 2.4, 1.3 * SC, 0, 0, 6.284);
      ctx.stroke();
    }

    function drawGlass() {
      var x = BX, y = BY, w = BE_W * SC, h = BE_H * SC;
      var xo = x - WALL * SC, yo = y - WALL * SC, wo = w + 2 * WALL * SC, ho = h + 2 * WALL * SC;
      ctx.save();
      ctx.beginPath();
      var rr = 17 * SC, fl = 5 * SC;
      ctx.moveTo(xo - fl, yo);
      ctx.lineTo(xo, yo + fl);
      ctx.lineTo(xo, yo + ho - rr);
      ctx.quadraticCurveTo(xo, yo + ho, xo + rr, yo + ho);
      ctx.lineTo(xo + wo - rr, yo + ho);
      ctx.quadraticCurveTo(xo + wo, yo + ho, xo + wo, yo + ho - rr);
      ctx.lineTo(xo + wo, yo + fl);
      ctx.lineTo(xo + wo + fl, yo);
      ctx.closePath();
      var g = ctx.createLinearGradient(xo, 0, xo + wo, 0);
      g.addColorStop(0, 'rgba(120,130,136,0.30)');
      g.addColorStop(0.05, 'rgba(255,255,255,0.34)');
      g.addColorStop(0.5, 'rgba(255,255,255,0.05)');
      g.addColorStop(0.95, 'rgba(255,255,255,0.26)');
      g.addColorStop(1, 'rgba(120,130,136,0.32)');
      ctx.fillStyle = g; ctx.fill();
      ctx.strokeStyle = 'rgba(28,26,23,0.5)'; ctx.lineWidth = 1.2 * SC; ctx.stroke();
      ctx.restore();
      ctx.save();
      pathInterior(0.8); ctx.clip();
      var gl = ctx.createLinearGradient(x, 0, x + 16 * SC, 0);
      gl.addColorStop(0, 'rgba(255,255,255,0.42)');
      gl.addColorStop(1, 'rgba(255,255,255,0)');
      ctx.fillStyle = gl; ctx.fillRect(x, y, 18 * SC, h);
      var gr = ctx.createLinearGradient(x + w - 7 * SC, 0, x + w, 0);
      gr.addColorStop(0, 'rgba(255,255,255,0)');
      gr.addColorStop(1, 'rgba(255,255,255,0.30)');
      ctx.fillStyle = gr; ctx.fillRect(x + w - 8 * SC, y, 8 * SC, h);
      ctx.restore();
      ctx.strokeStyle = 'rgba(28,26,23,0.45)'; ctx.lineWidth = 1.1 * SC;
      ctx.beginPath();
      ctx.ellipse((x + w / 2), yo + fl * 0.4, w / 2 + fl, 4.6 * SC, 0, 0, 6.284);
      ctx.stroke();
      ctx.strokeStyle = 'rgba(255,255,255,0.55)';
      ctx.beginPath();
      ctx.ellipse((x + w / 2), yo + fl * 0.4, w / 2 + fl, 4.6 * SC, 0, 0.35, 2.79);
      ctx.stroke();
    }

    function drawGraduations() {
      var x = BX, w = BE_W * SC;
      ctx.font = (9.5 * SC) + 'px ' + FONT2;
      ctx.textBaseline = 'middle';
      ctx.textAlign = 'left';
      for (var v = 10; v <= 130; v += 10) {
        var rows = v * ROWS_PER_ML;
        if (rows > NY) break;
        var y = BY + (NY - rows) / NY * BE_H * SC;
        var big = (v % 50 === 0);
        var len = (big ? 13 : 6) * SC;
        ctx.lineWidth = 1 * SC;
        ctx.strokeStyle = 'rgba(255,255,255,0.45)';
        ctx.beginPath(); ctx.moveTo(x + 1 * SC, y); ctx.lineTo(x + len + 1 * SC, y); ctx.stroke();
        ctx.strokeStyle = big ? 'rgba(28,26,23,0.38)' : 'rgba(28,26,23,0.20)';
        ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + len, y); ctx.stroke();
        if (big) {
          ctx.shadowColor = 'rgba(255,255,255,0.75)'; ctx.shadowBlur = 2 * SC;
          ctx.fillStyle = 'rgba(28,26,23,0.45)';
          ctx.fillText(String(v), x + len + 4 * SC, y);
          ctx.shadowBlur = 0;
        }
      }
      ctx.save();
      ctx.translate(BX - 12 * SC, BY + 26 * SC);
      ctx.rotate(Math.PI / 2);
      ctx.textAlign = 'left';
      ctx.font = (10 * SC) + 'px ' + FONT2;
      ctx.fillStyle = 'rgba(28,26,23,0.40)';
      ctx.fillText('mL', 0, 0);
      ctx.restore();
    }

    function fitText(t, wantX, align) {
      var w = ctx.measureText(t).width;
      var maxX = cvW() - 4 * SC;
      var x = align === 'right' ? wantX - w : wantX;
      if (x < 4 * SC) x = 4 * SC;
      if (x + w > maxX) x = maxX - w;
      if (x < 4 * SC) x = 4 * SC;
      return { x: x, w: w };
    }

    function leader(dx, dy, tx, ty, t, side) {
      var dir = side === 'left' ? -1 : 1;
      var f = fitText(t, tx + dir * 3 * SC, side === 'left' ? 'right' : 'left');
      var ex = f.x;
      ctx.beginPath();
      ctx.moveTo(dx, dy);
      ctx.lineTo(dx + dir * 11 * SC, ty);
      ctx.lineTo(ex + (side === 'left' ? -3 * SC : 3 * SC), ty);
      ctx.stroke();
      ctx.fillStyle = 'rgba(28,26,23,0.7)';
      ctx.beginPath(); ctx.arc(dx, dy, 1.8 * SC, 0, 6.284); ctx.fill();
      ctx.textAlign = side === 'left' ? 'right' : 'left';
      ctx.shadowColor = 'rgba(255,255,255,0.92)'; ctx.shadowBlur = 3 * SC;
      ctx.fillText(t, ex, ty);
      ctx.shadowBlur = 0;
    }

    function drawLabels() {
      if (labelsOn <= 0.01) return;
      ctx.save();
      ctx.globalAlpha = Math.min(1, labelsOn);
      ctx.font = (11.5 * SC) + 'px ' + FONT;
      ctx.textBaseline = 'middle';
      var top = levelPx(), bot = BY + BE_H * SC;
      var rx = BX + BE_W * SC + 14 * SC;

      // 固相：取含量最多的那个命名
      var best = null, bestM = 2e-5;
      for (var i = 0; i < S.solids.length; i++) {
        var sp = S.solids[i];
        if ((sp.mol || 0) > bestM) { bestM = sp.mol; best = sp; }
      }
      if (best) {
        var px = BX + BE_W * SC * 0.42, py = bot - (bot - top) * 0.07;
        ctx.strokeStyle = 'rgba(28,26,23,0.5)'; ctx.lineWidth = 1 * SC;
        leader(px, py - 2 * SC, rx + 8 * SC, py - 26 * SC, best.flag || best.name, 'right');
      }

      // 溶液：只认真正让液体显色的物种（有吸收系数且浓度可观），
      // 不能按"含量最多"选——NH₃ 数量最多但无色
      var aq = null, aqM = 0;
      for (i = 0; i < S.solutes.length; i++) {
        var s2 = S.solutes[i];
        var A = s2.absorb;
        if (!A || A[0] + A[1] + A[2] <= 0) continue;
        var v = s2.mol || 0;
        var dark = v * (A[0] + A[1] + A[2]);   // 按显色能力加权
        if (dark > aqM) { aqM = dark; aq = s2; }
      }
      if (aq) {
        var over = S.ratio && isOverBase();
        var sx = BX + BE_W * SC * 0.58, sy = top + (bot - top) * 0.28;
        ctx.strokeStyle = 'rgba(28,26,23,0.42)';
        leader(sx, sy, rx + 8 * SC, sy - 20 * SC,
               '上清液：' + (aq.flag || aq.name) + (over ? '' : ''), 'right');
      }
      ctx.restore();
    }

    function isOverBase() {
      if (!S.base) return false;
      var b = S.base, m = S.metal;
      return (S.byKey[b].mol || 0) * S.ratio.stoich > (S.byKey[m].mol || 0);
    }

    function drawBurner() {
      if (!heating) return;
      var cx = cellX(NX / 2 - 1);
      var baseY = BY + (BE_H + 2 * WALL) * SC + 1 * SC;
      var lift = (temp - T_ROOM) / (T_BOIL - T_ROOM);
      ctx.save();
      var fh = (10 + 20 * lift) * SC;
      var wob = Math.sin(simT * 13) * 1.1 * SC + Math.sin(simT * 7.3) * 0.7 * SC;
      var f1 = ctx.createLinearGradient(cx, baseY - fh, cx, baseY);
      f1.addColorStop(0, 'rgba(120,150,180,0.42)');
      f1.addColorStop(1, 'rgba(120,150,180,0.06)');
      ctx.fillStyle = f1;
      ctx.beginPath();
      ctx.moveTo(cx - 5.5 * SC, baseY);
      ctx.quadraticCurveTo(cx - 6.5 * SC + wob, baseY - fh * 0.55, cx + wob * 0.7, baseY - fh);
      ctx.quadraticCurveTo(cx + 6.5 * SC + wob, baseY - fh * 0.55, cx + 5.5 * SC, baseY);
      ctx.closePath(); ctx.fill();
      var fh2 = fh * 0.6;
      var f2 = ctx.createLinearGradient(cx, baseY - fh2, cx, baseY);
      f2.addColorStop(0, 'rgba(198,132,86,0.60)');
      f2.addColorStop(1, 'rgba(150,96,60,0.14)');
      ctx.fillStyle = f2;
      ctx.beginPath();
      ctx.moveTo(cx - 3 * SC, baseY);
      ctx.quadraticCurveTo(cx - 3.4 * SC + wob * 0.6, baseY - fh2 * 0.55, cx + wob * 0.5, baseY - fh2);
      ctx.quadraticCurveTo(cx + 3.4 * SC + wob * 0.6, baseY - fh2 * 0.55, cx + 3 * SC, baseY);
      ctx.closePath(); ctx.fill();
      ctx.fillStyle = 'rgba(96,90,82,0.85)';
      ctx.fillRect(cx - 7 * SC, baseY, 14 * SC, 2.4 * SC);
      ctx.beginPath();
      ctx.moveTo(cx - 4.4 * SC, baseY + 2.4 * SC);
      ctx.lineTo(cx - 3.2 * SC, baseY + 13 * SC);
      ctx.lineTo(cx + 3.2 * SC, baseY + 13 * SC);
      ctx.lineTo(cx + 4.4 * SC, baseY + 2.4 * SC);
      ctx.closePath();
      ctx.fillStyle = 'rgba(74,70,64,0.9)';
      ctx.fill();
      ctx.strokeStyle = 'rgba(28,26,23,0.55)'; ctx.lineWidth = 1 * SC; ctx.stroke();
      ctx.fillStyle = 'rgba(40,36,32,0.9)';
      ctx.fillRect(cx - 1.4 * SC, baseY - 1.5 * SC, 2.8 * SC, 2 * SC);
      ctx.restore();
    }

    function drawBubbles() {
      if (boil <= 0.02 || V < 1) return;
      var x0 = BX + 5 * SC, x1 = BX + (BE_W - 5) * SC;
      var yTop = levelPx(), yBot = BY + BE_H * SC;
      ctx.save();
      ctx.strokeStyle = 'rgba(255,255,255,0.55)';
      ctx.lineWidth = 1 * SC;
      for (var k = 0; k < 16; k++) {
        var ph = ((simT * (0.7 + (k % 5) * 0.16) + k * 0.137) % 1);
        var bx = x0 + ((k * 7919) % 1000) / 1000 * (x1 - x0);
        var by = yBot - ph * (yBot - yTop);
        if (by < yTop) continue;
        var r = (0.7 + (k % 3) * 0.45) * SC * (0.5 + 0.5 * Math.sin(ph * Math.PI));
        ctx.globalAlpha = 0.25 + 0.4 * Math.sin(ph * Math.PI) * boil;
        ctx.beginPath(); ctx.arc(bx, by, r, 0, 6.284); ctx.stroke();
      }
      ctx.restore();
    }

    function draw() {
      ctx.setTransform(dpr * SC, 0, 0, dpr * SC, 0, 0);
      ctx.clearRect(0, 0, cvW(), cvH());
      drawBench();
      drawShadow();
      drawBurette();
      ctx.save();
      pathInterior(0); ctx.clip();
      pathLiquid(); ctx.clip();
      buildImage();
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(off, BX, BY, BE_W * SC, BE_H * SC);
      var w = BE_W * SC;
      var gl = ctx.createLinearGradient(BX, 0, BX + 13 * SC, 0);
      gl.addColorStop(0, 'rgba(20,26,32,0.20)'); gl.addColorStop(1, 'rgba(20,26,32,0)');
      ctx.fillStyle = gl; ctx.fillRect(BX, levelPx() - 4, 13 * SC, BE_H * SC);
      var gr = ctx.createLinearGradient(BX + w - 8 * SC, 0, BX + w, 0);
      gr.addColorStop(0, 'rgba(20,26,32,0)'); gr.addColorStop(1, 'rgba(20,26,32,0.14)');
      ctx.fillStyle = gr; ctx.fillRect(BX + w - 8 * SC, levelPx() - 4, 8 * SC, BE_H * SC);
      drawBubbles();
      ctx.restore();
      if (V > 0.4) {
        ctx.strokeStyle = 'rgba(255,255,255,0.5)'; ctx.lineWidth = 1.1 * SC;
        ctx.beginPath();
        var y = levelPx(), m = 2.1 * SC;
        ctx.moveTo(BX, y - m);
        ctx.quadraticCurveTo(BX + w / 2, y + m * 0.55, BX + w, y - m);
        ctx.stroke();
      }
      drawGlass();
      drawGraduations();
      drawStream();
      drawBurner();
      drawLabels();
    }

    /* ---------- 尺寸 ---------- */

    function fitScale(fixed, span, box) {
      var disc = fixed * fixed + 4 * span * box;
      if (disc <= 0) return 1;
      return (-fixed + Math.sqrt(disc)) / (2 * span);
    }

    // 画布内的水平居中量：固定留白按各自尺寸缩放，杯子居中
    function fitOff(fixed, span) { return (fixed + span * SC) / 2; }

    function resize() {
      var el = document.querySelector('main');
      if (!el) return;
      var cw = el.clientWidth - 48, ch = el.clientHeight - 18;
      var sw = BE_W + 2 * WALL, sh = BE_H + 2 * WALL;
      SC = Math.min(1.5, fitScale(ML + MR, sw, cw - 6), fitScale(TOP + BURNER_H, sh, ch - 6));
      if (!(SC > 0.05)) SC = 0.05;
      dpr = Math.min(2, window.devicePixelRatio || 1);
      var tw = ML + sw * SC + MR, th = TOP + sh * SC + BURNER_H;
      canvas.style.width = (tw * SC) + 'px';
      canvas.style.height = (th * SC) + 'px';
      canvas.width = Math.round(tw * SC * dpr);
      canvas.height = Math.round(th * SC * dpr);
      BX = fitOff(ML, sw);       // 水平居中
      BY = TOP;                  // 滴定管在杯子上方，顶部留白不能动
    }

    /* ---------- 对外 ---------- */

    var api = {
      resize: resize,
      load: load,
      restart: restart,
      step: step,
      draw: draw,
      currentOrder: currentOrder,
      startPour: startPour,
      setHeating: function (v) { heating = !!v; },
      setSwirl: function (v) { swirl = v; },
      tickLabels: function (dt) { if (phase === 'done') labelsOn = Math.min(1, labelsOn + dt * 1.2); },
      // 相位机
      advance: function (dt) {
        if (phase === 'pour' && !pour) { phase = 'wait'; wait = 2.0; }
        else if (phase === 'wait') {
          wait -= dt;
          if (wait <= 0) {
            if (queue.length) { phaseKey = queue.shift(); startPour(phaseKey); phase = 'pour'; }
            else phase = 'done';
          }
        }
      },
      state: function () {
        return { V: V, temp: temp, heating: heating, phase: phase,
                 pouring: pour ? pour.key : null,
                 species: S.species, tally: tally,
                 ratio: null };
      },
      // 供界面读取
      plan: planned,
      S: function () { return S; },
      F: function () { return F; },
      get temp() { return temp; },
      get V() { return V; },
      get phase() { return phase; },
      setReagent: function (k, conc, vol) {
        var rd = S.reagentByKey[k];
        rd.conc = conc; rd.vol = vol;
      },
      // token 台账：某 token 当前被占住多少 / 加入总量 / 逸出多少
      ledger: function () {
        var out = {};
        for (var t = 0; t < S.tokenList.length; t++) {
          var tk = S.tokenList[t];
          var added = 0, held = 0;
          for (var r = 0; r < S.reagents.length; r++) {
            var rd = S.reagents[r];
            var n = rd.conc * rd.vol / 1000;
            for (var k in rd.injects) if (S.byKey[k].holds && S.byKey[k].holds[tk]) added += n * S.byKey[k].holds[tk];
          }
          for (var i = 0; i < S.species.length; i++) {
            var sp = S.species[i];
            if (sp.holds && sp.holds[tk]) held += (sp.mol || 0) * sp.holds[tk];
          }
          out[tk] = { added: added, held: held, escaped: added - held };
        }
        return out;
      }
    };
    return api;
  }

  root.ReactionEngine = { create: createEngine, GRID: { NX: NX, NY: NY } };
})(window);