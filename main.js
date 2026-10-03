/* CuSO4 / NaOH —— 浓度对反应产物的影响
   场量模型：每一格记录 Cu2+、OH-、Cu(OH)2 固体、[Cu(OH)4]2- 的量（mol·L-1），
   用半拉格朗日平流 + 扩散 + 局部反应推进。 */

(function () {
  'use strict';

  /* ============================ 常量 ============================ */

  var NX = 56, NY = 92, N = NX * NY;
  var ROWS_PER_ML = 0.62;                  // 每毫升液体占的格数
  var CV = 1 / (ROWS_PER_ML * NX);         // 每格体积 mL
  var MAXV = 132;              // 烧杯容量上限 mL
  var POUR_RATE = 17;          // 滴定速度 mL/s

  var OHHI = 1.7;              // 生成 [Cu(OH)4]2- 的 OH- 阈值 mol·L-1
  var PPT_MAX = 2.6;            // 絮状物的堆积密度 mol·L-1
  var PPT_NEAR = 0.30;          // 团块挤开的门限

  var PATH = 1.7;
  var ABS_CU = [1.45, 0.87, 0.23];
  var ABS_CU4 = [8.50, 4.00, 1.10];
  var PPTX = 1.7;
  var SCAT_L = [0.72, 0.86, 1.00];        // 稀薄悬浮的絮状物
  var SCAT_D = [0.42, 0.62, 0.93];        // 沉积浓密处

  /* ============================ 场 ============================ */

  var cu = new Float32Array(N), oh = new Float32Array(N),
      ppt = new Float32Array(N), cu4 = new Float32Array(N),
      tmp = new Float32Array(N), dv = new Float32Array(N),
      velX = new Float32Array(N), velY = new Float32Array(N), dv = new Float32Array(N);

  // 粗粒度的絮状噪声（按 2×2 分块，避免出现像素噪点）
  var noise = new Float32Array(N);
  (function () {
    var s = 20240117;
    function rnd() { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; }
    for (var y = 0; y < NY; y += 2) {
      for (var x = 0; x < NX; x += 2) {
        var v = (rnd() + rnd() + rnd() + rnd()) / 4 - 0.5;
        for (var j = 0; j < 2; j++) for (var i = 0; i < 2; i++) {
          if (y + j < NY && x + i < NX) noise[(y + j) * NX + x + i] = v;
        }
      }
    }
  })();

  /* ============================ 标量账目 ============================ */

  var V = 0, nCuTot = 0, nOHTot = 0, nP = 0, nCu4 = 0;
  var nCu2 = 0, nOHf = 0;                // 由 tally() 从场统计
  var LV = 0;                  // 液面高度（格）

  function wet(y) { return y >= 0 && y < NY && (NY - 1 - y) < LV; }

  /* ============================ 参数 ============================ */

  var params = { cCu: 1, vCu: 40, cNa: 2, vNa: 40, order: 'cu_first' };

  var PRESETS = {
    eq:          { cCu: 1,    vCu: 40, cNa: 2,    vNa: 40 },
    base_lack:   { cCu: 0.5,  vCu: 50, cNa: 0.5,  vNa: 25 },
    base_excess: { cCu: 1,    vCu: 40, cNa: 3,    vNa: 60 },
    concentrated:{ cCu: 1.6,  vCu: 30, cNa: 8,    vNa: 40 },
    dilute:      { cCu: 0.05, vCu: 40, cNa: 0.06, vNa: 40 }
  };

  /* ============================ 运行状态 ============================ */

  var pour = null;             // {r, remain, total, pos, conc}
  var pourPos = NX / 2 - 8;
  var jetT = 0, jetPos = pourPos;
  var swirl = 0;
  var simT = 0;
  var phase = 'pour', phaseR = 'cu', phaseT = 0, queue = [];
  var labelsOn = 0;

  function restart() {
    cu.fill(0); oh.fill(0); ppt.fill(0); cu4.fill(0);
    V = 0; nCuTot = 0; nOHTot = 0; nP = 0; nCu4 = 0; nCu2 = 0; nOHf = 0;
    pour = null; jetT = 0; swirl = 0; simT = 0; labelsOn = 0;
    var first = params.order === 'cu_first' ? 'cu' : 'na';
    queue = [first === 'cu' ? 'na' : 'cu'];
    phase = 'pour'; phaseR = first;
    startPour(first);
    updatePanel();
  }

  function startPour(r) {
    var conc = r === 'cu' ? params.cCu : params.cNa;
    var vol = r === 'cu' ? params.vCu : params.vNa;
    pourPos = r === 'cu' ? NX / 2 - 9 : NX / 2 + 8;
    pour = { r: r, remain: vol, total: vol, pos: pourPos, conc: conc };
  }

  /* ============================ 物理步进 ============================ */

  var wbuf = new Float32Array(96);
  function gauss(x, s) { return Math.exp(-(x * x) / (2 * s * s)); }

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

  // 由流函数 ψ = A·sin(πfx)·sin(πfy) 导出，天然无散
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

  function advect(dt, f, vx, vy) {
    var yTop = NY - Math.ceil(LV);            // 最上一行湿格，取样不越过液面
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

  function sample(f, x, y) {
    var x0 = Math.floor(x), y0 = Math.floor(y);
    if (x0 < 0) x0 = 0; else if (x0 > NX - 2) x0 = NX - 2;
    if (y0 < 0) y0 = 0; else if (y0 > NY - 2) y0 = NY - 2;
    var tx = x - x0, ty = y - y0;
    var a = f[y0 * NX + x0], b = f[y0 * NX + x0 + 1];
    var c = f[(y0 + 1) * NX + x0], d = f[(y0 + 1) * NX + x0 + 1];
    return (a + (b - a) * tx) + ((c + (d - c) * tx) - (a + (b - a) * tx)) * ty;
  }

  function diffuse(dt, f, D) {
    var k = Math.min(D, 0.2 / dt) * dt;
    for (var y = 0; y < NY; y++) {
      var rowWet = wet(y), upWet = wet(y - 1), dnWet = wet(y + 1);
      for (var x = 0; x < NX; x++) {
        if (!rowWet) { tmp[y * NX + x] = 0; continue; }
        var i = y * NX + x, c = f[i];
        var l = (x > 0) ? f[i - 1] : c;
        var r = (x < NX - 1) ? f[i + 1] : c;
        var u = upWet ? f[i - NX] : c;
        var d2 = dnWet ? f[i + NX] : c;
        tmp[i] = c + k * (l + r + u + d2 - 4 * c);
      }
    }
    f.set(tmp);
  }

  // 湍流松弛：让溶液最终趋于均匀，但短时间内保留羽状层次
  function relax(f, rate, dt) {
    var s = 0, n = 0;
    for (var y = 0; y < NY; y++) {
      if (!wet(y)) continue;
      var o = y * NX;
      for (var x = 0; x < NX; x++) { s += f[o + x]; n++; }
    }
    if (!n) return;
    var m = s / n, k = 1 - Math.exp(-rate * dt);
    for (var y2 = 0; y2 < NY; y2++) {
      if (!wet(y2)) continue;
      var o2 = y2 * NX;
      for (var x2 = 0; x2 < NX; x2++) f[o2 + x2] += (m - f[o2 + x2]) * k;
    }
  }

  function react(dt) {
    var kR = 1 - Math.exp(-20 * dt);
    var kD = 1 - Math.exp(-0.25 * dt);
    for (var y = 0; y < NY; y++) {
      if (!wet(y)) continue;
      for (var x = 0; x < NX; x++) {
        var i = y * NX + x;
        var a = oh[i] * 0.5;
        if (a > cu[i]) a = cu[i];
        a *= kR;
        if (a > 0) { cu[i] -= a; oh[i] -= 2 * a; ppt[i] += a; }
        var hi = oh[i] - OHHI - 2 * cu4[i];
        if (ppt[i] > 0 && hi > 0) {
          var b = hi * 0.5;
          if (b > ppt[i]) b = ppt[i];
          b *= kD;
          if (b > 0) { ppt[i] -= b; oh[i] -= 2 * b; cu4[i] += b; }
        }
      }
    }
    tally();
  }

  // 湿区求和与整体缩放（半拉格朗日平流不守恒，用它补回总量）
  function sumWet(f) {
    var s = 0;
    for (var y = 0; y < NY; y++) {
      if (!wet(y)) continue;
      var o = y * NX;
      for (var x = 0; x < NX; x++) s += f[o + x];
    }
    return s;
  }

  // 场是唯一真相：每步统计各物种总量
  function tally() {
    var c = 0, o = 0, p = 0, s4 = 0;
    for (var y = 0; y < NY; y++) {
      if (!wet(y)) continue;
      var b = y * NX;
      for (var x = 0; x < NX; x++) {
        var i = b + x;
        c += cu[i]; o += oh[i]; p += ppt[i]; s4 += cu4[i];
      }
    }
    var f = CV / 1000;
    nCu2 = c * f; nOHf = o * f; nP = p * f; nCu4 = s4 * f;
  }

  // 平流。液面之上不取样（否则溶质会被"漏"进干区），再做一次有界的总量修正
  function advectKeep(dt, f, vx, vy) {
    var n = 0;
    for (var y = 0; y < NY; y++) { if (wet(y)) n += NX; }
    if (!n) return;
    var s = sumWet(f);
    advect(dt, f, vx, vy);
    var s2 = sumWet(f);
    if (s > 1e-12 && s2 > 1e-12) {
      var k = s / s2;                       // 取样有界后每帧误差有限，直接修正即可
      for (var y2 = 0; y2 < NY; y2++) {
        if (!wet(y2)) continue;
        var o = y2 * NX;
        for (var x = 0; x < NX; x++) f[o + x] *= k;
      }
    }
    for (var y3 = 0; y3 < NY; y3++) {          // 干区不留残值
      if (wet(y3)) continue;
      var o3 = y3 * NX;
      for (var x3 = 0; x3 < NX; x3++) f[o3 + x3] = 0;
    }
  }

  // 沉淀不参与平流，改用守恒的通量：水平铺展 + 重力沉降 + 底部堆积
  function pack(dt) {
    var x, y, i;
    // 1) 水平铺展（相邻格对称交换，总量守恒）
    var kh = 0.5 * (1 - Math.exp(-(7 + swirl * 26) * dt));
    dv.fill(0);
    for (y = 0; y < NY; y++) {
      if (!wet(y)) continue;
      var row = y * NX;
      for (x = 0; x < NX - 1; x++) {
        var a = ppt[row + x], b = ppt[row + x + 1];
        if (a <= 0 && b <= 0) continue;
        var m = (a - b) * kh * 0.5;
        dv[row + x] -= m; dv[row + x + 1] += m;
      }
    }
    for (i = 0; i < N; i++) if (dv[i] !== 0) { var t = ppt[i] + dv[i]; ppt[i] = t > 0 ? t : 0; }

    // 2) 沉降：下方装满就停住，于是自然堆出密实层；摇动时重新扬起
    var f = 8.5 * (1 - swirl * 0.9) * dt;
    var up = swirl * 11 * dt;
    for (y = 0; y < NY; y++) {
      if (!wet(y)) continue;
      var below = wet(y + 1);
      var above = wet(y - 1);
      var base = y * NX;
      for (x = 0; x < NX; x++) {
        i = base + x;
        if (above && up > 0) {
          var u = ppt[i] * up;
          if (u > 0) { ppt[i] -= u; ppt[i - NX] += u; }
        }
        if (ppt[i] <= 0) continue;
        if (below) {
          var mv = ppt[i] * f;
          var room = PPT_MAX - ppt[i + NX];
          if (room <= 0) continue;
          if (mv > room) mv = room;
          if (mv > 0) { ppt[i] -= mv; ppt[i + NX] += mv; }
          continue;
        }
        // 杯底：向左右分摊
        var mv2 = ppt[i] * f;
        if (mv2 <= 0) continue;
        var nl = (x > 0 ? 1 : 0) + (x < NX - 1 ? 1 : 0);
        if (!nl) continue;
        var share = mv2 / nl;
        ppt[i] -= mv2;
        if (x > 0) ppt[i - 1] += share;
        if (x < NX - 1) ppt[i + 1] += share;
        ppt[i] += mv2 - share * nl;
      }
    }
  }

  function step(dt) {
    var r = pour, mol = 0;
    if (r) {                       // 先加体积
      var room = MAXV - V;
      if (room <= 0) { pour = null; r = null; }
      else {
        var add = Math.min(r.remain, POUR_RATE * dt, room);
        r.remain -= add;
        V += add;
        mol = add * r.conc / 1000;
        if (r.r === 'cu') nCuTot += mol; else nOHTot += mol;
        if (r.remain <= 1e-6) pour = null;
      }
    }
    LV = Math.min(NY, V * ROWS_PER_ML);
    if (r && mol > 0) {            // 液面定下来之后再投放溶质
      var row = Math.min(NY - 0.6, V * ROWS_PER_ML);
      inject(r.pos, NY - row, mol, r.r === 'cu' ? cu : oh);
      jetT = 0.55; jetPos = r.pos;
    }
    if (!r && LV < 0.3) { jetT = 0; swirl *= Math.pow(0.22, dt); tally(); return; }
    simT += dt;
    buildVel();
    advectKeep(dt, cu, velX, velY);
    advectKeep(dt, oh, velX, velY);
    advectKeep(dt, cu4, velX, velY);
    var mix = (pour || jetT > 0) ? 2.2 : 0.9;
    diffuse(dt / 2, oh, 11); diffuse(dt / 2, oh, 11);
    diffuse(dt / 2, cu, 7); diffuse(dt / 2, cu, 7);
    diffuse(dt / 2, cu4, 8);
    react(dt);
    relax(cu, mix, dt); relax(oh, mix, dt); relax(cu4, mix, dt);
    pack(dt);
    tally();
    jetT = Math.max(0, jetT - dt * 2.2);
    swirl *= Math.pow(0.22, dt);
  }

  /* ============================ 画布 ============================ */

  var cv = document.getElementById('scene');
  var ctx = cv.getContext('2d');
  var BE_W = 176, BE_H = 290, WALL = 7;
  var ML = 56, MR = 126, TOP = 100;   // 画布留白（布局单位）
  var BX = 0, BY = 0, SC = 1, dpr = 1;
  var off = document.createElement('canvas');
  off.width = NX; off.height = NY;
  var octx = off.getContext('2d');
  var img = octx.createImageData(NX, NY);
  var BG = [233, 229, 219];

  // 画布变换是 setTransform(SC)，所以画布自身的布局尺寸里，杯子要按 SC 放大后再算，
// 求 SC 时 tw、th 都含 SC，是个一元二次方程。
  function fitScale(fixed, span, box) {
    var disc = fixed * fixed + 4 * span * box;
    if (disc <= 0) return 1;
    return (-fixed + Math.sqrt(disc)) / (2 * span);
  }

  function resize() {
    var el = document.querySelector('main');
    // clientWidth/Height 含 padding，要减掉才能得到真正能放画布的盒子
    var cw = el.clientWidth - 48, ch = el.clientHeight - 18;
    var sw = BE_W + 2 * WALL, sh = BE_H + 2 * WALL;
    SC = Math.min(1.5, fitScale(ML + MR, sw, cw - 6), fitScale(TOP + 22, sh, ch - 6));
    if (!(SC > 0.05)) SC = 0.05;
    dpr = Math.min(2, window.devicePixelRatio || 1);
    var tw = ML + sw * SC + MR, th = TOP + sh * SC + 22;
    cv.style.width = (tw * SC) + 'px';
    cv.style.height = (th * SC) + 'px';
    cv.width = Math.round(tw * SC * dpr);
    cv.height = Math.round(th * SC * dpr);
    BX = ML;              // 布局单位；setTransform 会再乘一次 SC，这里不能预先乘
    BY = TOP;
  }

  function levelPx() { return BY + (NY - Math.min(NY, V * ROWS_PER_ML)) / NY * BE_H * SC; }

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

  function buildImage() {
    var D = img.data;
    for (var i = 0, p = 0; i < N; i++, p += 4) {
      var a = PATH * (ABS_CU[0] * cu[i] + ABS_CU4[0] * cu4[i]);
      var b = PATH * (ABS_CU[1] * cu[i] + ABS_CU4[1] * cu4[i]);
      var c = PATH * (ABS_CU[2] * cu[i] + ABS_CU4[2] * cu4[i]);
      var tr = Math.exp(-a), tg = Math.exp(-b), tbb = Math.exp(-c);
      var pz = ppt[i];
      if (pz > 0) {
        var t = pz > PPT_NEAR ? Math.min(1, (pz - PPT_NEAR) / (PPT_MAX - PPT_NEAR)) : 0;
        var nz = 1 + noise[i] * 0.26 * Math.min(1, pz * 2.5);
        var sr = (SCAT_L[0] + (SCAT_D[0] - SCAT_L[0]) * t) * nz;
        var sg = (SCAT_L[1] + (SCAT_D[1] - SCAT_L[1]) * t) * nz;
        var sb = (SCAT_L[2] + (SCAT_D[2] - SCAT_L[2]) * t) * nz;
        var al = 1 - Math.exp(-pz * PPTX);
        var w = 1 - al;
        D[p] = (sr * al * 255 + tr * BG[0] * w) | 0;
        D[p + 1] = (sg * al * 255 + tg * BG[1] * w) | 0;
        D[p + 2] = (sb * al * 255 + tbb * BG[2] * w) | 0;
      } else {
        D[p] = (tr * BG[0]) | 0;
        D[p + 1] = (tg * BG[1]) | 0;
        D[p + 2] = (tbb * BG[2]) | 0;
      }
      D[p + 3] = 255;
    }
    octx.putImageData(img, 0, 0);
  }

  function cellX(i) { return BX + ((i % NX) + 0.5) / NX * BE_W * SC; }

  function drawBench() {
    ctx.fillStyle = '#e9e5db';
    ctx.fillRect(0, 0, cv.width, cv.height);
    var yb = BY + (BE_H + 2 * WALL) * SC;
    ctx.strokeStyle = 'rgba(28,26,23,0.10)';
    ctx.lineWidth = 1 * SC;
    ctx.beginPath();
    ctx.moveTo(0, yb); ctx.lineTo(cv.width, yb); ctx.stroke();
  }

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

  /* 滴定管 */
  function drawBurette() {
    var tx = cellX(pourPos);
    var tipY = BY - 30 * SC;
    var w = 7.5 * SC;
    var topY = -10 * SC;
    var h = tipY - topY - 18 * SC;
    var fill = pour ? (1 - pour.remain / pour.total) : (doneBurette ? 1 : 0);
    if (fill > 0) {
      ctx.save();
      ctx.beginPath(); ctx.rect(tx - w + 1.2 * SC, topY, w * 2 - 2.4 * SC, h); ctx.clip();
      ctx.fillStyle = 'rgba(190,206,214,0.5)';
      ctx.fillRect(tx - w, topY, w * 2, h * (1 - fill));
      ctx.restore();
    }
    // 管身
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
    // 刻度
    ctx.strokeStyle = 'rgba(28,26,23,0.30)';
    for (var m = 0; m <= 8; m++) {
      var yy = topY + 12 * SC + m * 8 * SC;
      if (yy > topY + h - 8 * SC) break;
      var len = (m % 4 === 0 ? 6.5 : 3.2) * SC;
      ctx.beginPath(); ctx.moveTo(tx + w - len, yy); ctx.lineTo(tx + w, yy); ctx.stroke();
    }
    // 活塞
    ctx.fillStyle = '#3a352f';
    ctx.fillRect(tx - w - 2.5 * SC, tipY - 14 * SC, w * 2 + 5 * SC, 5 * SC);
    ctx.fillStyle = '#57504a';
    ctx.beginPath();
    ctx.arc(tx, tipY - 11.5 * SC, 2.6 * SC, 0, 6.284); ctx.fill();
    // 玻璃尖嘴
    ctx.beginPath();
    ctx.moveTo(tx - 3.4 * SC, tipY - 8 * SC);
    ctx.lineTo(tx - 1.3 * SC, tipY + 1.5 * SC);
    ctx.lineTo(tx + 1.3 * SC, tipY + 1.5 * SC);
    ctx.lineTo(tx + 3.4 * SC, tipY - 8 * SC);
    ctx.closePath();
    ctx.fillStyle = 'rgba(200,214,220,0.55)'; ctx.fill();
    ctx.strokeStyle = 'rgba(28,26,23,0.38)'; ctx.stroke();

    // 标签
    var r = pour ? pour.r : phaseR;
    var nm = r === 'cu' ? 'CuSO₄' : 'NaOH';
    var cc = r === 'cu' ? params.cCu : params.cNa;
    var lbl = nm + '  ' + sig(cc, 2) + ' mol·L⁻¹';
    ctx.textBaseline = 'middle';
    ctx.font = (11.5 * SC) + 'px ' + FONT;
    var ly = topY + h * 0.34;
    var f = fitText(lbl, tx + 14 * SC, 'left');
    var lw = f.w;
    var lx = f.x;
    if (lx + lw + 4 * SC > tx) {                 // 右边放不下就换到左边
      lx = fitText(lbl, tx - 14 * SC, 'right').x;
    }
    ctx.textAlign = 'left';
    ctx.fillStyle = '#4d4840';
    ctx.fillText(lbl, lx, ly);
    ctx.strokeStyle = 'rgba(28,26,23,0.35)'; ctx.lineWidth = 1 * SC;
    ctx.beginPath();
    if (lx > tx) { ctx.moveTo(tx + 8 * SC, ly); ctx.lineTo(lx - 3 * SC, ly); }
    else { ctx.moveTo(tx - 8 * SC, ly); ctx.lineTo(lx + lw + 3 * SC, ly); }
    ctx.stroke();
  }

  function drawStream() {
    if (!pour) return;
    var x = cellX(pourPos);
    var y0 = BY - 27 * SC;
    var y1 = levelPx();
    var len = Math.max(4, y1 - y0);
    var w = (1.5 + Math.min(2.4, len * 0.011)) * SC;
    // 液流本身带试剂的颜色
    var cc = pour.conc * (pour.r === 'cu' ? 1 : 0);
    var cr = Math.exp(-PATH * ABS_CU[0] * cc) * 255;
    var cg = Math.exp(-PATH * ABS_CU[1] * cc) * 255;
    var cb = Math.exp(-PATH * ABS_CU[2] * cc) * 255;
    var g = ctx.createLinearGradient(x - w, 0, x + w, 0);
    g.addColorStop(0, 'rgba(' + (cr | 0) + ',' + (cg | 0) + ',' + (cb | 0) + ',0.25)');
    g.addColorStop(0.38, 'rgba(' + (cr | 0) + ',' + (cg | 0) + ',' + (cb | 0) + ',0.62)');
    g.addColorStop(0.58, 'rgba(255,255,255,0.62)');
    g.addColorStop(1, 'rgba(' + (cr | 0) + ',' + (cg | 0) + ',' + (cb | 0) + ',0.3)');
    ctx.fillStyle = g;
    ctx.beginPath();
    var n = 16;
    for (var k = 0; k <= n; k++) {
      var f = k / n;
      var wob = Math.sin(f * 8 - simT * 12) * 0.45 * SC;
      ctx.lineTo(x - w * (0.8 + 0.5 * f) + wob, y0 + len * f);
    }
    for (var k2 = n; k2 >= 0; k2--) {
      var f2 = k2 / n;
      var wob2 = Math.sin(f2 * 8 - simT * 12 + 0.9) * 0.5 * SC;
      ctx.lineTo(x + w * (0.8 + 0.6 * f2) + wob2, y0 + len * f2);
    }
    ctx.closePath(); ctx.fill();
    // 入水处的波纹
    ctx.strokeStyle = 'rgba(255,255,255,0.32)';
    ctx.lineWidth = 1 * SC;
    ctx.beginPath();
    ctx.ellipse(x, y1 - 0.5 * SC, w * 2.4, 1.3 * SC, 0, 0, 6.284);
    ctx.stroke();
  }

  function drawGlass() {
    var x = BX, y = BY, w = BE_W * SC, h = BE_H * SC;
    var xo = x - WALL * SC, yo = y - WALL * SC, wo = w + 2 * WALL * SC, ho = h + 2 * WALL * SC;
    // 玻璃体
    ctx.save();
    ctx.beginPath();
    var rr = (17) * SC, fl = 5 * SC;
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
    // 内壁反光
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
    // 口沿
    ctx.strokeStyle = 'rgba(28,26,23,0.45)'; ctx.lineWidth = 1.1 * SC;
    ctx.beginPath();
    ctx.ellipse((x + w / 2), yo + fl * 0.4, w / 2 + fl, 4.6 * SC, 0, 0, 6.284);
    ctx.stroke();
    ctx.strokeStyle = 'rgba(255,255,255,0.55)';
    ctx.beginPath();
    ctx.ellipse((x + w / 2), yo + fl * 0.4, w / 2 + fl, 4.6 * SC, 0, 0.35, 2.79);
    ctx.stroke();
  }

  // 刻度印在杯壁上（前视时压在液体之上），所以画在杯内左侧
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

  var doneBurette = false;

  // 文字排布：先量宽度，再夹回画布内，绝不越界
  function fitText(t, wantX, align) {
    var w = ctx.measureText(t).width;
    var maxX = cv.width / (dpr * SC) - 4 * SC;
    var x = align === 'right' ? wantX - w : wantX;
    if (x < 4 * SC) x = 4 * SC;
    if (x + w > maxX) x = maxX - w;
    if (x < 4 * SC) x = 4 * SC;
    return { x: x, w: w };
  }

  // 一条引线 + 端点圆点 + 文字；文字左右可挤，画线跟着文字走
  function leader(dx, dy, tx, ty, t, side) {
    var dir = side === 'left' ? -1 : 1;
    var f = fitText(t, tx + dir * 3 * SC, side === 'left' ? 'right' : 'left');
    var ex = side === 'left' ? f.x : f.x;
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
    var a = Math.min(1, labelsOn);
    var ratio = nCuTot > 1e-9 ? nOHTot / nCuTot : null;
    ctx.save();
    ctx.globalAlpha = a;
    ctx.font = (11.5 * SC) + 'px ' + FONT;
    ctx.textBaseline = 'middle';

    var top = levelPx(), bot = BY + BE_H * SC;
    var rx = BX + BE_W * SC + 14 * SC;          // 标注统一排在杯外右侧
    if (nP > 2e-5) {
      var px = BX + BE_W * SC * 0.42, py = bot - (bot - top) * 0.07;
      ctx.strokeStyle = 'rgba(28,26,23,0.5)'; ctx.lineWidth = 1 * SC;
      leader(px, py - 2 * SC, rx + 8 * SC, py - 26 * SC, '蓝色沉淀', 'right');
    }
    if (nOHTot > 1e-9) {
      var over = ratio < 1.96;
      var deep = nCu4 > 2e-4;
      var sx = BX + BE_W * SC * 0.58, sy = top + (bot - top) * 0.28;
      var t = over ? '上清液：Cu²⁺'
            : deep ? '上清液：[Cu(OH)₄]²⁻'
            : '上清液：无色';
      ctx.strokeStyle = 'rgba(28,26,23,0.42)';
      leader(sx, sy, rx + 8 * SC, sy - 20 * SC, t, 'right');
    }
    ctx.restore();
  }

  var FONT = '"Songti SC","SimSun","Noto Serif SC",Georgia,serif';
  var FONT2 = '"Cascadia Mono",Consolas,monospace';

  function draw() {
    ctx.setTransform(dpr * SC, 0, 0, dpr * SC, 0, 0);
    ctx.clearRect(0, 0, cv.width / dpr, cv.height / dpr);
    drawBench();
    drawShadow();
    drawBurette();
    ctx.save();
    pathInterior(0); ctx.clip();
    pathLiquid(); ctx.clip();
    buildImage();
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(off, BX, BY, BE_W * SC, BE_H * SC);
    // 壁面暗角
    var w = BE_W * SC;
    var gl = ctx.createLinearGradient(BX, 0, BX + 13 * SC, 0);
    gl.addColorStop(0, 'rgba(20,26,32,0.20)'); gl.addColorStop(1, 'rgba(20,26,32,0)');
    ctx.fillStyle = gl; ctx.fillRect(BX, levelPx() - 4, 13 * SC, BE_H * SC);
    var gr = ctx.createLinearGradient(BX + w - 8 * SC, 0, BX + w, 0);
    gr.addColorStop(0, 'rgba(20,26,32,0)'); gr.addColorStop(1, 'rgba(20,26,32,0.14)');
    ctx.fillStyle = gr; ctx.fillRect(BX + w - 8 * SC, levelPx() - 4, 8 * SC, BE_H * SC);
    ctx.restore();
    // 液面
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
    drawLabels();
  }

  /* ============================ 面板 ============================ */

  function sig(v, n) {
    var s = v.toFixed(n);
    if (+s === 0) return '0';
    return String(+s);
  }
  function fmt(v, n) { return v.toFixed(n); }

  var elNumbers = document.getElementById('numbers');
  var elNumbers2 = document.getElementById('numbers-2');
  var elNote = document.getElementById('note');
  var elRatio = document.getElementById('ratio-v');
  var elFill = document.getElementById('gfill');
  var elMark = document.getElementById('gmarker');
  var elTag = document.getElementById('regime-tag');

  function updatePanel() {
    var L = V / 1000;
    var c2 = nCu2 / Math.max(L, 1e-9);
    var cOH = nOHf / Math.max(L, 1e-9);
    var ratio = nCuTot > 1e-9 ? nOHTot / nCuTot : null;

    var r = [];
    r.push(['加入 Cu²⁺', (nCuTot * 1000).toFixed(2) + ' mmol', 0]);
    r.push(['加入 OH⁻', (nOHTot * 1000).toFixed(2) + ' mmol', 0]);
    r.push(['混合后体积', V.toFixed(1) + ' mL', 0]);
    r.push(['生成 Cu(OH)₂', (nP * 1000).toFixed(2) + ' mmol', 0]);
    r.push(['剩余 Cu²⁺', (nCu2 * 1000).toFixed(2) + ' mmol　' + c2.toFixed(3) + ' M', nCu2 < 1e-9]);
    r.push(['剩余 OH⁻', (nOHf * 1000).toFixed(2) + ' mmol　' + cOH.toFixed(2) + ' M', nOHf < 1e-9]);
    if (nCu4 > 1e-7) r.push(['生成 [Cu(OH)₄]²⁻', (nCu4 * 1000).toFixed(2) + ' mmol　' + (nCu4 / Math.max(L, 1e-9)).toFixed(3) + ' M', 0]);
    elNumbers.innerHTML = r.map(function (q) {
      return '<span class="k">' + q[0] + '</span><span class="v' + (q[2] ? ' dim' : '') + '">' + q[1] + '</span>';
    }).join('');

    var r2 = [];
    r2.push(['剩余 OH⁻ 是否过量', (ratio !== null && ratio > 2.04) ? '是' : '否', 0]);
    r2.push(['剩余 Cu²⁺ 是否过量', (ratio !== null && ratio < 1.96 && nCu2 > 1e-9) ? '是' : '否', 0]);
    r2.push(['[Cu(OH)₄]²⁻ 是否生成', nCu4 > 1e-7 ? '是（溶液呈深蓝色）' : '否', 0]);
    elNumbers2.innerHTML = r2.map(function (q) {
      return '<span class="k">' + q[0] + '</span><span class="v' + (q[2] ? ' dim' : '') + '">' + q[1] + '</span>';
    }).join('');

    if (ratio === null) {
      elRatio.textContent = '—';
      elFill.style.width = '0%'; elMark.style.left = '0%';
    } else {
      var f = Math.max(0, Math.min(1, ratio / 3)) * 100;
      elRatio.textContent = fmt(ratio, 3);
      elFill.style.width = f + '%';
      elMark.style.left = f + '%';
      elFill.className = 'fill' + (ratio > 2 ? ' excess' : '');
    }

    var tag, note;
    if (nCuTot < 1e-9 && nOHTot < 1e-9) {
      tag = '等待加入试剂';
      note = '尚未开始。';
    } else {
      if (ratio === null) {
        tag = '滴加 CuSO₄ 溶液';
      } else if (ratio < 1.96) {
        tag = '<b>NaOH 不足</b> — CuSO₄ 过量';
      } else if (ratio <= 2.04) {
        tag = '<b>恰好完全反应</b> — 无过量离子';
      } else {
        tag = '<b>NaOH 过量</b>';
      }
      if (phase === 'pour') {
        tag += '　<span style="color:#7d766a">滴加中</span>';
      } else if (phase === 'wait') {
        tag += '　<span style="color:#7d766a">静置</span>';
      }
      note = makeNote(ratio, nCu2, nOHf, c2, cOH, nCu4);
    }
    elTag.innerHTML = tag;
    elNote.innerHTML = note;
  }

  function makeNote(ratio, nCu2, nOHf, c2, cOH, nC4) {
    if (ratio === null) {
      return '滴加 <b>CuSO₄ 溶液</b>。' + (params.order === 'cu_first'
        ? '蓝色 <span class="blue">Cu²⁺</span> 溶液充满烧杯底部。'
        : '此时烧杯内是 <b>NaOH 溶液</b>，无色。');
    }
    var col = c2 > 0.30 ? '呈蓝色' : c2 > 0.04 ? '呈浅蓝色' : '几乎无色';
    var s = '';
    if (ratio < 1.96) {
      s = 'n(NaOH) : n(CuSO₄) = <b>' + fmt(ratio, 2) + '</b>，NaOH 不足，CuSO₄ 过量。'
        + 'OH⁻ 全部用于沉淀 Cu²⁺，但 Cu²⁺ 仍有剩余（' + fmt(nCu2 * 1000, 2) + ' mmol，' + fmt(c2, 3) + ' mol·L⁻¹），'
        + '所以悬浊液' + col + '，沉淀分散在过量 CuSO₄ 溶液中。';
    } else if (ratio <= 2.04) {
      s = 'n(NaOH) : n(CuSO₄) = <b>2.00</b>，恰好按化学计量比反应。'
        + 'Cu²⁺ 与 OH⁻ 几乎全部转为 Cu(OH)₂ 沉淀，上清液是<b>无色</b>的 Na₂SO₄ 溶液。';
    } else {
      s = 'n(NaOH) : n(CuSO₄) = <b>' + fmt(ratio, 2) + '</b>，NaOH 过量。'
        + 'Cu²⁺ 已全部沉淀，多余的 OH⁻ 留在溶液中，混合后 c(OH⁻) = <b>' + fmt(cOH, 2) + ' mol·L⁻¹</b>。';
      if (nC4 > 1e-5) {
        s += ' 该浓度已足以使部分沉淀继续溶解：Cu(OH)₂ + 2OH⁻ → [Cu(OH)₄]²⁻，'
          + '含 <span class="blue">[Cu(OH)₄]²⁻</span> 的溶液呈<b>深蓝色</b>。';
      } else {
        s += ' 此浓度下 Cu(OH)₂ 不再溶解，沉淀保持<b>蓝色</b>，上清液无色。';
      }
    }
    return s;
  }

  /* ============================ 控件 ============================ */

  var sCu = document.getElementById('c-cu'), vCu = document.getElementById('v-cu');
  var sNa = document.getElementById('c-na'), vNa = document.getElementById('v-na');

  function sliderC(slider) {
    var t = +slider.value / 1000;
    return 0.05 * Math.pow(240, t);
  }

  function setSliderC(slider, c) {
    slider.value = Math.round(1000 * Math.log(c / 0.05) / Math.log(240));
  }

  function sync() {
    document.getElementById('c-cu-v').textContent = sig(params.cCu, 2);
    document.getElementById('v-cu-v').textContent = fmt(params.vCu, 1);
    document.getElementById('c-na-v').textContent = sig(params.cNa, 2);
    document.getElementById('v-na-v').textContent = fmt(params.vNa, 1);
  }

  function readAll() {
    params.cCu = sliderC(sCu);
    params.vCu = +vCu.value;
    params.cNa = sliderC(sNa);
    params.vNa = +vNa.value;
    sync(); restart();
  }
  [sCu, vCu, sNa, vNa].forEach(function (s) { s.addEventListener('input', readAll); });

  Array.prototype.forEach.call(document.querySelectorAll('[data-order]'), function (b) {
    b.addEventListener('click', function () {
      params.order = b.getAttribute('data-order');
      Array.prototype.forEach.call(document.querySelectorAll('[data-order]'), function (o) {
        o.className = (o === b) ? 'on' : '';
      });
      restart();
    });
  });

  Array.prototype.forEach.call(document.querySelectorAll('[data-preset]'), function (b) {
    b.addEventListener('click', function () {
      var p = PRESETS[b.getAttribute('data-preset')];
      params.cCu = p.cCu; params.vCu = p.vCu; params.cNa = p.cNa; params.vNa = p.vNa;
      setSliderC(sCu, params.cCu); setSliderC(sNa, params.cNa);
      vCu.value = params.vCu; vNa.value = params.vNa;
      sync(); restart();
    });
  });

  document.getElementById('btn-replay').addEventListener('click', restart);
  document.getElementById('btn-stir').addEventListener('click', function () { swirl = 1; });

  window.addEventListener('resize', resize);
  if (window.ResizeObserver) {
    new ResizeObserver(resize).observe(document.querySelector('main'));
  }

  /* ============================ 主循环 ============================ */

  var last = 0;
  function frame(t) {
    var dt = last ? Math.min(0.05, (t - last) / 1000) : 0.016;
    last = t;
    if (phase === 'pour' && !pour) {
      phase = 'wait'; phaseT = 2.0; doneBurette = true;
    } else if (phase === 'wait') {
      phaseT -= dt;
      if (phaseT <= 0) {
        if (queue.length) { phaseR = queue.shift(); startPour(phaseR); phase = 'pour'; doneBurette = false; }
        else { phase = 'done'; }
      }
    }
    if (phase === 'done') labelsOn = Math.min(1, labelsOn + dt * 1.2);

    step(dt);
    updatePanel();
    draw();
    requestAnimationFrame(frame);
  }

  /* 初始状态 = 恰好反应 */
  (function init() {
    resize();
    params.cCu = PRESETS.eq.cCu; params.vCu = PRESETS.eq.vCu;
    params.cNa = PRESETS.eq.cNa; params.vNa = PRESETS.eq.vNa;
    setSliderC(sCu, params.cCu); setSliderC(sNa, params.cNa);
    vCu.value = params.vCu; vNa.value = params.vNa;
    sync(); restart();
    requestAnimationFrame(frame);
  })();
})();
