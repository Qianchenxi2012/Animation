/* 界面层：按 systems.js 的数据生成控件、读数和结论，模拟交给引擎。 */

(function () {
  'use strict';

  var eng = window.ReactionEngine.create(document.getElementById('scene'));
  var sys = null;

  var el = {
    sysList: document.getElementById('sys-list'),
    blurb: document.getElementById('sys-blurb'),
    presets: document.getElementById('presets'),
    reagents: document.getElementById('reagents'),
    orders: document.getElementById('orders'),
    eq: document.getElementById('eq'),
    tag: document.getElementById('regime-tag'),
    numbers: document.getElementById('numbers'),
    numbers2: document.getElementById('numbers-2'),
    ratioV: document.getElementById('ratio-v'),
    gfill: document.getElementById('gfill'),
    gmark: document.getElementById('gmark'),
    note: document.getElementById('note'),
    btnHeat: document.getElementById('btn-heat'),
    btnReplay: document.getElementById('btn-replay'),
    btnStir: document.getElementById('btn-stir'),
    btnInfo: document.getElementById('btn-info'),
    info: document.getElementById('info'),
    molview: document.getElementById('molview'),
    molCap: document.getElementById('mol-cap')
  };

  function sig(v, n) { var s = v.toFixed(n); return String(+s); }

  /* ---------- 系统选择 ---------- */

  function buildSysList() {
    el.sysList.innerHTML = '';
    window.SYSTEMS.forEach(function (s) {
      var b = document.createElement('button');
      b.textContent = s.name;
      b.setAttribute('data-sys', s.id);
      b.addEventListener('click', function () { pick(s.id); });
      el.sysList.appendChild(b);
    });
  }

  function pick(id) {
    window.SYSTEMS.forEach(function (s) { if (s.id === id) sys = s; });
    Array.prototype.forEach.call(el.sysList.children, function (b) {
      b.className = (b.getAttribute('data-sys') === id) ? 'on' : '';
    });
    eng.load(sys, sys.presets[0]);
    el.blurb.textContent = sys.blurb;
    document.getElementById('sys-title').textContent = sys.name;
    document.getElementById('sys-sub').textContent = sys.sub;
    document.title = sys.name + ' — 沉淀反应模拟';
    document.getElementById('ratio-label').textContent = sys.ratio.label;
    buildPresets();
    buildReagents();
    buildOrders();
    buildEq();
    el.molCap.textContent = '左起依次是本实验涉及的粒子：' + (sys.molCap || '');
    buildInfo();
    document.getElementById('model-note').textContent =
      '这不是播放动画，而是一个逐格计算的场模型：烧杯分成 56 × 92 的网格，'
      + '每格记录各物种浓度，每帧依次做平流、扩散、局部反应与沉降。'
      + '物质的量由场直接统计，所以读数与画面始终一致，且严格守恒。'
      + '（本系统有 ' + sys.species.length + ' 个物种、' + sys.reactions.length + ' 个反应。）';
    el.btnHeat.parentNode.style.display = sys.heated ? '' : 'none';
    eng.setHeating(false);
    el.btnHeat.textContent = '加热';
    el.btnHeat.className = '';
    update();
  }

  function buildPresets() {
    el.presets.innerHTML = '';
    sys.presets.forEach(function (p) {
      var b = document.createElement('button');
      b.textContent = p.label;
      b.addEventListener('click', function () {
        for (var k in p.set) {
          var rd = sys.reagentByKey[k];
          rd.conc = p.set[k][0];
          rd.vol = p.set[k][1];
        }
        syncSliders();
        eng.restart();
      });
      el.presets.appendChild(b);
    });
  }

  var sliders = {};

  function buildReagents() {
    el.reagents.innerHTML = '';
    sliders = {};
    sys.reagents.forEach(function (rd) {
      rd.conc = rd.conc || 1;
      rd.vol = rd.vol || 40;
      var wrap = document.createElement('div');
      wrap.className = 'group';
      var head = document.createElement('div');
      head.className = 'chem-head';
      head.innerHTML = '<span class="name">' + rd.name + '</span><span class="mono">' + rd.formula + '</span>';
      wrap.appendChild(head);

      var f1 = document.createElement('div');
      f1.className = 'field';
      f1.innerHTML = '<label><span>浓度</span><span><span class="val" data-v="c"></span> ' +
                     '<span class="unit">mol·L⁻¹</span></span></label>';
      var s1 = document.createElement('input');
      s1.type = 'range'; s1.min = '0'; s1.max = '1000'; s1.step = '1';
      f1.appendChild(s1);
      wrap.appendChild(f1);

      var f2 = document.createElement('div');
      f2.className = 'field';
      f2.innerHTML = '<label><span>取用体积</span><span><span class="val" data-v="v"></span> ' +
                     '<span class="unit">mL</span></span></label>';
      var s2 = document.createElement('input');
      s2.type = 'range'; s2.min = '5'; s2.max = '80'; s2.step = '0.5';
      f2.appendChild(s2);
      wrap.appendChild(f2);

      el.reagents.appendChild(wrap);

      var f = s1, v = s2;
      f.addEventListener('input', function () {
        rd.conc = fromSlider(f.value, rd);
        syncSliders(); eng.restart();
      });
      v.addEventListener('input', function () {
        rd.vol = +v.value;
        syncSliders(); eng.restart();
      });
      sliders[rd.key] = { c: f, v: v, cv: f1.querySelector('[data-v=c]'), vv: f2.querySelector('[data-v=v]') };
    });
    syncSliders();
  }

  function fromSlider(t, rd) {
    var x = +t / 1000;
    return rd.cMin * Math.pow(rd.cMax / rd.cMin, x);
  }

  function toSlider(c, rd) {
    return Math.round(1000 * Math.log(c / rd.cMin) / Math.log(rd.cMax / rd.cMin));
  }

  function syncSliders() {
    sys.reagents.forEach(function (rd) {
      var s = sliders[rd.key];
      if (!s) return;
      s.c.value = toSlider(rd.conc, rd);
      s.v.value = rd.vol;
      s.cv.textContent = sig(rd.conc, 2);
      s.vv.textContent = rd.vol.toFixed(1);
    });
  }

  function buildOrders() {
    el.orders.innerHTML = '';
    sys.orders.forEach(function (o, i) {
      var b = document.createElement('button');
      b.textContent = o.label;
      if (i === 0) b.className = 'on';
      b.setAttribute('data-ord', o.key);
      b.addEventListener('click', function () {
        sys.orderKey = o.key;
        Array.prototype.forEach.call(el.orders.children, function (x) {
          x.className = (x === b) ? 'on' : '';
        });
        eng.restart();
      });
      el.orders.appendChild(b);
    });
    sys.orderKey = sys.orders[0].key;
  }

  function buildEq() {
    el.eq.innerHTML = sys.eq.map(function (line) {
      if (line.indexOf('Δ') >= 0) {
        line = line.replace('--Δ→', '<span class="heat">Δ</span> →');
      }
      return line;
    }).join('<br>');
  }

  // 说明浮层里随系统变化的部分
  function buildInfo() {
    document.getElementById('sys-aim').innerHTML = sys.aim || '';
    document.getElementById('sys-extra').innerHTML = sys.aimExtra || '';
    document.getElementById('sys-tabnote').innerHTML = sys.tabNote || '';
    var eqs = sys.eq.map(function (line) {
      if (line.indexOf('Δ') >= 0) line = line.replace('--Δ→', '<span class="heat">Δ</span> →');
      return line;
    }).join('<br>');
    document.getElementById('sys-eq').innerHTML = eqs;
    var rows = '<tr><th>' + sys.ratio.numName + ' : ' + sys.ratio.denName +
               '</th><th>产物</th><th>现象</th></tr>';
    (sys.table || []).forEach(function (r) {
      rows += '<tr><td>' + r[0] + '</td><td>' + r[1] + '</td><td>' + r[2] + '</td></tr>';
    });
    document.getElementById('sys-tab').innerHTML = rows;
  }

  /* ---------- 读数 ---------- */

  // 物种 → 它携带的 token（投料比用）
  function tokenOf(spKey) {
    var sp = sys.byKey[spKey];
    return sp && sp.holds ? Object.keys(sp.holds)[0] : spKey;
  }

  function ledgerTally() {
    var t = {
      mm: function (k) { var sp = sys.byKey[k]; return sp ? (sp.mol || 0) * 1000 : 0; },
      conc: function (k) {
        var sp = sys.byKey[k];
        var L = eng.V / 1000;
        return sp && L > 0 ? (sp.mol || 0) / L : 0;
      }
    };
    return t;
  }

  function update() {
    var S = eng.S();
    if (!S) return;
    var t = ledgerTally();
    var st = eng.state();
    var L = eng.V / 1000;

    // 各试剂带入的量
    var rows = [];
    S.reagents.forEach(function (rd) {
      var n = rd.conc * rd.vol / 1000 * 1000;
      rows.push(['加入 ' + rd.name.replace('溶液', ''), n.toFixed(2) + ' mmol', 0]);
    });
    rows.push(['混合后体积', eng.V.toFixed(1) + ' mL', 0]);
    if (S.heated) {
      rows.push(['温度', eng.temp.toFixed(0) + ' °C' + (st.heating ? '' : ''), 0]);
    }
    S.solids.forEach(function (sp) {
      var m = t.mm(sp.key);
      if (m > 0.005) rows.push([sp.name + ' 固体', m.toFixed(2) + ' mmol', 0]);
    });
    S.solutes.forEach(function (sp) {
      var m = t.mm(sp.key);
      if (m <= 0.005) return;
      if (sp.absorb && (sp.absorb[0] + sp.absorb[1] + sp.absorb[2]) > 0) {
        rows.push(['剩余 ' + sp.name, m.toFixed(2) + ' mmol　' + t.conc(sp.key).toFixed(3) + ' M', 0]);
      } else {
        rows.push([sp.name + ' 游离量', m.toFixed(2) + ' mmol', 0]);
      }
    });
    el.numbers.innerHTML = rows.map(function (q) {
      return '<span class="k">' + q[0] + '</span><span class="v' + (q[2] ? ' dim' : '') + '">' + q[1] + '</span>';
    }).join('');

    // 投料比按试剂的计划投料量算：分母若用会被消耗掉的自由离子，
    // 沉淀完全后会翻成 null，结论就断了。
    var pl = eng.plan();
    var den = pl[tokenOf(S.ratio.den)] || 0;
    var num = pl[tokenOf(S.ratio.num)] || 0;
    var ratio = den > 1e-9 ? num / den : null;
    if (ratio === null) {
      el.ratioV.textContent = '—';
      el.gfill.style.width = '0%'; el.gmark.style.left = '0%';
    } else {
      var f = Math.max(0, Math.min(1, ratio / 3)) * 100;
      el.ratioV.textContent = ratio.toFixed(3);
      el.gfill.style.width = f + '%';
      el.gmark.style.left = f + '%';
      el.gfill.className = 'fill' + (ratio > S.ratio.stoich ? ' excess' : '');
    }

    // 判定
    var flags = [];
    S.solids.forEach(function (sp) {
      if (sp.flag) flags.push(['是否生成 ' + sp.name, t.mm(sp.key) > 0.01 ? '是' : '否', 0]);
    });
    S.solutes.forEach(function (sp) {
      if (sp.flag) flags.push(['是否生成 ' + sp.name, t.mm(sp.key) > 0.01 ? '是' : '否', 0]);
    });
    el.numbers2.innerHTML = flags.map(function (q) {
      return '<span class="k">' + q[0] + '</span><span class="v">' + q[1] + '</span>';
    }).join('');

    // 标题与结论
    var tag;
    if (st.V < 0.5) {
      tag = '等待加入试剂';
    } else if (ratio === null) {
      tag = st.pouring ? ('滴加 ' + S.reagentByKey[st.pouring].name) : '静置';
    } else if (ratio < S.ratio.stoich - 0.04) {
      tag = '<b>' + S.ratio.denName + ' 不足</b>';
    } else if (ratio <= S.ratio.stoich + 0.04) {
      tag = '<b>恰好完全反应</b>';
    } else {
      tag = '<b>' + S.ratio.numName + ' 过量</b>';
    }
    if (st.phase === 'pour') tag += '　<span style="color:#7d766a">滴加中</span>';
    else if (st.phase === 'wait') tag += '　<span style="color:#7d766a">静置</span>';
    el.tag.innerHTML = tag;
    el.note.innerHTML = sys.note(t, ratio);
  }

  /* ---------- 按钮 ---------- */

  el.btnReplay.addEventListener('click', function () { eng.restart(); });
  el.btnStir.addEventListener('click', function () { eng.setSwirl(1); });
  el.btnHeat.addEventListener('click', function () {
    var on = el.btnHeat.textContent === '加热';
    eng.setHeating(on);
    el.btnHeat.textContent = on ? '停止加热' : '加热';
    el.btnHeat.className = on ? 'on' : '';
  });
  el.btnInfo.addEventListener('click', function () {
    var open = el.info.classList.toggle('open');
    el.btnInfo.setAttribute('aria-expanded', open ? 'true' : 'false');
    el.btnInfo.textContent = open ? '收起说明' : '实验说明';
  });

  /* ---------- 主循环 ---------- */

  var last = 0;
  function frame(t) {
    var dt = last ? Math.min(0.05, (t - last) / 1000) : 0.016;
    last = t;
    if (sys) {
      eng.advance(dt);
      eng.step(dt);
      eng.tickLabels(dt);
      update();
      eng.draw();
      // 说明栏展开时才画，省开销
      if (el.info.classList.contains('open')) {
        window.MoleculeView.draw(el.molview, sys, t / 1000);
      }
    }
    requestAnimationFrame(frame);
  }

  window.addEventListener('resize', function () { eng.resize(); });
  if (window.ResizeObserver) {
    new ResizeObserver(function () { eng.resize(); }).observe(document.querySelector('main'));
  }

  buildSysList();
  eng.resize();
  pick(window.SYSTEMS[0].id);
  requestAnimationFrame(frame);

  // 便于自动化检查
  window.__app = {
    eng: eng,
    pick: pick,
    list: function () { return window.SYSTEMS.map(function (s) { return s.id; }); },
    sys: function () { return sys; },
    state: function () { return eng.state(); },
    ledger: function () { return eng.ledger(); }
  };
})();