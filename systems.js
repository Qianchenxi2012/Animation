/* 反应系统定义。每个系统都是纯数据：引擎照着执行。
 *
 * 物种 species
 *   key      内部标识
 *   name     显示名
 *   phase    'aq' 溶质（随液流输运，参与扩散）| 'solid' 固相（只沉降）
 *   absorb   溶质：[r,g,b] 每 mol·L-1 每单位光程的吸收，近似比尔-朗伯
 *   scat     固相：稀薄时散射色 [r,g,b]（0~1）
 *   scatD    固相：堆积浓密时的颜色
 *   holds    该物种结合的 token 数。每 mol 此物种占住这些 token。
 *            token 守恒式：Σ(溶质token) + Σ(固相 物种量 × holds) = 加入总量
 *            分解反应里没被任何人 hold 的 token 就是逸出了（水、气体）。
 *
 * 反应 reactions
 *   r        消耗的物种与系数
 *   p        生成的物种与系数
 *   k        速率常数，数字或函数 rate(temp, 自由基浓度)
 *   gate     门控，决定这一步在什么条件下发生，见引擎里的 GATES
 *   t0/thresh 门控与速率的阈值参数
 */

var SYSTEMS = [

/* ------------------------------------------------------------------ */
{
  id: 'cuso4-naoh',
  aim: '把硫酸铜溶液和氢氧化钠溶液先后倒入同一只烧杯，观察蓝色氢氧化铜沉淀的生成。',
  tabNote: '每摩尔铜做成 Cu(OH)₂ 要 2 mol OH⁻，做成碱式硫酸铜只要 1 mol（代价是占用 SO₄²⁻），所以总碱不足时产物就变了。',
  aimExtra: '碱不够时优先生成碱式硫酸铜；碱富裕时它又被转化回 Cu(OH)₂，并把硫酸根放回溶液。',
  molCap: '从左到右：Cu²⁺（溶解的铜，使溶液呈浅蓝）、OH⁻（无色，画成淡灰小球）、Cu(OH)₂ 固相（蓝色絮状物）。Cu²⁺ 与 OH⁻ 一相遇就结合成沉淀，所以蓝色总是先出现在两者的交界处，再慢慢散开。',
  table: [,
    ['&lt; 2', 'Cu₂(OH)₂SO₄', '淡蓝绿色沉淀，Cu²⁺ 有剩余'],
    ['= 2', 'Cu(OH)₂', '蓝色沉淀，上清液无色'],
    ['&gt; 2', 'Cu(OH)₂', '上清液仍无色'],
    ['更浓', '[Cu(OH)₄]²⁻', '溶液转深蓝，固体减少'],
  ],

  name: 'CuSO₄ 与 NaOH',
  sub: 'CuSO4 / NaOH',
  blurb: '最经典的氢氧化铜沉淀实验。碱不足时产物是碱式硫酸铜，浓碱下沉淀又会溶解。',
  metal: 'cu2',
  base: 'oh',
  heated: true,

  reagents: [
    { key: 'cu', name: '硫酸铜溶液', formula: 'CuSO₄ · aq', cMin: 0.05, cMax: 3,
      injects: { cu2: 1, so4: 1 } },
    { key: 'na', name: '氢氧化钠溶液', formula: 'NaOH · aq', cMin: 0.05, cMax: 15,
      injects: { oh: 1 } }
  ],
  orders: [
    { key: 'cu_na', label: '先 CuSO₄ 后 NaOH', seq: ['cu', 'na'] },
    { key: 'na_cu', label: '先 NaOH 后 CuSO₄', seq: ['na', 'cu'] }
  ],
  presets: [
    { key: 'eq', label: '恰好反应 2∶1', set: { cu: [1, 40], na: [2, 40] } },
    { key: 'lack', label: 'NaOH 不足', set: { cu: [0.5, 50], na: [0.5, 25] } },
    { key: 'excess', label: 'NaOH 过量', set: { cu: [1, 40], na: [3, 60] } },
    { key: 'conc', label: '浓碱溶解', set: { cu: [1.6, 30], na: [8, 40] } },
    { key: 'dilute', label: '极稀溶液', set: { cu: [0.05, 40], na: [0.06, 40] } }
  ],

  species: [
    { key: 'cu2', name: 'Cu²⁺', phase: 'aq', absorb: [1.45, 0.87, 0.23], holds: { cu: 1 } },
    { key: 'so4', name: 'SO₄²⁻', phase: 'aq', absorb: [0, 0, 0], holds: { so4: 1 } },
    { key: 'oh', name: 'OH⁻', phase: 'aq', absorb: [0, 0, 0], holds: { oh: 1 } },
    { key: 'ppt', name: 'Cu(OH)₂', phase: 'solid', density: 2.6,
      scat: [0.72, 0.86, 1.00], scatD: [0.42, 0.62, 0.93], holds: { cu: 1, oh: 2 },
      flag: '蓝色沉淀' },
    { key: 'bs', name: 'Cu₂(OH)₂SO₄', phase: 'solid', density: 1.3,
      scat: [0.56, 0.72, 0.70], scatD: [0.40, 0.55, 0.56], holds: { cu: 2, oh: 2, so4: 1 },
      flag: '碱式硫酸铜' },
    { key: 'cuo', name: 'CuO', phase: 'solid', density: 2.6,
      scat: [0.15, 0.15, 0.16], scatD: [0.09, 0.09, 0.10], holds: { cu: 1 },
      flag: 'CuO（黑色）' },
    { key: 'cu4', name: '[Cu(OH)₄]²⁻', phase: 'aq', absorb: [8.50, 4.00, 1.10],
      holds: { cu: 1, oh: 4 }, flag: '[Cu(OH)₄]²⁻（深蓝）' }
  ],

  reactions: [
    // 碱式硫酸铜这一组只在"总碱不足 2:1"时才有意义，
    // 由引擎按计划投料量统一开关（grp:'basic'）
    { r: { cu2: 2, so4: 1, oh: 2 }, p: { bs: 1 }, k: 5, gate: 'metal_rich', grp: 'basic' },
    { r: { bs: 1, oh: 2 }, p: { ppt: 2, so4: 1 }, k: 1.5, gate: 'metal_poor', grp: 'basic' },
    { r: { ppt: 1, so4: 0.5 }, p: { bs: 0.5, oh: 1 }, k: 8, gate: 'metal_rich', when: 'base_scarce', scarce: 1.6, grp: 'basic' },
    { r: { cu2: 1, oh: 2 }, p: { ppt: 1 }, k: 20, gate: 'metal_poor' },
    { r: { ppt: 1, oh: 2 }, p: { cu4: 1 }, k: 0.25, gate: 'conc', thresh: 1.7 },
    { r: { ppt: 1 }, p: { cuo: 1 }, k: 2.2, gate: 'hot', t0: 55 },
    { r: { bs: 1 }, p: { cuo: 2, so4: 1 }, k: 2.2, gate: 'hot', t0: 55 }
  ],

  eq: ['CuSO₄ + 2NaOH → Cu(OH)₂↓ + Na₂SO₄',
       '2CuSO₄ + 2NaOH → Cu₂(OH)₂SO₄↓ + Na₂SO₄',
       'Cu(OH)₂ + 2NaOH ⇌ [Cu(OH)₄]²⁻',
       'Cu(OH)₂ --Δ→ CuO↓ + H₂O'],

  ratio: { num: 'oh', den: 'cu2', stoich: 2, numName: 'NaOH', denName: 'CuSO₄', label: '投料比 n(NaOH) : n(CuSO₄)' },

  note: function (t, r) {
    if (r === null) return '滴加试剂中。';
    var s = 'n(NaOH) : n(CuSO₄) = <b>' + r.toFixed(2) + '</b>。';
    if (r < 1.96) {
      s += 'NaOH 不足。碱不足以把 Cu²⁺ 全部沉淀成 Cu(OH)₂，于是生成<b>碱式硫酸铜</b> Cu₂(OH)₂SO₄'
         + '（每摩尔 Cu 只用 1 个 OH⁻，同时占用 SO₄²⁻）。沉淀呈<b>淡蓝绿色</b>。';
    } else if (r <= 2.04) {
      s += '恰好按化学计量比反应，Cu²⁺ 与 OH⁻ 几乎全部转为 Cu(OH)₂，上清液是<b>无色</b>的 Na₂SO₄ 溶液。';
    } else {
      s += 'NaOH 过量，混合后 c(OH⁻) = <b>' + t.conc('oh').toFixed(2) + ' mol·L⁻¹</b>。';
      if (t.mm('cu4') > 0.01) {
        s += '浓度已足以使部分沉淀继续溶解生成 [Cu(OH)₄]²⁻，溶液呈<b>深蓝色</b>。';
      } else {
        s += '此浓度下 Cu(OH)₂ 不再溶解，沉淀保持<b>蓝色</b>。';
      }
    }
    if (t.mm('bs') > 0.01) s += '　碱式硫酸铜 ' + t.mm('bs').toFixed(1) + ' mmol。';
    if (t.mm('cuo') > 0.01) {
      s += '<br><b>加热后：</b>蓝色沉淀受热分解为 CuO：Cu(OH)₂ <span class="heat">Δ</span> → CuO↓ + H₂O，'
         + '已生成 CuO ' + t.mm('cuo').toFixed(1) + ' mmol。此过程不可逆，停止加热也不会变回蓝色。';
    }
    return s;
  }
},

/* ------------------------------------------------------------------ */
{
  id: 'mgcl2-naoh',
  molCap: 'Mg²⁺ 与 OH⁻ 本身都无色，画成淡灰小球；生成的 Mg(OH)₂ 是白色固相。溶液始终没有颜色——白色来自固相对光的散射，而不是溶质的颜色。',

  aim: '把氯化镁溶液和氢氧化钠溶液先后倒入同一只烧杯，观察白色氢氧化镁沉淀。',
  tabNote: 'Mg(OH)₂ 属于难溶碱，和 Cu(OH)₂ 不同，它不溶于过量碱。',
  aimExtra: '所以无论把碱加到多过量，沉淀量都不会再增加——多出来的碱只是留在溶液里。',
  table: [,
    ['&lt; 2', 'Mg(OH)₂', '部分沉淀，Mg²⁺ 有剩余'],
    ['= 2', 'Mg(OH)₂', '白色沉淀，上清液无色'],
    ['&gt; 2', 'Mg(OH)₂', '沉淀量不再增加'],
  ],

  name: 'MgCl₂ 与 NaOH',
  sub: 'MgCl2 / NaOH',
  blurb: '白色絮状沉淀，不溶于过量碱；灼烧后变成白色的氧化镁。',
  metal: 'mg2',
  base: 'oh',
  heated: true,

  reagents: [
    { key: 'mg', name: '氯化镁溶液', formula: 'MgCl₂ · aq', cMin: 0.05, cMax: 3,
      injects: { mg2: 1, cl: 1 } },
    { key: 'na', name: '氢氧化钠溶液', formula: 'NaOH · aq', cMin: 0.05, cMax: 15,
      injects: { oh: 1 } }
  ],
  orders: [
    { key: 'mg_na', label: '先 MgCl₂ 后 NaOH', seq: ['mg', 'na'] },
    { key: 'na_mg', label: '先 NaOH 后 MgCl₂', seq: ['na', 'mg'] }
  ],
  presets: [
    { key: 'eq', label: '恰好反应 2∶1', set: { mg: [0.5, 40], na: [1, 40] } },
    { key: 'lack', label: 'NaOH 不足', set: { mg: [1, 50], na: [0.5, 25] } },
    { key: 'excess', label: 'NaOH 过量', set: { mg: [0.5, 40], na: [4, 60] } },
    { key: 'conc', label: '很浓的碱', set: { mg: [1, 30], na: [12, 50] } },
    { key: 'dilute', label: '极稀溶液', set: { mg: [0.05, 40], na: [0.1, 40] } }
  ],

  species: [
    { key: 'mg2', name: 'Mg²⁺', phase: 'aq', absorb: [0, 0, 0], holds: { mg: 1 } },
    { key: 'cl', name: 'Cl⁻', phase: 'aq', absorb: [0, 0, 0], holds: { cl: 1 } },
    { key: 'oh', name: 'OH⁻', phase: 'aq', absorb: [0, 0, 0], holds: { oh: 1 } },
    { key: 'ppt', name: 'Mg(OH)₂', phase: 'solid', density: 2.6,
      scat: [0.94, 0.95, 0.96], scatD: [0.87, 0.88, 0.90], holds: { mg: 1, oh: 2 },
      flag: '白色沉淀' },
    { key: 'mgo', name: 'MgO', phase: 'solid', density: 2.6,
      scat: [0.91, 0.91, 0.92], scatD: [0.84, 0.84, 0.86], holds: { mg: 1 },
      flag: 'MgO（白色）' }
  ],

  reactions: [
    { r: { mg2: 1, oh: 2 }, p: { ppt: 1 }, k: 20 },
    { r: { ppt: 1 }, p: { mgo: 1 }, k: 2.2, gate: 'hot', t0: 55 }
  ],

  eq: ['MgCl₂ + 2NaOH → Mg(OH)₂↓ + 2NaCl',
       'Mg(OH)₂ --Δ→ MgO + H₂O'],

  ratio: { num: 'oh', den: 'mg2', stoich: 2, numName: 'NaOH', denName: 'MgCl₂', label: '投料比 n(NaOH) : n(MgCl₂)' },

  note: function (t, r) {
    if (r === null) return '滴加试剂中。';
    var s = 'n(NaOH) : n(MgCl₂) = <b>' + r.toFixed(2) + '</b>。';
    if (r < 1.96) {
      s += 'NaOH 不足，只有部分 Mg²⁺ 沉淀，剩余 Mg²⁺ '
         + t.mm('mg2').toFixed(1) + ' mmol 留在溶液里。溶液无色，沉淀是<b>白色</b>絮状物。';
    } else if (r <= 2.04) {
      s += '恰好完全反应，生成白色 Mg(OH)₂ 沉淀，上清液是<b>无色</b>的 NaCl 溶液。';
    } else {
      s += 'NaOH 过量，混合后 c(OH⁻) = <b>' + t.conc('oh').toFixed(2) + ' mol·L⁻¹</b>。'
         + 'Mg(OH)₂ <b>不溶于过量碱</b>，所以多出来的碱只是留在溶液中，沉淀量不再增加。';
    }
    if (t.mm('mgo') > 0.01) {
      s += '<br><b>加热后：</b>Mg(OH)₂ 受热分解，生成白色的 MgO：Mg(OH)₂ <span class="heat">Δ</span> → MgO + H₂O，'
         + '已生成 MgO ' + t.mm('mgo').toFixed(1) + ' mmol。Mg(OH)₂ 与 MgO 都是白色，'
         + '所以<b>肉眼看不出颜色变化</b>——判断分解是否发生要看读数。';
    }
    return s;
  }
},

/* ------------------------------------------------------------------ */
{
  id: 'fecl3-naoh',
  molCap: 'Fe³⁺ 本身就使溶液呈黄色；与 OH⁻ 结合后变成红褐色固相。滴加过程中黄色逐渐褪去，因为 Fe³⁺ 被沉淀带走了。',

  aim: '把氯化铁溶液和氢氧化钠溶液先后倒入同一只烧杯，观察红褐色氢氧化铁沉淀。',
  tabNote: 'Fe(OH)₃ 的生成比例是 1 : 3，和两价的 Cu(OH)₂、Mg(OH)₂ 不同。',
  aimExtra: 'Fe(OH)₃ 不溶于过量碱；受热脱水生成红褐色的 Fe₂O₃，两者颜色接近，肉眼看不出分解前后差别。',
  table: [,
    ['&lt; 3', 'Fe(OH)₃', '部分沉淀，Fe³⁺ 有剩余'],
    ['= 3', 'Fe(OH)₃', '红褐色沉淀，上清液几乎无色'],
    ['&gt; 3', 'Fe(OH)₃', '沉淀量不再增加'],
  ],

  name: 'FeCl₃ 与 NaOH',
  sub: 'FeCl3 / NaOH',
  blurb: '生成红褐色的氢氧化铁，是检验三价铁的特征反应。',
  metal: 'fe3',
  base: 'oh',
  heated: true,

  reagents: [
    { key: 'fe', name: '氯化铁溶液', formula: 'FeCl₃ · aq', cMin: 0.02, cMax: 2,
      injects: { fe3: 1, cl: 1 } },
    { key: 'na', name: '氢氧化钠溶液', formula: 'NaOH · aq', cMin: 0.05, cMax: 15,
      injects: { oh: 1 } }
  ],
  orders: [
    { key: 'fe_na', label: '先 FeCl₃ 后 NaOH', seq: ['fe', 'na'] },
    { key: 'na_fe', label: '先 NaOH 后 FeCl₃', seq: ['na', 'fe'] }
  ],
  presets: [
    { key: 'eq', label: '恰好反应 3∶1', set: { fe: [0.5, 40], na: [1.5, 40] } },
    { key: 'lack', label: 'NaOH 不足', set: { fe: [1, 50], na: [0.5, 25] } },
    { key: 'excess', label: 'NaOH 过量', set: { fe: [0.5, 40], na: [3, 60] } },
    { key: 'conc', label: '浓的 FeCl₃', set: { fe: [1.5, 40], na: [1.5, 20] } },
    { key: 'dilute', label: '极稀溶液', set: { fe: [0.03, 40], na: [0.09, 40] } }
  ],

  species: [
    // Fe³⁺ 水溶液本身呈黄色：主要吸收蓝光
    { key: 'fe3', name: 'Fe³⁺', phase: 'aq', absorb: [0.90, 1.50, 3.20], holds: { fe: 1 } },
    { key: 'cl', name: 'Cl⁻', phase: 'aq', absorb: [0, 0, 0], holds: { cl: 1 } },
    { key: 'oh', name: 'OH⁻', phase: 'aq', absorb: [0, 0, 0], holds: { oh: 1 } },
    { key: 'ppt', name: 'Fe(OH)₃', phase: 'solid', density: 2.6,
      scat: [0.59, 0.27, 0.16], scatD: [0.45, 0.19, 0.10], holds: { fe: 1, oh: 3 },
      flag: '红褐色沉淀' },
    { key: 'fe2o3', name: 'Fe₂O₃', phase: 'solid', density: 1.3,
      scat: [0.50, 0.22, 0.13], scatD: [0.38, 0.15, 0.08], holds: { fe: 2 },
      flag: 'Fe₂O₃（红褐）' }
  ],

  reactions: [
    { r: { fe3: 1, oh: 3 }, p: { ppt: 1 }, k: 20 },
    { r: { ppt: 2 }, p: { fe2o3: 1 }, k: 2.2, gate: 'hot', t0: 55 }
  ],

  eq: ['FeCl₃ + 3NaOH → Fe(OH)₃↓ + 3NaCl',
       '2Fe(OH)₃ --Δ→ Fe₂O₃ + 3H₂O'],

  ratio: { num: 'oh', den: 'fe3', stoich: 3, numName: 'NaOH', denName: 'FeCl₃', label: '投料比 n(NaOH) : n(FeCl₃)' },

  note: function (t, r) {
    if (r === null) return '滴加试剂中。';
    var s = 'n(NaOH) : n(FeCl₃) = <b>' + r.toFixed(2) + '</b>。';
    if (r < 2.85) {
      s += 'NaOH 不足，剩余 Fe³⁺ ' + t.mm('fe3').toFixed(1) + ' mmol 留在溶液中，'
         + '溶液仍呈<b>黄色</b>。沉淀是<b>红褐色</b>的 Fe(OH)₃。';
    } else if (r <= 3.15) {
      s += '恰好完全反应，生成红褐色 Fe(OH)₃ 沉淀。FeCl₃ 溶液的黄色被沉淀带走，上清液几乎无色。';
    } else {
      s += 'NaOH 过量，c(OH⁻) = <b>' + t.conc('oh').toFixed(2) + ' mol·L⁻¹</b>。'
         + 'Fe(OH)₃ <b>不溶于过量碱</b>（与两性氢氧化物不同），沉淀量不再增加。';
    }
    if (t.mm('fe2o3') > 0.01) {
      s += '<br><b>加热后：</b>Fe(OH)₃ 脱水生成 Fe₂O₃：2Fe(OH)₃ <span class="heat">Δ</span> → Fe₂O₃ + 3H₂O，'
         + '已生成 ' + t.mm('fe2o3').toFixed(1) + ' mmol。两者都是红褐色，颜色上看不出差别。';
    }
    return s;
  }
},

/* ------------------------------------------------------------------ */
{
  id: 'bacl2-na2so4',
  molCap: 'Ba²⁺ 与 SO₄²⁻ 都无色，画成淡灰小球；两者结合生成白色 BaSO₄。溶解度极低，因此把某一方加到过量，沉淀量也不再增加。',

  aim: '把氯化钡溶液和硫酸钠溶液先后倒入同一只烧杯，观察白色硫酸钡沉淀。',
  tabNote: 'BaSO₄ 的溶解度比上面几种都小得多，是检验硫酸根的常用反应。',
  aimExtra: '它既不溶于过量试剂，也不受加热影响，所以调比例或加热都不会改变沉淀量。',
  table: [,
    ['&lt; 1', 'BaSO₄', '部分沉淀，Ba²⁺ 有剩余'],
    ['= 1', 'BaSO₄', '白色沉淀'],
    ['&gt; 1', 'BaSO₄', '沉淀量不再增加'],
  ],

  name: 'BaCl₂ 与 Na₂SO₄',
  sub: 'BaCl2 / Na2SO4',
  blurb: '生成极难溶的白色硫酸钡，是检验硫酸根的常用反应。',
  metal: 'ba2',
  base: null,
  heated: false,

  reagents: [
    { key: 'ba', name: '氯化钡溶液', formula: 'BaCl₂ · aq', cMin: 0.02, cMax: 2,
      injects: { ba2: 1, cl: 1 } },
    { key: 'na', name: '硫酸钠溶液', formula: 'Na₂SO₄ · aq', cMin: 0.02, cMax: 2,
      injects: { so4: 1 } }
  ],
  orders: [
    { key: 'ba_na', label: '先 BaCl₂ 后 Na₂SO₄', seq: ['ba', 'na'] },
    { key: 'na_ba', label: '先 Na₂SO₄ 后 BaCl₂', seq: ['na', 'ba'] }
  ],
  presets: [
    { key: 'eq', label: '恰好反应 1∶1', set: { ba: [1, 40], na: [1, 40] } },
    { key: 'lack', label: 'Na₂SO₄ 不足', set: { ba: [1, 50], na: [0.5, 25] } },
    { key: 'excess', label: 'Na₂SO₄ 过量', set: { ba: [0.5, 40], na: [2, 60] } },
    { key: 'conc', label: '浓溶液', set: { ba: [2, 30], na: [2, 30] } },
    { key: 'dilute', label: '极稀溶液', set: { ba: [0.05, 40], na: [0.05, 40] } }
  ],

  species: [
    { key: 'ba2', name: 'Ba²⁺', phase: 'aq', absorb: [0, 0, 0], holds: { ba: 1 } },
    { key: 'cl', name: 'Cl⁻', phase: 'aq', absorb: [0, 0, 0], holds: { cl: 1 } },
    { key: 'so4', name: 'SO₄²⁻', phase: 'aq', absorb: [0, 0, 0], holds: { so4: 1 } },
    { key: 'baso4', name: 'BaSO₄', phase: 'solid', density: 2.6,
      scat: [0.96, 0.96, 0.97], scatD: [0.89, 0.89, 0.91], holds: { ba: 1, so4: 1 },
      flag: '白色沉淀' }
  ],

  reactions: [
    { r: { ba2: 1, so4: 1 }, p: { baso4: 1 }, k: 20 }
  ],

  eq: ['BaCl₂ + Na₂SO₄ → BaSO₄↓ + 2NaCl'],

  ratio: { num: 'so4', den: 'ba2', stoich: 1, numName: 'Na₂SO₄', denName: 'BaCl₂', label: '投料比 n(Na₂SO₄) : n(BaCl₂)' },

  note: function (t, r) {
    if (r === null) return '滴加试剂中。';
    var s = 'n(Na₂SO₄) : n(BaCl₂) = <b>' + r.toFixed(2) + '</b>。';
    if (r < 0.96) {
      s += 'Na₂SO₄ 不足，剩余 Ba²⁺ ' + t.mm('ba2').toFixed(1) + ' mmol 留在溶液中。';
    } else if (r <= 1.04) {
      s += '恰好完全反应，生成白色 BaSO₄ 沉淀。';
    } else {
      s += 'Na₂SO₄ 过量，剩余 ' + t.mm('so4').toFixed(1) + ' mmol。';
    }
    return s + 'BaSO₄ <b>既不溶于过量试剂，也不受加热影响</b>，'
         + '所以无论怎么调比例或加热，沉淀量都不再增加——这正是它被用来检验硫酸根的原因。';
  }
},

/* ------------------------------------------------------------------ */
{
  id: 'cuso4-nh3',
  molCap: 'Cu²⁺ 与 NH₃ 结合先得到蓝色 Cu(OH)₂；氨过量时沉淀继续溶解，NH₃ 配位到 Cu²⁺ 周围，形成深蓝色的铜氨配离子。',

  aim: '把硫酸铜溶液和氨水先后倒入同一只烧杯。氨水是弱碱，先看到蓝色沉淀。',
  tabNote: '氨过量时 Cu(OH)₂ 会继续溶解，NH₃ 配位到 Cu²⁺ 周围，形成深蓝色的铜氨配离子。',
  aimExtra: '这和「碱过量」不同：Cu(OH)₂ 在过量氨水里是溶解的——所以沉淀会越来越少，溶液反而更深。',
  table: [,
    ['&lt; 2', 'Cu(OH)₂', '蓝色沉淀，Cu²⁺ 有剩余'],
    ['= 2', 'Cu(OH)₂', '蓝色沉淀'],
    ['&gt; 2 且较浓', '[Cu(NH₃)₄]²⁺', '沉淀减少，溶液转深蓝'],
  ],

  name: 'CuSO₄ 与氨水',
  sub: 'CuSO4 / NH3',
  blurb: '先出现蓝色沉淀，氨水过量时沉淀溶解成深蓝色的铜氨配离子。',
  metal: 'cu2',
  base: 'nh3',
  heated: false,

  reagents: [
    { key: 'cu', name: '硫酸铜溶液', formula: 'CuSO₄ · aq', cMin: 0.05, cMax: 2,
      injects: { cu2: 1, so4: 1 } },
    { key: 'am', name: '氨水', formula: 'NH₃ · H₂O', cMin: 0.05, cMax: 12,
      injects: { nh3: 1 } }
  ],
  orders: [
    { key: 'cu_am', label: '先 CuSO₄ 后氨水', seq: ['cu', 'am'] },
    { key: 'am_cu', label: '先氨水后 CuSO₄', seq: ['am', 'cu'] }
  ],
  presets: [
    { key: 'eq', label: '恰好生成沉淀', set: { cu: [1, 40], am: [2, 40] } },
    { key: 'lack', label: '氨水不足', set: { cu: [1, 50], am: [1, 25] } },
    { key: 'excess', label: '氨水过量', set: { cu: [1, 40], am: [6, 40] } },
    { key: 'conc', label: '很浓的氨水', set: { cu: [1, 30], am: [12, 40] } },
    { key: 'dilute', label: '极稀溶液', set: { cu: [0.05, 40], am: [0.1, 40] } }
  ],

  species: [
    { key: 'cu2', name: 'Cu²⁺', phase: 'aq', absorb: [1.45, 0.87, 0.23], holds: { cu: 1 } },
    { key: 'so4', name: 'SO₄²⁻', phase: 'aq', absorb: [0, 0, 0], holds: { so4: 1 } },
    { key: 'nh3', name: 'NH₃', phase: 'aq', absorb: [0, 0, 0], holds: { nh3: 1 } },
    // 氨水是弱碱，反应的另一头是铵离子
    { key: 'nh4', name: 'NH₄⁺', phase: 'aq', absorb: [0, 0, 0], holds: { nh3: 1 } },
    { key: 'ppt', name: 'Cu(OH)₂', phase: 'solid', density: 2.6,
      scat: [0.72, 0.86, 1.00], scatD: [0.42, 0.62, 0.93], holds: { cu: 1 },
      flag: '蓝色沉淀' },
    // 铜氨配离子的蓝比铜(II)本身鲜亮得多
    { key: 'am', name: '[Cu(NH₃)₄]²⁺', phase: 'aq', absorb: [6.00, 2.60, 0.40],
      holds: { cu: 1, nh3: 4 }, flag: '[Cu(NH₃)₄]²⁺（深蓝）' }
  ],

  reactions: [
    { r: { cu2: 1, nh3: 2 }, p: { ppt: 1, nh4: 2 }, k: 20 },
    { r: { ppt: 1, nh3: 4 }, p: { am: 1 }, k: 0.6, gate: 'conc', thresh: 1.2 }
  ],

  eq: ['CuSO₄ + 2NH₃·H₂O → Cu(OH)₂↓ + (NH₄)₂SO₄',
       'Cu(OH)₂ + 4NH₃ ⇌ [Cu(NH₃)₄]²⁺ + 2OH⁻'],

  ratio: { num: 'nh3', den: 'cu2', stoich: 2, numName: 'NH₃', denName: 'CuSO₄', label: '投料比 n(NH₃) : n(CuSO₄)' },

  note: function (t, r) {
    if (r === null) return '滴加试剂中。';
    var s = 'n(NH₃) : n(CuSO₄) = <b>' + r.toFixed(2) + '</b>。';
    if (r < 1.96) {
      s += '氨水不足，剩余 Cu²⁺ ' + t.mm('cu2').toFixed(1) + ' mmol，溶液呈蓝色。';
    } else if (r <= 2.04) {
      s += '氨水刚好够生成蓝色 Cu(OH)₂ 沉淀。';
    } else {
      s += '氨水过量。Cu(OH)₂ 会<b>溶于过量氨水</b>，生成深蓝色的铜氨配离子：'
         + 'Cu(OH)₂ + 4NH₃ ⇌ [Cu(NH₃)₄]²⁺ + 2OH⁻。'
         + '蓝色沉淀逐渐减少，溶液转为<b>纯深蓝色</b>。';
    }
    if (t.mm('am') > 0.01) s += '　已生成铜氨配离子 ' + t.mm('am').toFixed(1) + ' mmol。';
    return s;
  }
}

];