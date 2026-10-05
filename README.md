# 沉淀反应模拟器

氢氧化钠、氨水等碱与金属盐溶液反应的交互模拟。调节浓度和体积，观察沉淀的生成、
沉降、受热分解，以及过量试剂带来的产物变化。

纯前端，无构建步骤，无运行时依赖。

## 在线演示

<https://user.github.io/animation/> （部署后把链接换成你自己的）

## 内置实验

左上角可以切换。当前五个：

| 实验 | 沉淀 | 溶液 | 特点 |
|---|---|---|---|
| CuSO₄ 与 NaOH | 蓝色 Cu(OH)₂ | 浅蓝 | 碱不足生成**碱式硫酸铜**；浓碱下溶解成深蓝 `[Cu(OH)₄]²⁻`；加热变黑 CuO |
| MgCl₂ 与 NaOH | 白色 Mg(OH)₂ | 无色 | 不溶于过量碱；灼烧得白色 MgO |
| FeCl₃ 与 NaOH | 红褐色 Fe(OH)₃ | 黄色 | 1 : 3 生成；受热脱水成红褐 Fe₂O₃ |
| BaCl₂ 与 Na₂SO₄ | 白色 BaSO₄ | 无色 | 极难溶，过量试剂和加热都不改变沉淀量 |
| CuSO₄ 与氨水 | 蓝色 Cu(OH)₂ | 深蓝 | 氨过量时沉淀**溶解**成铜氨配离子 |

每个系统都有五组典型配比：恰好反应、试剂不足、试剂过量、浓溶液、极稀溶液。

## 本地运行

直接双击 `index.html` 即可，无需服务器。

需要本地服务器时：

```sh
python -m http.server 8000
# 或
npx serve .
```

## 模拟方式

不是播放动画，而是逐格计算的场模型。烧杯内部分成 56 × 92 的网格，每格记录各物种
浓度，每帧依次做：

- **平流** — 半拉格朗日法输运溶质，取样不越过液面，再做有界总量修正
- **扩散** — 显式格式，分两次子步保证稳定
- **反应** — 通用的化学计量内核，按配方逐条算反应进度
- **湍流松弛** — 让溶液最终趋于均匀，同时短时间内保留羽状层次
- **沉降** — 固相用守恒通量下沉，装满堆积密度后停住，于是自然堆出密实层

环流由流函数导出，天然无散。

### 物质守恒怎么保证

每个物种声明自己**结合了哪些 token**（Cu、OH、SO₄、Mg、NH₃……），例如：

```js
{ key: 'ppt', name: 'Cu(OH)₂', phase: 'solid',
  holds: { cu: 1, oh: 2 } }        // 每摩尔沉淀占住 1 个铜、2 个氢氧根
```

反应只搬运物种，不直接改账。守恒式由 `holds` 自动给出：

```
加入总量 = Σ(游离 token) + Σ(物种量 × holds) + 逸出量
```

分解反应里没被任何人 hold 的 token 就是逸出了（生成水）。所以不存在"反应写错
导致凭空多出物质"的可能——一旦多算，界面的 token 台账会立刻显示误差。

### 反应由数据描述

`systems.js` 里每个实验就是一份数据，没有分支代码：

```js
{ r: { cu2: 2, so4: 1, oh: 2 }, p: { bs: 1 }, k: 5, gate: 'metal_rich', grp: 'basic' }
```

引擎按配平式算反应进度（受最紧的 reactant 限制），乘上速率常数和门控。
可用的门控：`metal_rich` / `metal_poor`（按金属与碱的局部比例）、`hot`（温度阈值）、
`conc`（浓度阈值）。`grp` 用来成组开关——比如碱式硫酸铜那一组只在**总碱不足 2:1**
时才有意义。

新增实验只需往 `systems.js` 里加一份数据，界面、微观示意、读数、结论文字都会自动生成。

## 微观示意

说明浮层里的粒子图由同一份物种数据画出：金属离子、碱、固相、配离子各自一格，
小球按固相/溶液分层缓慢浮动。无色离子画成淡灰色，有颜色的按吸收系数反推色深。

## 浏览器支持

需要 `canvas`、`ResizeObserver`、CSS Grid。Chrome / Edge / Firefox / Safari 近两年的版本均可。

## 部署到 GitHub Pages

仓库是纯静态站点，推上去就能用：

```sh
git init
git add .
git commit -m "沉淀反应模拟器"
git branch -M main
git remote add origin https://github.com/<用户名>/<仓库名>.git
git push -u origin main
```

然后在仓库 **Settings → Pages** 里把 Source 选成 **Deploy from a branch**，
分支选 `main`、目录选 `/ (root)`，保存后等一两分钟即可访问
`https://<用户名>.github.io/<仓库名>/`。

`.nojekyll` 是为了让 GitHub Pages 不走 Jekyll 处理。页面里所有资源都是相对路径，
所以放在子路径、换域名、中间再加一层目录都不用改代码。

<details>
<summary>其他两种部署方式</summary>

- **`docs/` 目录**：把那些文件放进 `docs/`，Pages 的 Source 里选 `/docs`。
- **用户主页站点**：仓库名必须是 `<用户名>.github.io`，访问 `https://<用户名>.github.io/`。

</details>

## 文件

```
index.html    页面结构
style.css     样式
systems.js    五个实验的数据定义
engine.js     模拟与绘制引擎（由数据驱动）
molecule.js   微观粒子示意
main.js       界面层
favicon.svg   图标
.nojekyll     关闭 Jekyll 处理
LICENSE       MIT
```

`.gitignore` 会挡掉 `opencode.json` 和截图等本地产物。

## 已知简化

- 溶液颜色用比尔-朗伯吸收加散射近似，吸收系数按观感调校，**并非从真实吸收光谱计算**。
- 沉淀的絮状纹理、沉降速度、沸腾翻滚都做了视觉上的简化。
- 不含活度系数与离子强度的影响，浓度很高时与真实溶液有偏差。

## 许可

MIT，见 [LICENSE](LICENSE)。
