# CuSO₄ 与 NaOH — 浓度对反应的影响

氢氧化钠与硫酸铜在不同浓度下反应的交互模拟。调节两种溶液的浓度和体积，观察蓝色氢氧化铜沉淀的生成、沉降，以及浓碱条件下生成深蓝色 `[Cu(OH)₄]²⁻` 的过程。

纯前端，无构建步骤，无运行时依赖。

## 在线演示

<https://user.github.io/animation/> （部署后把链接换成你自己的）

## 本地运行

直接双击 `index.html` 即可，无需服务器。

需要本地服务器时（比如想避免某些浏览器的 `file://` 限制）：

```sh
python -m http.server 8000
# 或
npx serve .
```

然后打开 <http://localhost:8000>。

## 演示的化学

```
CuSO₄ + 2NaOH → Cu(OH)₂↓ + Na₂SO₄
Cu(OH)₂ + 2NaOH ⇌ [Cu(OH)₄]²⁻
```

界面内置五种典型配比，点击即可依次播放：

| 配比 | 现象 |
|---|---|
| 恰好反应 2∶1 | 蓝色沉淀，上清液无色 |
| NaOH 不足 | 悬浊液呈浅蓝色，Cu²⁺ 有剩余 |
| NaOH 过量 | 沉淀不变，上清液仍无色 |
| 浓碱溶解 | 溶液转深蓝，部分沉淀溶解成 `[Cu(OH)₄]²⁻` |
| 极稀溶液 | 沉淀少到几乎看不见 |

浓度、体积、加入顺序都可手动调节；「摇动烧杯」会把沉淀重新扬起。

## 模拟方式

不是播放动画，而是一个逐格计算的场模型。烧杯内部分成 56 × 92 的网格，每格记录
`Cu²⁺`、`OH⁻`、`Cu(OH)₂` 固体和 `[Cu(OH)₄]²⁻` 的浓度，每帧依次做：

- **平流** — 半拉格朗日法输运溶质，取样不越过液面，再做有界总量修正
- **扩散** — 显式格式，分两次子步保证稳定
- **反应** — 局部接触即反应；`Cu²⁺ + 2OH⁻ → Cu(OH)₂`，浓碱下 `Cu(OH)₂ + 2OH⁻ → [Cu(OH)₄]²⁻`
- **湍流松弛** — 让溶液最终趋于均匀，同时短时间内保留羽状层次
- **沉降** — 沉淀用守恒通量下沉，装满堆积密度后停住，于是自然堆出密实层

滚动的环流由流函数导出，天然无散。物质的量由场直接统计，因此读数与画面始终一致。

溶液颜色用比尔-朗伯吸收加散射近似，吸收系数按观感调校，并非从真实吸收光谱计算。

## 浏览器支持

需要 `canvas`、`ResizeObserver`、CSS Grid。Chrome / Edge / Firefox / Safari 近两年的版本均可。

## 部署到 GitHub Pages

仓库是纯静态站点，推上去就能用：

```sh
git init
git add .
git commit -m "CuSO4 与 NaOH 反应模拟"
git branch -M main
git remote add origin https://github.com/<用户名>/<仓库名>.git
git push -u origin main
```

然后在仓库 **Settings → Pages** 里把 Source 选成 **Deploy from a branch**，
分支选 `main`、目录选 `/ (root)`，保存后等一两分钟即可访问
`https://<用户名>.github.io/<仓库名>/`。

`.nojekyll` 是为了让 GitHub Pages 不走 Jekyll 处理。仓库根目录下有 `index.html`，
所以选 `/ (root)` 就能直接跑，不需要任何构建或配置。

<details>
<summary>其他两种部署方式</summary>

- **`docs/` 目录**：把上面那些文件放进 `docs/`，Pages 的 Source 里选
  `/docs`。适合想让仓库根目录放 README、截图、代码说明的情况。
- **用户主页站点**：仓库名必须是 `<用户名>.github.io`，此时访问
  `https://<用户名>.github.io/`（不带仓库名那一层）。

因为页面里所有资源都是相对路径（`style.css`、`main.js`、`favicon.svg`），
放在子路径、换域名、或者中间再加一层目录都不用改代码。

</details>

## 文件

```
index.html    页面结构
style.css     样式
main.js       模拟与绘制
favicon.svg   图标
.nojekyll     关闭 Jekyll 处理
LICENSE       MIT
```

`.gitignore` 会挡掉 `opencode.json` 和截图等本地产物，不用手动清理。

## 许可

MIT，见 [LICENSE](LICENSE)。