# 稀释配液谱系与原样浓度回算仓（PREPQUANT）

用 Python 3.11 + SQLite 实现的**独立**配液账与样品定量仓。两个组件只通过
“导出产品”这一外部格式耦合，全部数值为 `fractions.Fraction` 精确有理数，
不连接任何设备，不依赖第三方求解器或统计推断。

```
prepquant/
  rational.py   # 有理数解析/格式化（拒绝 float）
  expr.py       # 线性表达式 const + Σ a_i·x_i
  lp.py         # 精确两阶段单纯形（Bland 规则）
  ledger.py     # 配液账：有序原子批、来源账本、产品版本、更正重放
  vault.py      # 定量仓：测量导入、联合线性区间约束、精确端点
  scenarios.py  # 合成场景生成器（含手算真值）
  demo.py       # 端到端演示脚本（写 out/）
  cli.py        # 命令行
tests/          # pytest（30 项）
```

## 运行

```bash
python3 -m prepquant.demo          # 全流程演示 + 手算核对，写 out/demo_report.json
python3 -m pytest -q               # 单元/CLI 测试
python3 -m prepquant.scenarios --out out/scenarios.json   # 仅生成场景
```

## 模型与表示

* 体积为正常数 mL，溶质质量为非负一次式，浓度 = 质量/体积，均精确有理。
* 每个容器保存 `mass = const + Σ a_i·x_i`（mg），**x_i 是各原样的未知浓度**。
  转移/混合按体积比例搬运一次式，蒸发只改体积不改质量，定容只加零质量溶剂，
  废弃清零。因此浓度回算是精确的，混合来源与中间移除都不会算错。
* 每个容器另有**来源账本**：`origin -> (仍贡献的体积, 仍贡献的质量一次式)`。
  转移按体积占比分摊来源；混合合并账本；蒸发按比例缩放体积（质量表达式不变）；
  废弃清空。`explain` 可单独交付来源与余额，不依赖任何测量。
* 标准液（已知 const）与原样（x_i）同为来源，可任意混合；系数为零的测量
  仍然参与约束（可直接导致矛盾）。

### 表达式 / 重放 / 求解方法的取舍及影响

| 选择 | 影响 |
| --- | --- |
| 以**质量**为线性货币而非浓度 | 每次操作都是一次线性运算，无需在中间步骤做除法；浓度只在查询时算一次 |
| 容器内系数**稠密**（全部已声明原样） | 实现简单、`rebase` 可在后续新增原样时对齐；导出时**只导出非零原样**，避免为无关样品编造上下界 |
| **重放**而非增量修改 | 更正 = 从目标步之前的状态重放后续全部步骤；失败则整批不落库。代价是更正 O(n)，但题面规模（≤80 容器/≤20 测量）可忽略 |
| 产品 id = `P<seq>-<容器>-<内容哈希>` | 来源与读数绑定到**内容**。未消费时可更正并让旧 id 消失；已消费时更正被拒绝，旧分析快照原样保留 |
| 精确两阶段单纯形 + Bland | 端点与可行性可被本题算术逐行复核；无需浮点容差。规模小，耗时 <1 s |
| 每次分析对所有测量**联合**求解 | 共享原样的不同混合不会被拆开解；多解时只报告精确区间，绝不取均值当作真值 |

## 命令

```bash
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db init
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db prep ops.json     # 有序原子批
python3 -m prepquant.cli --ledger lab.db explain M [--json]                     # 体积/质量式/来源
python3 -m prepquant.cli --ledger lab.db export M --out m.json                  # 产品导出
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db \
        measure --name r1 --product-file m.json --observed 28 --uncertainty 0    # 测量导入
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db \
        analyze a1 --measurement r1 --bound s1:0:100                            # 求解
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db show-analysis a1  # 历史结果
python3 -m prepquant.cli --ledger lab.db correct 3 fix.json --note "wrong top-up"
python3 -m prepquant.cli --ledger lab.db log | corrections | consumed | containers
```

操作 JSON 为对象列表，`op` 取 `standard / sample / transfer / mix / solvent /
evaporate / discard`。数值用 `"12"`、`"1/3"`、`"0.25"` 这类字符串（拒绝二进制
float）。示例：

```json
[
 {"op":"standard","container":"STD","volume":"10","concentration":"100"},
 {"op":"sample","container":"RAW","volume":"10","sample":"s1"},
 {"op":"transfer","source":"RAW","volume":"2","destination":"D"},
 {"op":"solvent","container":"D","volume":"8"},
 {"op":"mix","sources":[{"container":"D","volume":"5"},{"container":"STD","volume":"1"}],"destination":"M"},
 {"op":"evaporate","container":"M","to":"4"}
]
```

产品导出格式（定量仓只认这个，不读配液私有表；`product_id` 可省略）：

```json
{"product_id":"P6-M-…", "container":"M", "volume":"4", "unknowns":["s1"],
 "mass":{"const":"100","terms":{"s1":"1"}},
 "origins":[{"origin":"STD","volume":"1","mass":{"const":"100","terms":{}}},
            {"origin":"RAW","volume":"2/3","mass":{"const":"0","terms":{"s1":"1"}}}]}
```

## 规则与不变量

* 建立标准液/原样必须落在空容器上；容器身份稳定，已装液容器不能再被
  `standard`/`sample`/`transfer`/`mix` 当作新容器写入。
* 取出体积必须 ≤ 真实剩余体积；同一份液不能供两次；`mix` 中同一容器只可出现一次。
* 蒸发后溶质不变且体积必须仍为正、必须严格减小。
* 任一步失败 → 整批回滚（`apply_batch` 先全量重放验证，成功才写库）。
* 更正未消费步骤：从其之前状态重放；已消费产品 id 涉及的更正被拒绝，
  旧分析结果与表达式快照永不改写。
* 测量 `observed ± uncertainty` ⇒ 真浓度 ∈ `[max(0,o-u), o+u]`；乘以产品体积
  后变成质量上的两条精确线性行。分析对所有测量联立，并叠加各原样给定的
  非负上下界。输出：`feasible`/`contradiction`、各原样精确 `[min,max]`、
  是否唯一（所有原样 min==max 才标记唯一）。多解只给区间和端点处的其他
  原样取值，不取平均。
* 外部产品只要符合上述导出格式即可被定量仓直接导入。

## 合成场景与手算核对

`prepquant.scenarios` 生成 6 个场景，`demo.py` 逐一核对质量式与最终区间：

1. **两级稀释 + 标准加入 + 蒸发**：`B` 中 `100 + s1` mg / 4 mL，读数 28 ⟹ `s1=12`。
2. **共享原样的不同混合**：`M1`、`M3` 两条读数联合 ⟹ `s1=15`（联合而非分别解）。
3. **区间多解**：读数 `6±1` 于 `s1/5` ⟹ `s1∈[25,35]`，标记为非唯一。
4. **不相容读数**：同一产品 `6` 与 `7` ⟹ 矛盾。
5. **两个未知量唯一解**：两条独立读数 ⟹ `s1=10, s2=40`，均唯一。
6. **零系数矛盾**：只有标准液的产品读数 50，任何原样取值都无法满足 ⟹ 矛盾。

另在同一构建里验证：取液不足整批回滚、蒸发保持质量、更正后被消费分支拒绝、
未消费更正重放、重开数据库状态一致。

`out/demo_report.json` 是运行产物，`tests/` 中还包含：
`test_lp.py`（与 3 维暴力顶点枚举对比）、`test_ledger.py`、`test_vault.py`、
`test_cli.py`。

## 精确性保证

* 输入禁止 `float`；内部只用 `Fraction`。
* 单纯形返回的每个顶点都会用原始约束逐行 `Fraction` 复核（`_verify`），
  不一致即抛错而不是给出近似值。
* `tests/test_lp.py` 在二维/三维随机系统上把可行性、每个变量的最小/最大
  端点与暴力枚举全部顶点（Cramer 法则）比对。
