# Stanford RNA 3D Folding：项目说明

**Team GZSL｜76 / 1,516｜最终私榜 0.42295｜铜牌**

我是 Alex Guan（Yuzhen Guan）。这里整理了我们参加 2025 年 Stanford RNA 3D Folding
比赛时的推理实验，以及我赛后补充的坐标处理和提交验证工具。我在队内主要参与
协调、按序列长度比较实验结果，以及模型输出的整合。

[技术报告 PDF](evidence/solution_notes.pdf) · [获奖证书](evidence/certificate.jpg)
· [榜单截图](evidence/private_leaderboard.png) · [本机运行指南](docs/VERIFY_ZH.md)

## 这个项目做什么

比赛要求根据 RNA 序列预测三维结构。长度为 `L` 的序列需要提交五个候选结构，
每个结构包含 `L` 个 C1′ 原子的三维坐标。

我们尝试了 Protenix、DRfold2 和 trRNA 相关的预训练模型，重点比较不同长度下的
推理设置和候选组合。保留下来的代码包括 Protenix + trRNA notebook 和一份独立的
DRfold2 notebook，它们记录的是不同实验，没有合成一份最终三模型运行入口。

我把历史备份放在 [`historical/`](historical/README.md)，赛后整理的工具放在
[`rna_folding/`](rna_folding)。当前工具接收已有预测坐标，完成残基检查、重叠片段
对齐、五候选选择和提交导出。它可以在普通电脑上运行；完整重跑当年的模型推理，
还需要对应的外部源码、权重和最终提交文件。

## 先运行一个例子

下载仓库后，在项目根目录使用 Python 3.10 或以上版本，建议先创建虚拟环境：

```bash
python -m pip install .
python -m rna_folding demo --output-dir outputs/demo
python -m rna_folding validate --sequences outputs/demo/sequences.csv --submission outputs/demo/submission.csv
python scripts/verify_repo.py --output outputs/verification.json
```

`outputs/demo/` 中会生成三个测试目标、候选坐标、长度分配计划、提交 CSV 和选择记录。
每次演示请使用新的或空的输出目录。演示坐标是人工生成的几何样例，用来检查代码，
没有运行神经网络，也不代表真实 RNA 预测精度。

最后一条命令会执行测试、检查安装后的代码，以及导入和重新导出的数值一致性。
Windows、macOS 和 Linux 的环境设置步骤见 [本机运行指南](docs/VERIFY_ZH.md)；
[演示 notebook](notebooks/01_workflow_walkthrough.ipynb) 展示了坐标处理过程。

## 成绩与实验记录

| 内容 | 结果 |
| --- | --- |
| 最终 private leaderboard | **0.42295** |
| Team GZSL 最终排名 | **76 / 1,516，铜牌** |
| 后续测试时的 public leaderboard | **0.409** |

提交时我们只能看到公榜分数，最终私榜公布后得到的成绩是 **0.42295**。
公榜和私榜对应不同评估范围，所以不能把 `0.409 → 0.42295` 解释成同一测试集上的提升。
早期报告中的其他分数及其上下文保留在 [成绩记录](docs/RESULTS.md)；
`0.381.ipynb` 只是备份文件名，没有找到与之对应的评分记录。

历史 notebook 的配置、保存的输出以及已知问题见 [历史实验说明](docs/HISTORICAL_NOTES.md)。
当前 CPU 工具没有重新跑出比赛分数。

## 接入自己的预测

如果已有五候选 `submission.csv` 和对应的序列 CSV，可以导入为候选数组：

```bash
python -m examples.import_submission --sequences data/test_sequences.csv --submission data/recovered_submission.csv --model historical_mixed --output-dir outputs/recovered
```

程序检查目标和残基 ID，将坐标导出为 `.npy` 与候选清单。`historical_mixed` 用于
每个候选来源尚未区分的混合提交。输入格式、组装命令和外部模型接入要求见
[运行说明](docs/REPRODUCIBILITY.md)，片段对齐与候选选择见 [技术方法](docs/METHOD.md)。

比赛成绩由团队共同取得，基础模型与公开实现来自各自的上游作者。项目整理过程见
[项目历史](docs/PROVENANCE.md)，相关工作及使用说明见 [参考资料](docs/REFERENCES.md)
和 [NOTICE](NOTICE.md)。
