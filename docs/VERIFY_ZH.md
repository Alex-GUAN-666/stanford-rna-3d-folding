# 安装与验证

[English guide](VERIFY.md) · [完整命令与输入格式](REPRODUCIBILITY.md)

这份指南介绍如何在新环境中安装 CPU 工具、运行合成示例和演示 Notebook。
不需要 GPU 或模型权重。我们的比赛成绩和相关记录见[比赛结果](RESULTS.md)。

## 1. 下载并安装

先安装 Git 和 **Python 3.10 或以上版本**。下载仓库和安装依赖需要联网。

### Windows PowerShell

直接使用虚拟环境中的 Python，无需运行激活脚本：

```powershell
git clone https://github.com/Alex-GUAN-666/stanford-rna-3d-folding.git
cd stanford-rna-3d-folding
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install .
.\.venv\Scripts\python.exe scripts/verify_repo.py --output outputs/verification.json
```

如果没有 `py` 命令，先检查 `python --version`，再使用
`python -m venv .venv` 创建环境。

### macOS / Linux

```bash
git clone https://github.com/Alex-GUAN-666/stanford-rna-3d-folding.git
cd stanford-rna-3d-folding
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install .
.venv/bin/python scripts/verify_repo.py --output outputs/verification.json
```

如果 Linux 提示缺少 `venv`，先安装发行版对应的 Python 虚拟环境组件。
修改包内代码或切换提交版本后，在同一环境中再次执行 `pip install .`，然后验证。

## 2. 查看结果

验证脚本会运行单元测试，并在仓库目录之外调用安装好的命令行工具。检查范围包括
长度规划、序列统计、模型输入导出、候选组装、提交格式，以及导入后重新组装的数值核对。

| 检查项 | 预期结果 |
| --- | --- |
| 单元测试 | 45 项通过 |
| 命令检查 | 11 项通过 |
| 示例序列 | 3 条合成序列，长度为 72、160、520 |
| 提交文件 | 752 行、18 列，每条序列 5 组坐标 |
| 导入与重新组装 | 按目标和残基 ID 核对，坐标值保持一致 |

成功时退出码为 0，`outputs/verification.json` 中的 `passed` 为 `true`。
报告包含各条命令的输出、运行环境和源码哈希，也会检查安装的包是否与下载的源码
一致。失败时退出码非零，报告中会记录 `error`。

每次验证使用临时示例数据。如果想保留一份示例供自己查看，在仓库根目录运行：

```powershell
.\.venv\Scripts\python.exe -m rna_folding demo --output-dir outputs/my_demo
.\.venv\Scripts\python.exe -m rna_folding validate --sequences outputs/my_demo/sequences.csv --submission outputs/my_demo/submission.csv
```

macOS / Linux 用户将 `.\.venv\Scripts\python.exe` 替换成 `.venv/bin/python`。
随后查看 `outputs/my_demo` 中的 `submission.csv`、`plan.json` 和 `audit.json`，
分别了解坐标输出、按长度规划的结果和候选结构选择。再次运行时使用新目录或空目录。

## 3. 运行 Notebook

安装可选依赖，并启用 Notebook 检查：

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-notebooks.txt
.\.venv\Scripts\python.exe scripts/verify_repo.py --notebook --output outputs/verification_notebook.json
```

macOS / Linux 同样替换 Python 路径。此选项在 Jupyter 内核中执行演示 Notebook 的
**7 个代码单元格**，并在报告旁保存 `executed_walkthrough.ipynb`，包含输出和合成几何图。

如果想自己逐格运行，可以启动 JupyterLab：

```powershell
.\.venv\Scripts\python.exe -m jupyterlab
```

打开 `notebooks/01_workflow_walkthrough.ipynb`，选择安装了本项目的 Python 环境，
按顺序执行即可。

## 4. 对照 GitHub 自动检查

[GitHub Actions](https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions)
在 Ubuntu、Windows 上分别检查 Python 3.10 和 3.12；Ubuntu / Python 3.12 还会执行
演示 Notebook。选择与你使用的提交版本对应的记录，查看日志或下载验证报告。

反馈问题时，可以用下面的命令取得代码版本：

```bash
git rev-parse HEAD
```

附上这串提交 ID 和验证 JSON，便于定位问题。已有运行记录见
[VALIDATION.md](VALIDATION.md)。

## 检查范围

以上检查覆盖 CPU 工具和合成几何示例，不会生成神经网络预测或重新计算比赛成绩。
历史 GPU Notebook 需要另外配置模型环境和资源，依赖与已知问题见
[历史实现笔记](HISTORICAL_NOTES.md)。

可选的 `evaluate` 命令还需要单独安装 USalign，并提供真实预测结构与参考结构。
自动测试使用替代程序检查调用和解析逻辑，不会执行完整的比赛评估。
