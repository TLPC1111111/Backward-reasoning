# AB1 File Analyzer - Sequencing Analysis Reverse Engineering Tool

A comprehensive Python tool for deep analysis and comparison of AB1 (ABIF format) sequencing files to reverse-engineer the data enhancement algorithms used by Applied Biosystems Sequencing Analysis 5.2.0 software.

## 概述 (Overview)

本工具通过对比原始和增强的 `.ab1` 文件，深入分析并推导 Sequencing Analysis 5.2.0 使用的数据增强和优化算法。

This tool performs in-depth analysis by comparing original and enhanced `.ab1` files to reverse-engineer the data enhancement and optimization algorithms used by Sequencing Analysis 5.2.0.

## 功能特性 (Features)

### 1. AB1文件解析 (AB1 File Parsing)
- 完整解析 ABIF (Applied Biosystems 格式) 文件结构
- 提取关键数据字段:
  - `DATA9-12`: A, T, G, C 四个通道的原始信号 trace 数据
  - `PLOC1/2`: 峰位置数据
  - `PBAS1/2`: Base calling 结果
  - `PCON1/2`: Quality scores (Phred分数)
  - `FWO_1`: 染料顺序
  - 元数据（机器型号、运行参数等）

### 2. 对比分析功能 (Comparative Analysis)

#### a) 信号处理分析 (Signal Processing Analysis)
- 四通道信号对比绘制
- 信号变化模式计算（平滑、基线校正、峰值增强）
- 滤波器类型检测（移动平均、Savitzky-Golay、高斯滤波）
- 信噪比（SNR）改善分析

#### b) 峰检测算法推断 (Peak Detection Algorithm Inference)
- 峰位置变化对比
- 峰高度/宽度/形状调整分析
- 峰分离/合并策略识别

#### c) Base Calling优化 (Base Calling Optimization)
- 碱基调用结果差异对比
- 质量分数（Phred score）变化分析
- 模糊碱基（N）处理策略识别

#### d) 文件结构分析 (File Structure Analysis)
- 文件大小变化分析
- 数据压缩或裁剪策略检测
- 元数据修改分析

### 3. 可视化 (Visualizations)
自动生成以下对比图表:
- `trace_comparison.png`: 四通道信号对比
- `peak_comparison.png`: 峰位置和高度对比
- `quality_comparison.png`: 质量分数分布对比
- `snr_comparison.png`: 信噪比对比

### 4. 详细报告 (Comprehensive Report)
生成包含完整分析结果的 Markdown 报告，包括:
- 文件信息和大小变化
- 信号处理算法推断
- 峰检测优化策略
- Base calling 改进分析
- 算法复现建议

## 安装 (Installation)

### 环境要求 (Requirements)
- Python 3.8+
- 依赖库在 `requirements.txt` 中列出

### 安装步骤 (Installation Steps)

```bash
# 克隆仓库
git clone https://github.com/TLPC1111111/Backward-reasoning.git
cd Backward-reasoning

# 安装依赖
pip install -r requirements.txt
```

## 使用方法 (Usage)

### 基本用法 (Basic Usage)

```bash
python ab1_analyzer.py \
  --original Origin/pGem-PC_Hongene_2025-09-05_E09.ab1 \
  --augmented Software_augmented/pGem-PC_Hongene_2025-09-05_E09.ab1 \
  --output analysis_report.md
```

### 参数说明 (Arguments)

- `--original`: 原始 AB1 文件路径 (Required)
- `--augmented`: 增强 AB1 文件路径 (Required)
- `--output`: 输出报告文件名 (默认: `analysis_report.md`)
- `--output-dir`: 输出目录路径 (默认: 当前目录)

### 示例 (Example)

```bash
# 分析示例文件对
python ab1_analyzer.py \
  --original Origin/pGem-PC_Hongene_2025-09-05_E09.ab1 \
  --augmented Software_augmented/pGem-PC_Hongene_2025-09-05_E09.ab1 \
  --output-dir ./results
```

## 输出文件 (Output Files)

运行后将生成以下文件:

```
.
├── analysis_report.md          # 详细分析报告
└── visualizations/              # 可视化图表目录
    ├── trace_comparison.png     # 四通道信号对比图
    ├── peak_comparison.png      # 峰位置对比图
    ├── quality_comparison.png   # 质量分数分布图
    └── snr_comparison.png       # 信噪比对比图
```

## 分析结果示例 (Example Results)

### 文件大小变化
- 原始文件: 295,336 bytes
- 增强文件: 288,022 bytes
- 减少: 7,314 bytes (2.48%)

### SNR改善
- G通道: +9.7%
- A通道: +11.8%
- T通道: +6.0%
- C通道: +12.5%
- **平均改善: 10.0%**

### 质量分数提升
- 原始平均质量分数: 47.20
- 增强平均质量分数: 50.31
- **提升: 3.11 points**

### Base Calling优化
- 碱基调用差异率: 6.04%
- 模糊碱基(N)减少: 已优化

## 技术细节 (Technical Details)

### 核心算法 (Core Algorithms)

1. **ABIF文件解析**: 使用 BioPython 的 SeqIO 模块
2. **信号处理分析**: 
   - 移动平均滤波器检测
   - Savitzky-Golay滤波器检测
   - 高斯滤波器检测
3. **基线校正**: 百分位数法
4. **SNR计算**: 峰值功率 / 噪声标准差
5. **相关性分析**: Pearson相关系数

### 依赖库 (Dependencies)

- `numpy`: 数值计算
- `matplotlib`: 数据可视化
- `scipy`: 信号处理算法
- `biopython`: ABIF文件解析

## 算法复现建议 (Algorithm Replication Guidelines)

基于分析结果，可以复现以下算法:

### 1. 信号预处理
```python
import numpy as np
from scipy.signal import savgol_filter

# 应用Savitzky-Golay滤波器
filtered = savgol_filter(signal, window_length=7, polyorder=2)

# 基线校正
baseline = np.percentile(filtered, 10)
corrected = filtered - baseline
```

### 2. 峰检测优化
```python
from scipy.signal import find_peaks

# 自适应阈值峰检测
threshold = np.mean(signal) + 2 * np.std(signal)
peaks, properties = find_peaks(signal, height=threshold, distance=5)
```

### 3. 质量分数重新校准
```python
# 基于信号强度和峰清晰度重新计算Phred分数
def recalculate_quality(peak_height, noise_level):
    snr = peak_height / noise_level
    phred = -10 * np.log10(1 / (1 + snr))
    return int(phred)
```

## 项目结构 (Project Structure)

```
Backward-reasoning/
├── Origin/                      # 原始AB1文件目录
│   └── pGem-PC_Hongene_2025-09-05_E09.ab1
├── Software_augmented/          # 增强AB1文件目录
│   └── pGem-PC_Hongene_2025-09-05_E09.ab1
├── ab1_analyzer.py             # 主分析脚本
├── requirements.txt            # Python依赖列表
├── README.md                   # 本文档
├── analysis_report.md          # 生成的分析报告
└── visualizations/             # 生成的可视化图表
    ├── trace_comparison.png
    ├── peak_comparison.png
    ├── quality_comparison.png
    └── snr_comparison.png
```

## 常见问题 (FAQ)

### Q: 为什么增强文件比原始文件小?
A: 增强过程中可能进行了以下优化:
- 移除冗余数据字段
- 降低数值精度
- 裁剪低质量区域数据

### Q: 相关系数较低是什么原因?
A: 这是正常现象。原始和增强文件的信号长度不同（12631 vs 14054），且经过了不同的处理流程。相关系数主要用于识别最可能的滤波器类型，而不是表示整体相似度。

### Q: 如何处理其他AB1文件对?
A: 只需替换 `--original` 和 `--augmented` 参数的文件路径即可:
```bash
python ab1_analyzer.py \
  --original Origin/your_file.ab1 \
  --augmented Software_augmented/your_file.ab1
```

## 贡献 (Contributing)

欢迎提交问题报告和改进建议！

## 许可证 (License)

本项目用于研究和教育目的。

## 联系方式 (Contact)

如有问题或建议，请在 GitHub Issues 中提出。

## 致谢 (Acknowledgments)

- BioPython 项目提供的 ABIF 文件解析支持
- SciPy 和 NumPy 提供的科学计算工具
- Applied Biosystems 的测序技术

---

**注意**: 本工具仅用于研究和学习目的，用于理解测序数据处理算法。请尊重软件许可和知识产权。
