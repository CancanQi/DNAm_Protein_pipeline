# 注意事项
### 1. 数据使用
* 需要去除X, Y染色体上的探针
* 去除质控不合格的探针
* 结果显著的探针需要通过画图double check

### 2. 模型的校正
* 性别，年龄
* 吸烟情况（如吸烟缺失值多，则使用epigenetic smoking score替代） \
  https://htmlpreview.github.io/?https://github.com/sailalithabollepalli/EpiSmokEr/blob/master/vignettes/epismoker.html
* 细胞类型比例估算 \
  https://www.bioconductor.org/packages/release/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html
* 注意计算lambda值观察是否有inflation，再继续调整模型
