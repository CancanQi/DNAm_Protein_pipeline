# DNAm_Protein_pipeline

## 研究思路
* EWAS 识别与表型相关的CpG位点
* pQTM 识别与CpG位点相关的proteins
* 中介 识别表型-CpG-protein关联通路
* SMR/coloc 证实表型-CpG-protein间的causal通路
  
## 参考文献
https://www.nature.com/articles/s41467-022-32319-8 \
https://www.nature.com/articles/s41467-025-57288-6

## 参考代码
### 1. EWAS
* 使用线性/逻辑回归识别相关CpG位点，[R语言参考代码](./1.EWAS/EWAS.demo.R) 注意根据结局变量的类型调整function中的模型
* 


