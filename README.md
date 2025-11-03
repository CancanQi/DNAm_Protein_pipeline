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
* 使用线性/逻辑回归识别相关CpG位点 \
  [R语言参考代码](./1.EWAS/EWAS.demo.R)，注意根据结局变量的类型调整function中的模型
* 提交任务至集群运行 \
  [任务提交参考代码](./1.EWAS/Job.submit.sh)，注意调整具体参数
* CpG位点的注释 \
  注释文件：https://github.com/zhou-lab/InfiniumAnnotationV1/tree/main/Anno/MSA

### 2. pQTM
* 使用线性模型\相关性分析等方法识别CpG-protein之间的关联,[R语言参考代码](./2.pQTM/2.pQTM.R) \
  注意需要先准备蛋白对应的编码[基因的注释](./2.pQTM/1.Annotation_gene.R) \
* 使用MatrixEqtl包速度更快，参考https://github.com/andreyshabalin/MatrixEQTL \

### 3. 中介分析
* 首先，使用使用线性模型\相关性分析等方法，找到上述pQTM蛋白与表型之间的关联 \
* 使用mediation包进行中介分析，[参考代码](./3.Mediation/mediation.R) \

### 4. SMR/coloc
  


