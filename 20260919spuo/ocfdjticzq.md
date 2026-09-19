# SeuratTutorial使用Seurat分析多模态数据

> 更新时间：2026-09-19 (UTC+8)

**写在前面**

学习一个软件最好的方法就是啃它的官方文档。本着自己学习、分享他人的态度，分享官方文档的中文教程。软件可能随时更新，建议配合官方文档一起阅读。推荐先按顺序阅读往期内容：

文献篇：

1.文献阅读：(Seurat V1) 单细胞基因表达数据的空间重建 

2.文献阅读：(Seurat V2) 整合跨越不同条件、技术、物种的单细胞转录组数据 

3.文献阅读：(Seurat V3) 单细胞数据综合整合 

4.文献阅读：(Seurat V4) 整合分析多模态单细胞数据 

5.文献阅读：(Seurat V5) 用于集成、多模态和可扩展单细胞分析的字典学习 

教程篇：

1.Seurat Tutorial 1：常见分析工作流程，基于 PBMC 3K 数据集
目录
1 导入数据

2 设置 Seurat 对象，添加 RNA 和蛋白质数据

3 根据 scRNA-seq 谱对细胞进行聚类

4 并排可视化多种模态

5 识别 scRNA-seq clusters 的细胞表面标记

6 多模态数据的附加可视化

7 从 10X multi-modal 实验加载数据

8 Seurat 中多模态数据的附加功能

官网教程：
1 导入数据
同时测量同一细胞的多种数据类型的能力，称为多模态分析（multimodal analysis），代表了单细胞基因组学的一个令人兴奋的新前沿。例如，CITE-seq 能够同时测量同一细胞的转录组和细胞表面蛋白。其他令人兴奋的多模态技术，例如 10x multiome kit，可以对细胞转录组和染色质可及性进行配对测量（i.e scRNA-seq+scATAC-seq）。其他可以与细胞转录组一起测量的模式包括遗传扰动、细胞甲基化和来自细胞哈希的标签寡核苷酸。我们设计 Seurat4 来实现各种多模态单细胞数据集的无缝存储、分析和探索。

在此小节中，我们介绍了创建多模态 Seurat 对象并执行初始分析的介绍性工作流程。例如，我们演示了如何根据测量的细胞转录组对 CITE-seq 数据集进行聚类，并随后发现每个聚类中富集的细胞表面蛋白。我们注意到 Seurat4 还支持更先进的技术来分析多模态数据，特别是我们的 Weighted Nearest Neighbors (WNN) 方法的应用，该方法能够基于两种模态的加权组合同时对细胞进行聚类。

在这里，我们分析了 8,617 个脐带血单核细胞 (CBMCs) 的数据集，其中转录组测量与 11 种表面蛋白的丰度估计值配对，其水平通过 DNA-barcoded antibodies 进行量化。首先，我们加载两个计数矩阵：一个用于 RNA 测量，另一个用于 antibody-derived tags (ADT)。您可以在官网下载 ADT 文件，以及 RNA 文件。

library(Seurat)library(ggplot2)library(patchwork)
# 导入 RNA UMI 矩阵# 请注意，该数据集还包含约 5% 的小鼠细胞，我们可以将其用作蛋白质测量的负样本对照。# 因此，基因表达矩阵有 HUMAN_ 或 MOUSE_ 附加到每个基因的开头。cbmc.rna - as.sparse(read.csv( file = "data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", header = TRUE, row.names = 1))# 为了让以后的生活更轻松一些，我们将丢弃除前 100 个高表达的小鼠基因之外的所有基因，并从 CITE-seq 前缀中删除 “HUMAN_”。cbmc.rna - CollapseSpeciesExpressionMatrix(cbmc.rna)# 导入 ADT UMI 矩阵cbmc.adt - as.sparse(read.csv( file = "data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", header = TRUE, row.names = 1))# 请注意，由于测量是在相同的细胞中进行的，因此两个矩阵具有相同的列名称all.equal(colnames(cbmc.rna), colnames(cbmc.adt))
2 设置 Seurat 对象，添加 RNA 和蛋白质数据

现在我们创建一个 Seurat 对象，并添加 ADT 数据作为第二个 assay

# 基于 scRNA-seq 数据创建 Seurat 对象cbmc - CreateSeuratObject(counts = cbmc.rna)# 我们可以看到，默认情况下，cbmc 对象包含一个存储 RNA 测量值的 assayAssays(cbmc)# 创建一个新的 assay 来存储 ADT 信息adt_assay - CreateAssayObject(counts = cbmc.adt)# 将此 assay 添加到之前创建的 Seurat 对象中cbmc[["ADT"]] - adt_assay# 验证对象现在包含多个 assaysAssays(cbmc)# 提取 ADT assay 中测量的特征列表rownames(cbmc[["ADT"]])# 请注意，我们可以轻松地在两种检测之间来回切换以指定默认值用于可视化和分析# 列出当前默认 assayDefaultAssay(cbmc)# 将默认值切换为 ADTDefaultAssay(cbmc) - "ADT"DefaultAssay(cbmc)
3 根据 scRNA-seq 谱对细胞进行聚类

以下步骤代表基于 scRNA-seq 数据的 PBMCs 快速聚类。有关各个步骤或更高级选项的更多详细信息，请参阅此处的 PBMC 聚类指导教程。

# 请注意，以下所有操作均在 RNA assay 上设置上，并验证默认 assay 是 RNADefaultAssay(cbmc) - "RNA"DefaultAssay(cbmc)
# 执行可视化和聚类步骤cbmc - NormalizeData(cbmc)cbmc - FindVariableFeatures(cbmc)cbmc - ScaleData(cbmc)cbmc - RunPCA(cbmc, verbose = FALSE)cbmc - FindNeighbors(cbmc, dims = 1:30)cbmc - FindClusters(cbmc, resolution = 0.8, verbose = FALSE)cbmc - RunUMAP(cbmc, dims = 1:30)DimPlot(cbmc, label = TRUE)
4 并排可视化多种模态

现在我们已经从 scRNA-seq profiles 中获得了 clusters，我们可以可视化数据集中蛋白质或 RNA 分子的表达。重要的是，Seurat 提供了几种在模态之间切换的方法，并指定您有兴趣分析或可视化的模态。这一点特别重要，因为在某些情况下，相同的特征可以以多种方式出现 - 例如，该数据集包含 B cell marker CD19的独立测量结果（both protein and RNA levels）。

# Normalize ADT data,DefaultAssay(cbmc) - "ADT"cbmc - NormalizeData(cbmc, normalization.method = "CLR", margin = 2)DefaultAssay(cbmc) - "RNA"# 请注意，以下命令是替代命令，但返回相同的结果# Note that the following command is an alternative but returns the same resultcbmc - NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")# 现在，我们将可视化 RNA 和 protein 水平的 CD14 通过设置默认 assay，我们可以可视化其中一个或另一个DefaultAssay(cbmc) - "ADT"p1 - FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")DefaultAssay(cbmc) - "RNA"p2 - FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")# place plots side-by-sidep1 | p2
# 或者，我们可以使用特定的 assay key 来指定特定的模态识别密钥用于 RNA 和蛋白质 assaysKey(cbmc[["RNA"]])Key(cbmc[["ADT"]])
# 现在，我们可以在功能名称中包含 key，这会覆盖默认 assayp1 - FeaturePlot(cbmc, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")p2 - FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")p1 | p2
5 识别 scRNA-seq clusters 的细胞表面标记

我们可以利用配对的 CITE-seq 测量来帮助注释来自 scRNA-seq 的 clusters，并识别 protein 和 RNA markers。

# 因为我们知道 CD19 是 一个 B cell marker，我们可以识别 cluster 6 表面表达 CD19VlnPlot(cbmc, "adt_CD19")
# 我们还可以通过差异表达识别该 cluster 的替代 protein 和 RNA markersadt_markers - FindMarkers(cbmc, ident.1 = 6, assay = "ADT")rna_markers - FindMarkers(cbmc, ident.1 = 6, assay = "RNA")head(adt_markers)## p_val avg_log2FC pct.1 pct.2 p_val_adj## CD19 2.067533e-215 1.2787751 1 1 2.687793e-214## CD45RA 8.106076e-109 0.4117172 1 1 1.053790e-107## CD4 1.123162e-107 -0.7255977 1 1 1.460110e-106## CD14 7.212876e-106 -0.5060496 1 1 9.376739e-105## CD3 1.639633e-87 -0.6565471 1 1 2.131523e-86## CD8 1.042859e-17 -0.3001131 1 1 1.355716e-16head(rna_markers)## p_val avg_log2FC pct.1 pct.2 p_val_adj## BANK1 0 1.963277 0.456 0.015 0## CD19 0 1.563124 0.351 0.004 0## CD22 0 1.503809 0.284 0.007 0## CD79A 0 4.177162 0.965 0.045 0## CD79B 0 3.774579 0.944 0.089 0## FCRL1 0 1.188813 0.222 0.002 0
6 多模态数据的附加可视化

# 绘制 ADT 散点图（如 FACS 的双轴图）。请注意，您甚至可以“门控”细胞，如果使用 HoverLocator 和 FeatureLocator 所需的值FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
# 查看 protein 和 RNA 之间的关系FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")
# 让我们看看原始（non-normalized）ADT counts。# 可以看到数值相当高，特别是与 RNA 值相比。# 这是由于细胞中蛋白质拷贝数明显较高，这显着减少了 ADT 数据中的“drop-out”FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")
7 从 10X multi-modal 实验加载数据

Seurat 还能够分析使用 CellRanger v3 处理的多模态 10X 实验的数据；例如，我们使用 7,900 个外周血单核细胞 (PBMC) 的数据集重新创建了上面的图，该数据集可从 10X Genomics 免费获取。

pbmc10k.data - Read10X(data.dir = "../data/pbmc10k/filtered_feature_bc_matrix/")rownames(x = pbmc10k.data[["Antibody Capture"]]) - gsub( pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = pbmc10k.data[["Antibody Capture"]]))pbmc10k - CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 3, min.features = 200)pbmc10k - NormalizeData(pbmc10k)pbmc10k[["ADT"]] - CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])pbmc10k - NormalizeData(pbmc10k, assay = "ADT", normalization.method = "CLR")plot1 - FeatureScatter(pbmc10k, feature1 = "adt_CD19", feature2 = "adt_CD3", pt.size = 1)plot2 - FeatureScatter(pbmc10k, feature1 = "adt_CD4", feature2 = "adt_CD8a", pt.size = 1)plot3 - FeatureScatter(pbmc10k, feature1 = "adt_CD3", feature2 = "CD3E", pt.size = 1)(plot1 + plot2 + plot3) NoLegend()
8 Seurat 中多模态数据的附加功能

Seurat v4 还包括用于分析、可视化和整合多模态数据集的附加功能：

使用 Seurat v4 中的 WNN 分析从多模态数据定义细胞身份

将 scRNA-seq 数据映射到 CITE-seq references

空间转录组学分析简介

使用 WNN 分析进行 10x 多组分析 (paired scRNA-seq + ATAC)

Signac：单细胞染色质数据集的分析、解释和探索

Mixscape：用于汇总单细胞遗传筛选的分析工具包

这些内容将在后续推文中介绍

**结束**

## 相关阅读

- [月经期间可以泡脚吗月经期间的注意事项有哪些](https://github.com/ntyvivo01u/baby-feeding-guide/blob/main/20260911rckp/idhecokkea.md)
- [赛增副作用危害真不小？其实只要3！](https://github.com/uo8lrun64a/pregnancy-care-hub/blob/main/20260910lkfw/xysothhzcp.md)
- [美国第三代试管婴儿终极手册：流程、费用、成功率全解析](https://github.com/sa1ec5y0bz/baby-care-journal/blob/main/20260910aruz/ijtcluursi.md)
- [三代试管婴儿7种促排方法成功率高吗，附成功案例分享！](https://github.com/iebkyzpjrn/pregnancy-care-hub/blob/main/20260910xmph/alacorbapd.md)
- [山西十大试管婴儿医院排行榜来了](https://github.com/bx6ti255zt/family-parenting-notes/blob/main/20260915ophe/bazlhfidcn.md)
- [仁济和九院做试管婴儿哪个比较好！钱花到了哪！](https://github.com/vmlbl9r4m3/toddler-food-ideas/blob/main/20260917kfqh/pehhhgwlhd.md)
- [泰国三代试管婴儿攻略](https://github.com/yoz4ykilda/mother-baby-diary/blob/main/20260911jixj/fvodxjlavf.md)
- [39做试管用拮抗剂方案取卵较少吗？](https://github.com/na1l60kg9l/family-parenting-notes/blob/main/20260915rycd/ybpsqxvzht.md)
- [如何提高高龄女性试管婴儿成功概率？高龄产妇产后如何坐好月子？](https://github.com/ddk2koak3u/child-care-essays/blob/main/20260910zasd/yadisugxih.md)
- [常州三代试管婴儿成功案例纪实？常州三代试管婴儿成功案例纪实？](https://github.com/whprpfn9bc/baby-care-journal/blob/main/20260910nbuy/hsiuiyitlr.md)
- [试管移植后便秘会有哪些影响(听马医生在线科普便秘的危害)](https://github.com/wggadvmpg6/child-care-essays/blob/main/20260916sxbu/wereigrrjs.md)
- [育龄期妇女有高血压，可以怀孕吗？](https://github.com/nnhgjqxjg6/family-life-notes/blob/main/20260911trpz/ovarrvnozr.md)
- [试管胚胎10b质量怎么样？试管婴儿技术移植胚胎的时间怎么决定](https://github.com/rnf9cvz5iw/mother-baby-diary/blob/main/20260915lqgc/fhobklymtm.md)
- [猫咪胎动后多久生产？详解猫咪妊娠期](https://github.com/bnab3b3j5y/pregnancy-diary-hub/blob/main/20260911tivr/fegnzauwtk.md)
- [国内单身试管私立医院哪家好？黑河排名榜单值得参考](https://github.com/rzchuf6kdk/baby-feeding-guide/blob/main/20260918dykr/rbcgirihhp.md)
- [广州三代试管婴儿好的医院是哪一家](https://github.com/cwz1rtzls4/child-care-essays/blob/main/20260916thcm/fihrnhndkc.md)
- [做过三代试管需要多少钱〖三代试管大概要多久〗](https://github.com/ij0s3j0vss/child-care-essays/blob/main/20260916qsdw/nspzkawjnb.md)
- [做试管一个月能不能成功，一个月能完成不可信](https://github.com/bjpnmb0r46/parenting-skills-log/blob/main/20260915lsdn/whnjtpnzrk.md)
- [HIV艾滋病去山西煤炭医院做三代试管婴儿流程参考！](https://github.com/l9lvqnbe4d/baby-growth-journal/blob/main/20260917iriw/itouxhjxvl.md)
- [试管婴儿成功后，保胎攻略全解析](https://github.com/helxwyn5td/infant-nutrition-hub/blob/main/20260918qibv/rnopssspwp.md)
- [赤峰宁城县试管婴儿哪个医院好！千万不要错过！](https://github.com/vedmkiygf6/maternal-care-journal/blob/main/20260915qreh/gbhmzdxvdj.md)
- [中国试管婴儿医院排名及推荐](https://github.com/y9qvvxks1i/family-health-notes/blob/main/20260916wqgh/hamvkdhwjb.md)
- [没有结婚证可以去南昌康健生殖医院做三代试管生子吗？成功率高吗](https://github.com/qws8inv2p1/family-health-notes/blob/main/20260916hpst/mkdynafnru.md)
- [试管促排卵血值高的原因？试管促排血值低是什么原因？](https://github.com/bx6ti255zt/child-education-notes/blob/main/20260911gmww/bewlmyowoj.md)
- [三代试管一次成功多吗？一篇文章详细解释了123代试管婴儿的成功概率](https://github.com/sa1ec5y0bz/mommy-baby-notes/blob/main/20260916hdyx/hapzrpcgek.md)
- [北京花生医疗试管婴儿服务：咨询热线与价格解析，尽](https://github.com/nnhgjqxjg6/family-life-notes/blob/main/20260911trpz/ypnzpcatvy.md)
- [自己从国内去俄罗斯做三代试管多少钱啊(去俄罗斯做试管移植费用)](https://github.com/oizha1rquq/newborn-parenting-log/blob/main/20260915deeh/hicuzfkwfa.md)
- [四川试管婴儿成功率最高的医院可以签约吗！四川能做试管婴儿的医院！](https://github.com/nih9jzz6yi/baby-care-journal/blob/main/20260916uube/uijdrvctvn.md)
- [问题性肌肤有哪几种](https://github.com/achf8mo3od/mommy-care-diary/blob/main/20260911phew/wjbcohkcdy.md)
- [避开费用坑！染色体异常家庭在河池市人民医院做试管的真实花费](https://github.com/ualf0k98cv/pregnancy-care-hub/blob/main/20260910nzzm/qfunaftdqc.md)
- [乌克兰试管婴儿直营机构排行榜(乌克兰试管婴儿直营机构排行榜最新)](https://github.com/mxtw9dwa7v/baby-sleep-tips/blob/main/20260915ahfy/zsrkuealrv.md)
- [广州正规试管机构排名哪家最好？广州市试管婴儿排名？](https://github.com/g6iv5x0e8m/child-care-essays/blob/main/20260916xyfx/zehbmqsqzk.md)
- [怀不上教你三招助孕方法 想怀怀不上要做检查](https://github.com/nnhgjqxjg6/family-life-notes/blob/main/20260917ncfm/qyffyqzmcq.md)
- [泰国私立医院试管成功率一览：详细流程与成功案例解析！](https://github.com/sa1ec5y0bz/pregnancy-care-hub/blob/main/20260910iycg/twbznvouot.md)
- [齐心向党，医路前行](https://github.com/achf8mo3od/newborn-parenting-log/blob/main/20260911hner/pwdmplowih.md)
- [【名医坐诊】8月15日（周四）鱼台县总医院（县人民医院）妇产科专家王洪玲在县妇幼保健院坐诊](https://github.com/znp78by4gt/toddler-activity-ideas/blob/main/20260911rjls/hoovajmzxr.md)
- [只有4个c级婴性精子可以做一代试管吗？](https://github.com/vjd2jnnrxj/infant-nutrition-hub/blob/main/20260911pvbe/dahqkeipcb.md)
- [国内三代试管比较好的私立医院名单汇总](https://github.com/h5z4rt20ta/parenting-daily-tips/blob/main/20260919fani/hiasnabkvw.md)
- [美国单身女性做一次试管婴儿大概要多少费用？](https://github.com/s6nb3rgjk9/baby-care-journal/blob/main/20260910iseb/lqgqkbxhyj.md)
- [苏州试管婴儿价格多少？10万费用贵吗？](https://github.com/mxtw9dwa7v/maternal-care-journal/blob/main/20260915sahx/cuhbfinzbb.md)

## 推荐站点

- [['https://www.syldezdhkj.cn/18771352153535.html', '国外40岁以上试管婴儿借卵费用成功率有多大 国外做试管婴儿借卵费用成功率很高吗']](https://www.syldezdhkj.cn/18771352153535.html)
- [['https://www.skiguo.cn/20250927-179.html', '试管供卵排名-27岁绝经还能恢复吗']](https://www.skiguo.cn/20250927-179.html)
- [['https://hangzhou.ccxwlkx.cn/309.html', '泰国三代试管婴儿费用详解及详细步骤指南']](https://hangzhou.ccxwlkx.cn/309.html)
- [['https://www.uueamru.cn/20250222-89.html', '广州第三代试管：女性取卵会疼吗？有没有什么副作用？']](https://www.uueamru.cn/20250222-89.html)
- [['https://www.cmanrxrr.cn/2775115062531.html', '2026年成都西囡医院试管借卵费用要多少钱？,代孕供卵价格']](https://www.cmanrxrr.cn/2775115062531.html)
- [['https://www.luruihang.com/2333.html', '徐州妇幼保健院泰国代生哪家好成功率高不高']](https://www.luruihang.com/2333.html)
- [['https://www.chengyanghg.cn/336.html', '三代试管婴儿医院费用预算与选择指南']](https://www.chengyanghg.cn/336.html)
- [['https://www.wqxmm.cn/308844489472.html', '【2026实测】河南家圆医院供卵多少钱？费用构成与公立医院对比']](https://www.wqxmm.cn/308844489472.html)
- [['https://www.dygsdyw.com/206661657145.html', '苏州供卵试管套餐，苏州供卵生小孩,苏州市立医院试管婴儿要住院吗试管婴儿需要请假吗']](https://www.dygsdyw.com/206661657145.html)
- [['https://www.satghenga.cn/228171675169.html', '武汉康健妇婴医院怎么联系？武汉康健妇婴医院地址在哪？']](https://www.satghenga.cn/228171675169.html)
- [['https://www.afa2019.com/125943447577.html', '2026年石家庄做三代试管婴儿费用及成功率']](https://www.afa2019.com/125943447577.html)
- [['https://www.sjb493.cn/30996460907122.html', '代生公司官网-供卵代怀网价格表,国内做第三代试管比较厉害的医院大全']](https://www.sjb493.cn/30996460907122.html)
- [['https://www.gzgudadl.cn/2299049882357.html', '二代哪家代生靠谱费用是多少贵吗？要多少费用！']](https://www.gzgudadl.cn/2299049882357.html)
- [['https://www.ppmaas.com/guoneishiguanjigou/184.html', '代孕供卵费用：乌鲁木齐试管婴儿专家有哪些(乌鲁木齐市妇幼保健医院试管医生有哪些)']](https://www.ppmaas.com/guoneishiguanjigou/184.html)
- [['https://www.monpun.com/8309207216965.html', '卵巢早衰诊疗专家指南']](https://www.monpun.com/8309207216965.html)
- [['https://www.bjwdzxkj.cn/2553773130805.html', '代生那儿最权威促排期间腰疼怎么回事']](https://www.bjwdzxkj.cn/2553773130805.html)
- [['https://www.apkbwvg.cn/danshenshiguanfangan/126.html', '二胎备孕科学调理与排卵监测方法']](https://www.apkbwvg.cn/danshenshiguanfangan/126.html)
- [['https://www.hflrwzhs.cn/175.html', 'XY和XX的奥秘：除了XY看性别，染色体里还藏着哪些遗传病密码？']](https://www.hflrwzhs.cn/175.html)
- [['https://www.fmngst.com/1524239148997.html', '2026安徽试管男孩多少钱？试管成功率高吗？']](https://www.fmngst.com/1524239148997.html)
- [['https://www.hg00fj88.com/2099.html', '胚胎移植后会不会掉出来胚胎移植后什么情况会掉出来']](https://www.hg00fj88.com/2099.html)
- [['https://www.mimi567.com/389.html', '在南宁二医院做试管婴儿需要审核结婚证吗？']](https://www.mimi567.com/389.html)
- [['https://www.ewdboe.cn/223805117489.html', '海南借卵试管生男孩医院排名：助孕优选指南与成功经验']](https://www.ewdboe.cn/223805117489.html)
- [['https://www.cndcxc.com/daiyunliucheng/20251021/17058.html', '孕5周胚胎着床了吗']](https://www.cndcxc.com/daiyunliucheng/20251021/17058.html)
- [['https://www.zrbbavaq.cn/12288769724157.html', '上海私立供卵机构名单统计，2026高龄供卵试管代生报价医院指南']](https://www.zrbbavaq.cn/12288769724157.html)
- [['https://www.cd-hssf.com/318570643024.html', '过来人分享性病八项检查多少钱（省钱版）,试管代孕借卵收费吗,哪里有靠谱代孕公司']](https://www.cd-hssf.com/318570643024.html)
- [['https://www.vecsi.cn/shanxizhuyunjiage/2719.html', '做人工授精女人要打针吗？女人人工授精手术需要打麻药吗？']](https://www.vecsi.cn/shanxizhuyunjiage/2719.html)
- [['https://www.xnnpbhdz.cn/12039834425293.html', '西安三甲🏥预约挂号攻略,代孕哪里技术好&国内的供卵机构有哪些']](https://www.xnnpbhdz.cn/12039834425293.html)
- [['https://www.cecigou.cn/daihuaiyunfuwu/20250928/15205.html', '月经干净十天后又有褐色分泌物']](https://www.cecigou.cn/daihuaiyunfuwu/20250928/15205.html)
- [['https://www.sdshunhezb.cn/617471309190.html', '2026年山东供卵代生医院Top 5权威榜单：高成功率与代怀生男孩机会解析']](https://www.sdshunhezb.cn/617471309190.html)
- [['https://www.jszgyh.com/126252483409.html', '南通试管婴儿助孕费用明细，2026助孕全流程总花费预估']](https://www.jszgyh.com/126252483409.html)
- [['https://www.esc45.com/57.html', '腺肌瘤试管婴儿初诊多少钱 子宫腺肌瘤做试管婴儿长方案']](https://www.esc45.com/57.html)
- [['https://www.phetpalace.com/224.html', '卵泡发育不良的症状（卵泡发育不良的治疗）']](https://www.phetpalace.com/224.html)
- [['https://www.chdhaishendq.cn/227460267355.html', '苏州正规助孕公司靠谱吗？试管婴儿生男孩服务解析']](https://www.chdhaishendq.cn/227460267355.html)
- [['https://www.3899234.com/20250927-154.html', '代生服务平台&空腹备孕运动，空腹备孕运动有影响吗']](https://www.3899234.com/20250927-154.html)
- [['https://www.sasksjob.com/416513646063.html', '2026版医保目录执行指南：如何找到合适的生育辅助治疗方案']](https://www.sasksjob.com/416513646063.html)
- [['https://www.bjfhyly.com/1225.html', '国内代孕排名好的公司_代孕供卵机构,囊胚11天试纸颜色很浅有希望吗？是失败']](https://www.bjfhyly.com/1225.html)
- [['https://www.njxxwcr.cn/daishengdaihuaishengzi/156.html', '海口试管攻略：海医附一院生殖中心技术水平、医生团队及费用参考']](https://www.njxxwcr.cn/daishengdaihuaishengzi/156.html)
- [['https://www.weywjei.cn/20250826-179.html', 'AMH 0.02的绝地求生：拦截早衰结局，通过DHEA与中药联合调理方案']](https://www.weywjei.cn/20250826-179.html)
- [['https://www.anyhdlyb.cn/2528104630377.html', '有没有找人代生&靠谱代生机构,不孕症的分型治疗方法（女性不孕不育预防）']](https://www.anyhdlyb.cn/2528104630377.html)
- [['https://www.sandwnot.com/104782313433.html', '辽宁三甲试管婴儿医院Top10排行-辽宁省试管婴儿的医院排名！']](https://www.sandwnot.com/104782313433.html)
- [['https://www.toothree006.cn/126623122267.html', '有腺肌症做试管后更痛了']](https://www.toothree006.cn/126623122267.html)
- [['https://www.cddyunw.com/502503519383.html', '辅助生殖技术医保可报销额度与流程解析']](https://www.cddyunw.com/502503519383.html)
- [['https://www.zhangruiqing.cn/103631920301.html', '供卵代怀费用-代生在线咨询,试管婴儿做第三代多少钱-试管婴儿第三代费用大约多少']](https://www.zhangruiqing.cn/103631920301.html)
- [['https://www.szanguangkeji.cn/tongxingshiguanzhuyun/136.html', '世纪供卵试管公司排名：解析行业巨头在助孕流程风控上的优势']](https://www.szanguangkeji.cn/tongxingshiguanzhuyun/136.html)
- [['https://www.hs52.cc/tesefuwu/68.html', 'tsh高需要中断促排卵针吗？tsh偏高能打疫苗吗？']](https://www.hs52.cc/tesefuwu/68.html)
- [['https://www.bjjinyukechuangzdh.cn/219.html', '闭经后梦见月经来潮']](https://www.bjjinyukechuangzdh.cn/219.html)
- [['https://www.gaodunxinkj.cn/20250509-168.html', '三代试管婴儿的好处']](https://www.gaodunxinkj.cn/20250509-168.html)
- [['https://www.haojiezhishi.cn/111.html', '试管婴儿不成功能否退款']](https://www.haojiezhishi.cn/111.html)
- [['https://www.dymgp.com/7930.html', '最好试管代怀-促排卵药副作用']](https://www.dymgp.com/7930.html)
- [['https://www.jzcwjz.net/237.html', '毕节试管婴儿全下来多少钱,贵州试管婴儿费用']](https://www.jzcwjz.net/237.html)
- [['https://www.hghbjm.com/252.html', '做试管内膜薄移植成功率高吗？子宫内膜薄试管移植一定不能成功吗？']](https://www.hghbjm.com/252.html)
- [['https://www.dgshengxigongchengsl.cn/1773851428741.html', '第六章-6 有了麻醉，取卵没有想象中的可怕,2026代孕网站']](https://www.dgshengxigongchengsl.cn/1773851428741.html)
- [['https://www.fyluanpu.cn/325644463170.html', '漳州试管婴儿哪家好？本地口碑最好的助孕公司排名前三']](https://www.fyluanpu.cn/325644463170.html)
- [['https://www.sdhuabenhuanbao.cn/zhenshijingli/78.html', '上海优孕助孕靠不靠谱？教你通过三个细节分辨机构真假']](https://www.sdhuabenhuanbao.cn/zhenshijingli/78.html)
- [['https://www.sdxxy.cn/20250606-494.html', '如何提高试管代生咨询成功率！最有效的方法还是得从自身入手']](https://www.sdxxy.cn/20250606-494.html)
- [['https://www.dyqlsu.com/20251014-407.html', '代孕宝宝特征：昆明三代试管婴儿可以做吗？']](https://www.dyqlsu.com/20251014-407.html)
- [['https://www.sdwmtgccl.cn/22444569459145.html', '2026年国内做借卵试管生男孩价格多少？附费明细？ ,最著名代孕试管医院']](https://www.sdwmtgccl.cn/22444569459145.html)
- [['https://www.bkudgf.cn/174.html', '国内三代试管婴儿医院推荐及费用解析']](https://www.bkudgf.cn/174.html)
- [['https://www.mymydz.cn/207631680258.html', '助孕供卵代生-40岁以上女性做试管婴儿容易畸形吗？要十万吗？']](https://www.mymydz.cn/207631680258.html)
- [['https://www.dhsuzouzy.cn/14679186129669.html', '广州大学城中医院双胞胎试管助孕费用解析：详细了解代生和供卵代怀']](https://www.dhsuzouzy.cn/14679186129669.html)
- [['https://www.chengdusokh.cn/312280515295.html', '汕头供卵产子机构电话：24小时汕头助孕在线咨询']](https://www.chengdusokh.cn/312280515295.html)
- [['https://www.qumengru.com/224073805356.html', '上海代生就找开心帼,上海中山医院生殖科周六日有门诊吗？上午几点开门？']](https://www.qumengru.com/224073805356.html)
- [['https://www.vhpowpj.cn/20250821-127.html', '宝妈亲述：广州三代试管助孕全攻略，详细流程解析']](https://www.vhpowpj.cn/20250821-127.html)
- [['https://www.dyokx.com/shiguandaihuaijiage/103.html', '试管供卵助-杭州无精症中医！杭州市富阳区中医院调理无精症比较好的医生有哪些']](https://www.dyokx.com/shiguandaihuaijiage/103.html)
- [['https://www.sgdaiyun.com/205884349046.html', '杭州做试管龙凤胎的费用-杭州代生价格表,在杭州做了输卵管切除术后会不会宫外孕！切除输卵管会宫外孕吗！']](https://www.sgdaiyun.com/205884349046.html)
- [['https://www.hbhuihaohb.cn/154.html', '三代试管婴儿对卵巢与子宫条件的具体要求解析']](https://www.hbhuihaohb.cn/154.html)
- [['https://www.huaiyunq.cn/221864770069.html', '海总医院供卵试管婴儿成功率如何？姐妹真实分享与经验解析']](https://www.huaiyunq.cn/221864770069.html)
- [['https://www.gyzhixiao.cn/451.html', '借卵需要流程：O型血做试管婴儿可以避免abo溶血吗？']](https://www.gyzhixiao.cn/451.html)
- [['https://www.eduency.com/120661905210.html', '代生公司价格']](https://www.eduency.com/120661905210.html)
- [['https://www.xmxinyhwzhs.cn/19337996619692.html', '试管人工周期怀孕孕酮低，试管人工周期孕酮低有没有关系？']](https://www.xmxinyhwzhs.cn/19337996619692.html)
- [['https://www.sjzgwfjwzhs.cn/21895367944840.html', '深圳哪些代生女孩价格医院成功率高本地人推荐这些医院？附代生女孩价格医院推荐！']](https://www.sjzgwfjwzhs.cn/21895367944840.html)
- [['https://www.tjsjyongsheng.cn/210254439405.html', '2026福州供卵的私立机构汇总，附供卵三代生男孩详细步骤']](https://www.tjsjyongsheng.cn/210254439405.html)
- [['https://www.qzmx56.com/534.html', '取卵后哪些人容易腹水，什么体质的人取卵容易有腹水？']](https://www.qzmx56.com/534.html)
- [['https://www.jmxmintuhg.cn/20250420-148.html', '第三代试管婴儿全程需多少天，期间需要做好哪些准备？']](https://www.jmxmintuhg.cn/20250420-148.html)
- [['https://www.sdjiaxin.net/544.html', '解冻复苏优质胚胎做三代移植多久会妊娠？']](https://www.sdjiaxin.net/544.html)
- [['https://www.syldezdhkj.cn/31746199740937.html', '从贵阳去马来西亚做试管婴儿代生价位(马来西亚做试管婴儿代生成功率和费用)']](https://www.syldezdhkj.cn/31746199740937.html)
- [['https://www.skiguo.cn/20250927-67.html', '吃榴莲有助于卵泡发育']](https://www.skiguo.cn/20250927-67.html)
- [['https://hangzhou.ccxwlkx.cn/382.html', '备孕期间胃不舒服怎么办（孕妇胃痛吃什么好）']](https://hangzhou.ccxwlkx.cn/382.html)
- [['https://www.uueamru.cn/20250821-32.html', '代孕生殖机构：试管促排卵的时间与成功率解析']](https://www.uueamru.cn/20250821-32.html)
- [['https://www.cmanrxrr.cn/3193726921108.html', '2026做借卵三代试管选性别成功率多少？附详细介绍？,试管供卵代怀生子机构']](https://www.cmanrxrr.cn/3193726921108.html)
- [['https://www.luruihang.com/2347.html', '遵义第三代试管婴儿哪些医院遵义医学院第三代试管婴儿']](https://www.luruihang.com/2347.html)
- [['https://www.chengyanghg.cn/333.html', '2026年郑州供卵试管婴儿费用及流程详解']](https://www.chengyanghg.cn/333.html)
- [['https://www.wqxmm.cn/310891522574.html', '代孕机构哪个靠谱,加拿大生孩子全过程 亲生体验分享']](https://www.wqxmm.cn/310891522574.html)
- [['https://www.dygsdyw.com/220093497103.html', '苏州助孕网助孕服务,苏州做试管婴儿哪几家医院是正规的？']](https://www.dygsdyw.com/220093497103.html)
- [['https://www.satghenga.cn/117320140575.html', '武汉锦欣医院骗局是真的吗？武汉锦欣医院评价汇总']](https://www.satghenga.cn/117320140575.html)
- [['https://www.afa2019.com/223772965577.html', '石家庄公立供卵-石家庄借卵试管价格,石家庄做试管促排选哪个生殖科医生好？梁莹怎么样？']](https://www.afa2019.com/223772965577.html)
- [['https://www.sjb493.cn/12531467517263.html', '代怀方法有哪些_代生正规的机构,婚礼为什么忌讳大腹部，和新娘犯冲只是其一']](https://www.sjb493.cn/12531467517263.html)
- [['https://www.gzgudadl.cn/2776310424794.html', '抽烟对做试管有影响吗,专业私人试管代孕机构,做第三代试管代孕']](https://www.gzgudadl.cn/2776310424794.html)
- [['https://www.ppmaas.com/baoshengnanhaishiguan/486.html', '看一下!福州多囊做代生子机构包男孩成功率?']](https://www.ppmaas.com/baoshengnanhaishiguan/486.html)
- [['https://www.monpun.com/8214055360771.html', '取卵术后护理与卵巢保护：应对腹水与不适']](https://www.monpun.com/8214055360771.html)
- [['https://www.bjwdzxkj.cn/1756910659051.html', '云南省试管助孕生宝宝费用总共多少,代孕试管包成功']](https://www.bjwdzxkj.cn/1756910659051.html)
- [['https://www.apkbwvg.cn/shiguanyingergonglve/178.html', '辅助生殖费用全面解析：试管婴儿各项费用明细']](https://www.apkbwvg.cn/shiguanyingergonglve/178.html)
- [['https://www.hflrwzhs.cn/129.html', '同性伴侣的三代试管技术：如何实现血缘纽带的最大化？']](https://www.hflrwzhs.cn/129.html)
- [['https://www.fmngst.com/2442252291958.html', 'amh0.06，卵巢早衰，33岁拼一胎,爱心代孕公司']](https://www.fmngst.com/2442252291958.html)
- [['https://www.mimi567.com/219.html', '能做借卵试管:做试管婴儿前水果能吃吗（试管婴儿移植后吃什么水果好）']](https://www.mimi567.com/219.html)
- [['https://www.ewdboe.cn/318562263550.html', '湖北三代借卵试管助孕费用解析：辅助生殖技术与性别选择的合规考量']](https://www.ewdboe.cn/318562263550.html)
- [['https://www.cndcxc.com/daiyunliucheng/20251023/17108.html', '试管供卵自怀，怀孕多久能用测孕棒测出来']](https://www.cndcxc.com/daiyunliucheng/20251023/17108.html)
- [['https://www.zrbbavaq.cn/16798478751056.html', '拉萨做试管的医院预约流程，第一名技术全面成功率高,试管代孕案例']](https://www.zrbbavaq.cn/16798478751056.html)
- [['https://www.cd-hssf.com/121640963223.html', '试管捐卵代怀：卵巢功能早衰会恢复吗，女生抽烟会影响卵泡发育吗']](https://www.cd-hssf.com/121640963223.html)

*本文整理自母婴健康资讯，仅供科普参考。*
