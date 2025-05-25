# lab-hw-2
# 算法Lab2 复杂 DNA 序列的比对实验报告


## 实验要求

编写代码，根据输入的query和reference 计算两者的匹配关系，输出的一组非重叠切割方式 {(q_st₁, q_en₁, r_st₁, r_en₁), (q_st₂, q_en₂, r_st₂, r_en₂), ..., (q_stₙ, q_enₙ, r_stₙ, r_enₙ)}，使得：每个切割片段query[q_stᵢ: q_enᵢ]与对应reference[r_stᵢ: r_enᵢ]的编辑距离最小。

编程语言不限。需使用图相关算法。

**正确性要求**：需要将算法的输出提交到评分网站，分数不得低于基线分数。

**复杂度要求**：时间复杂度不得高于平方量级O(mn)（不含平方量级）。其中，复杂度几 乎为线性量级可得满分，平方量级可得大部分分数

## 实验思路

首先将参考序列 `reference` 表示为图中的路径，query 序列的每个字符在图上尝试与 reference 匹配，图的每个节点表示 query 中的一个位置与 reference 中的一个位置之间的匹配状态。从 `query` 和 `reference` 的每个可能的起始对 `(i, j)` 开始，逐个尝试延伸匹配路径。如果 `query[i] == reference[j]`，则记录匹配并向后扩展，直到无法继续。在得到最长匹配路径后对其中的gap进行填充，最终返回连续的最长匹配路径。

## 算法伪代码

```pseudocode
Input:
    reference: 字符串，长度 m
    query: 字符串，长度 n
    k: k-mer 长度，人为调整
    ref_window_size: gap 比对时的参考窗口大小，人为调整
Output:
    matched_intervals: List of (start_in_query, end_in_query, start_in_ref, end_in_ref)

Procedure align_sequence(reference, query, k, ref_window_size):
    // 构建 reference 的 k-mer 索引
       ref_kmers = empty dict
       for i = 0 to m - k:
           kmer = reference[i:i+k]
           ref_kmers[kmer].append(i)

    // 找到 query 中匹配的 anchors（forward 和 reverse-complement 两个方向）
       anchors = []
       for strand in [query, reverse_complement(query)]:
           for i = 0 to n - k:
               kmer = strand[i:i+k]
               if kmer in ref_kmers:
                   for pos in ref_kmers[kmer]:
                       (q_start, q_end), (r_start, r_end) = extend_match(strand, reference, i, pos)
                       anchors.append(((q_start, q_end), (r_start, r_end), strand_direction))

    // 对 anchors 进行最长链式连接（DAG 上的最长路径）
       dp[i] = anchors[i].length
       prev[i] = None
       for i = 0 to len(anchors) - 1:
           for j = 0 to i - 1:
               if anchors[j] can chain to anchors[i]:
                   score = dp[j] + overlap_score(anchors[j], anchors[i])
                   if score > dp[i]:
                       dp[i] = score
                       prev[i] = j

       使用 prev[] 回溯得到 optimal anchor path

    // 对于每对相邻 anchors，检测 gap 区域：
       if gap exists between q1_end and q2_start:
           提取 query gap 区段
           在 reference 中滑动 window：
               for each ref_window:
                   计算 edit distance 与 gap 匹配得分
           选最优 ref window，对应为 alignment 补全段

    // 将 anchor 对齐段和 gap 对齐段合并，生成最终 matched_intervals
    matched_intervals = merge (optimal anchor path,alignment 补全段)

    Return matched_intervals

```



## 时空复杂度

本算法的总体复杂度为那么 **O(n⋅m)**，其中， `reference` 长度为 m，`query` 长度为 n 。

以下是主要函数的复杂度评估：

### find_kmers(seq, k)

```python
for i in range(len(seq) - k + 1):
    kmer = seq[i:i+k]
```

- **调用一次，复杂度：** O(m)

  遍历 reference 的长度 m，每次切片长度为 k，即 O(m⋅k)

### 匹配 anchor 的两轮循环（forward 和 reverse complement）

```python
for i in range(len(query) - k + 1):  # 最多 n 次
    kmer = query[i:i+k]
    if kmer in ref_kmers:  # 哈希查找 O(1)
        for j in ref_kmers[kmer]:  # 假设平均命中次数为 α
            # 向两边延伸匹配，最多 O(min(n, m))
```

- 最坏情况：每个 query kmer 命中多个 reference 位置，每个进行线性扩展：

- **复杂度：** O(n⋅α⋅L)，其中 L 是最大延伸长度，最坏是 O(n)

- 实际上这个部分是代码中最耗时的地方，**最坏情况下为：**

  O(n⋅m)如果 ref_kmers 命中率高或无重复限制。

### chain_anchors

这是对找到的 anchors 进行 DAG 最长路径构造：

```python
for i in range(n):
    for j in range(i):
        score_alignment(...)  # O(1)
```

- 假设找到的 anchor 数为 a，那么复杂度为：
- **复杂度：** O(a^2)

### create_contiguous_alignments

该函数包含了**gap 部分的局部比对（暴力 levenshtein）**：

```python
for j in range(ref_window_start, ref_window_end, step):  # 常为 50 步
    levenshtein_distance(gap_query, r_segment)  # O(g^2)
```

- 假设有 g 个 gap，每个 gap 长度为 l，window 搜索范围 50 步；
- **复杂度：** O(g⋅l^2)

### 总体时间复杂度总结

设：

- n= query 长度
- m = reference 长度
- a= 找到的 anchor 数
- l=平均 gap 区段长度
- g= gap 区段数

综合各部分：

${ O(m \cdot k) + O(n) + O(n \cdot m) + O(a^2) + O(g \cdot l^2) }$

在实际使用中，如果：

- k 足够大、重复少，则 $a \ll n \cdot m$，可以忽略 $a^2$
- gap 数不多，$g \cdot l^2 \ll n \cdot m$

那么 **主导复杂度为 $O(n \cdot m)$**

