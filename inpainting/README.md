# Protein Inpainting

Codes for running Protein Inpainting developed by Joe Watson and David Juergens 

## Description 

Protein inpainting is a method for "conditional joint protein
sequence/structure generation". This means that given some combination of
protein sequence and structure that you have, you can use this method to
simultaneously generate more sequence and structure conditioned on that input. 

**Things inpainting is good at:**

- Refining non-ideal parts of proteins 
- Resampling protein structures near a starting structure 
- Re-looping proteins (i.e., keep tertiary/secondary structure but changing the
  order in which elements appear in sequence space)
- Rigidly fusing two protein domains 
- Loop building 
- Scaffolding medium-sized motifs

**Things that are currently challenging with inpainting:**
- Generating large amounts of protein from very little or no protein structure
  (it's worth a try, but don't expect whole proteins to come out consistently) 
- Inpainting with excluded volumes (though it can be done)

## Usage 

In this section we will demonstrate how to run inpainting. 

The actual script you will execute is called `inpaint.py`. There are many ways
to run it, which we will explore through examples below. 

### Basic execution

A baseline execution of the script requires 3 pieces of
information, which are provided in the form of flags coming in from the command
line that you provide. The three things are 

1. A template protein structure/sequence, in the form of a pdb file. (`--pdb`)

2. Specification of which parts of the protein are being kept, removed, and
inpainted. This is provided as a "contigs string". (`--contigs`) 

3. A location where output from the script should be written. (`--out`)

```
python inpaint.py --pdb 2KL8.pdb --contigs A25-50,10,A61-79 --out pdbs_test/auto_out 
```

The first flag denotes that the pdb file `2KL8.pdb` will be used as the source for any template protein information we want to scaffold. 

The second flag is the contigs string, which says a few things:
- The A25-50 denotes the first contiguously *kept* region ("contig") - meaning
  that we are essentially going to copy the sequence and structure of all
  residues from A25 to A50 in `2KL8.pdb` (inclusive). Since `A25` comes first
  in the contigs string, it will be the first residue in the output protein.
  Note that even though `2KL8.pdb` also contains residues A1-24, they will not
  be considered during design because they are not in the contigs string.  
- The `,10` denotes that *directly* attached to `A1-25` there will be 10
  inpainted residues, where both sequence and structure will be generated. 
- The `,A61-79` denotes that *directly* attached to the 10-residue inpainted
  region we will have residues `A61-79` from `2KL8.pdb`.  

Finally, the last flag denotes where output will go. This is always in the form
`path/to/my/outdir/prefix` - that is, it begins with the path to a folder, and
the last bit denotes the prefix that all files will be tagged with when they go
into that folder. The standard files which are output is a `.pdb` file for each
design, allong with a metadata `.trb` file which contains information about how
your contigs ("kept regions") were mapped from the input protein to the output
protein, the RoseTTAFold pLDDT predictions, etc. 

**NOTE:** within a given design run with `inpaint.py`, any files being output
into your output direectory will not clobber each other because each individual
design will be tagged with a unique ID number. But designs from separate runs
can clobber each other. 

## Understanding the --contigs flag in depth.

1. Variable length inpainting. In many cases, you don't know the exact number
of residues needed to fill a gap. Therefore, you can specify a range, e.g.
`A25-50,7-13,A61-79`. This will allow between 7 and 13 (inclusive) residues to
be inpainted. However, there is no 'intelligent' sampling here - the length is
currently just randomly sampled from the given range. Therefore, you'll need to
run multiple inpainting runs to sample this range, with e.g. `--num_designs 3`

2. Inpainting in the presence of another, fixed chain. Sometimes, you'll want
to inpaint in the presence of another chain, e.g. a receptor that you want to
inpaint a target against. This can be specified by adding the receptor PDB indices
as a separate block in the contig, separated by a ',0' indicating a chain break.
E.g. if you wanted to inpaint in the presence of an 100-amino acid B chain (1-indexed),
the contig would be `B1-100,0 A25-50,7-13,A61-79`. This will put a 200-residue
index break into the input, which will allow the inpainting network to see this
B chain as a separate chain.

3. Fusing multiple chains. Your 'visible' parts of structure can be linked
together with inpainting, even if they're on separate chains in the input file.
E.g. `A25-50,15-25,B1-50` will fuse part of the 'A' chain to part of the 'B'
chain.

## Understanding the outputs
> 任何有支架的基序最终都会出现在链A上。如果你包含一个受体链，它将在链B上输出。
1. The PDB file. Any scaffolded motifs will end up on chain A. If you included
a receptor chain, that will be output on chain B.
> numpy 查看trb文件，包含lddt, inpaint部分的lddt,
2. The trb file. This contains metadata about the inpainting run. Open this
with np.load([File], allow_pickle=True).
    > 这是修复网络对输出结构‘好坏’的预测。这是按氨基酸残基计算的，并包括所有未被遮罩的区域
    - `lddt`- This is the inpainting network's prediction at how 'good' the
      output structure is. This is per residue, and includes all of the
      unmasked region.
    > 这只是 'lddt' 输出中对应于被修复区域的部分。通常，我们会根据该输出的平均值进行筛选，以获取‘最优’的大约 5-10% 的结果。
    > lddt：指的是 Local Distance Difference Test，是一种用于衡量蛋白质结构预测准确性的评分指标，值通常在 0 到 1（或 0 到 100）之间。数值越高，表示预测越接近真实结构。
    > take the 'best' ~5-10% of outputs（取最优的约5-10%结果）：在所有预测结果中，只保留 lddt 平均分数最高的大约前 5% 到 10%，作为最可靠的预测结构。
    - `inpaint_lddt` - This is just the part of the 'lddt' output corresponding
      to the inpainted region. Normally, we filter on the mean of this output,
      to take the 'best' ~5-10% of outputs.
    - details about mapping (i.e. how residues in the input map to residues in
      the output)
      > con_ref_pdb_idx/con_hal_pdb_idx - 这是两个数组，分别包含输入PDB中的索引（在con_ref_pdb_idx中）和它们在输出PDB中的位置（在con_hal_pdb_idx中）。
      > 这仅包含进行了修复（inpainting）的链（即不包含任何固定的受体/靶标链）。
        - `con_ref_pdb_idx`/`con_hal_pdb_idx` - These are two arrays including
          the input pdb indices (in con_ref_pdb_idx), and where they are in the
          output pdb (in con_hal_pdb_idx). This only contains the chains where
          inpainting took place (i.e. not any fixed receptor/target chains)
      > con_ref_idx0/con_hal_idx0 - 这些与上述参数相同，但采用从0开始的索引（0-indexed），且不包含链信息。这对提取坐标（用于评估比对等）很有用。
        - `con_ref_idx0`/`con_hal_idx0` - These are the same as above, but 0
          indexed, and without chain information. This is useful for splicing
          coordinates out (to assess alignment etc).
      > complex_con...和receptor_con... - 如果您包含了固定的受体/靶标链，这些链会单独出现在receptor_con...中，或者与修复后的链一起出现在complex_con...中。
        - `complex_con...` and `receptor_con...` - If you have included fixed
          receptor/target chains, these will be included either on their own
          (in receptor_con...) or with the inpainted chains (in
          complex_con...).
      > sampled_mask - 如果您指定了要修复的长度范围，这个sampled_mask将给出该次特定修复运行中实际使用的精确长度。
        - `sampled_mask` - if you've specified a range of lengths to be
          inpainted, this `sampled_mask` will give the precise length used
          during that specific inpainting run
> .npz文件 - 该文件包含来自RosettaFold预测的成对距离和取向角度。它将作为trfold_relax.sh脚本的输入，用于添加侧链并进一步优化结构。不过，如果您只是计划对设计进行AlphaFold2预测（我们经常这样做），则不需要这个npz文件。
3. The .npz file. This contains the predicted pairwise distances and
orientation angles from RosettaFold. This is used as an input to the
`trfold_relax.sh` script to add sidechains and refine the structure further.
However, if you simply plan to make an AlphaFold2 prediction of the design (as
we often do), the npz isn't needed.

## Other flags

- `--inpaint_seq`: This allows the sequence of residues to be masked, while the
  backbone coordinates are given. This is handy if you're, for example, making
  a fusion between two monomeric proteins, and part of what was the surface is
  now going to be in the core of the protein. This is specified similarly to
  the contig string, but without chain breaks (e.g. `A1-3,A4,A6,B10-100`)
- `--inpaint_str`: This allows the structure of a residue to be masked, but not
  its sequence. This essentially asks the network to predict the structure of
  this residue.
- `--res_translate`: Which residues to translate (randomly in x, y and z
  direction), with maximum distance to translate specified, e.g. `A35,2:B22,4`
  translates residue A35 up to 2A in a random direction, and B22 up to 4A. If
  specified residues are in masked --window, they will be unmasked. In --contig
  mode, residues must not be masked (as need to know where to put them. Default
  distance to translate is 2A.
- `--tie_translate`: For randomly translating multiple residues together (e.g.
  to move a whole secondary structure element). Syntax is e.g.
  `A22,A27,A30,4.0:A48,A50` which would randomly move residues A22, A27 and A30
  together up to 4A, and A48 and A50 together (but in a different random
  direction/distance to the first block) to a default distance of up to 2A.
  Alternatively, residues can be specifed like `A12-26,6.0:A40-52,A56`. This
  can be specified alongside --res_translate, so some residues are tied, and
  some are not, but if residues are specified in both, they will only be moved
  in their tied block (i.e. their --res_translate will be ignored)
- `--block_rotate`: Do you want to rotate a whole structural block (or single
  residue)? Syntax is same as tie_translate. Rotation is in degrees.
- `--num_designs`: Number of designs (inpaints) to generate. Note INPAINTING IS
  DETERMINISTIC. Therefore, you need at least as many possible input
  combinations (e.g. lengths of inpainted regions etc) to generate actual
  diversity. Otherwise, identical inputs in = identical inputs out (and a LOT
  of wasted compute)

All other flags can generally be left as default, but either dive into the code
or ask us if you have any questions.

## RF Joint for mutation effect prediction 

usage: python rfjoint_mutation_effect_prediction.py -msa input_msa.a3m -out_dir output_mut -out_file_name mut_effect_pred -templ input_pdb.pdb (optional) --checkpoint BFF_mix_epoch25.pt

For evaluation: 
    output numpy file is of shape Lx21 (20 amino acids + gap). The effect of a mutation at position i is calculated as: 
    
        alphabet = 'ARNDCQEGHILKMFPSTWYV-'
        
        wt_aa = alphabet.index('A') # wild type amino acid at position i
        
        mt_aa = alphabet.index('P') # mutated amino acid at position i
        # A突变为P的预测效果（如果是阳性，预测是有益的，如果是阴性，预测是有害的）
        y_hat = pred[i][mt_aa] - pred[i][wt_aa]  # predicted effect of mutating A to P (if positive, predicted to be beneficial, if negative, predicted to be deleterios) 
        

## FAQs

1. 'How much protein can inpainting inpaint?' 
> 这取决于具体问题，但一般而言，对于超过约60个残基的修复（inpainting），模型会较为困难。如果你有多个需要修复的区域（位于可见结构块之间），那么待修复的蛋白质总量可能会远超这个数字（例如，每段约50个残基时）。
This depends on the problem, but generally it will struggle with inpainting
more than around 60 residues. If you have multiple regions to be inpainted
(between visible blocks), the total amount of protein to be inpainted could be
quite a lot more than this though (if each segment was say, around 50
residues).

2. 'I see chain breaks/clashes in my pdb output'
> 修复（inpainting）过程经常失败，这可能是因为在该次运行中提供给网络的长度参数组合无法生成合理的蛋白质结构（至少在修复模式下如此）。不过，这些设计通常会有较低的“inpaint_lddt”平均分，因此如果只保留得分最高的前10%左右的设计，这些结果应该不会存在严重的原子冲突等问题。
Inpainting often fails, probably because the set of lengths given to the
network during that run were incompatible with making a good protein (at least
with inpainting). These designs will generally have a low mean 'inpaint_lddt'
metric though, so if you only take the top-scoring 10% or so of designs, these
shouldn't have bad clashes etc.

3. 'What is a good cutoff for the mean 'inpaint_lddt' score?
> 这个参数的原始值会受到多种因素影响，因此会因具体问题而异。通常，我会直接查看一批具有不同"inpaint_lddt"评分的输出结果，据此确定一个评分阈值。或者，我就直接选取得分最高的前10%的输出结果。
The raw value for this depends on a range of factors, so will vary
problem-by-problem. Normally, I just visually inspect a few outputs with a
range of 'inpaint_lddt' metrics, and choose a cutoff based on this. Or, I just
take the top-scoring 10% of outputs.

## Authors and acknowledgment

This work was developed by Joseph Watson (jwatson3@uw.edu), David Juergens
(davidcj@uw.edu), Jue Wang (jue@uw.edu) and Woody Ahern (ahern@uw.edu)

