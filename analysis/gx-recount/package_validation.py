#!/usr/bin/env python3
import hashlib,json,pathlib,shutil,tarfile,sys
root=pathlib.Path(sys.argv[1]);work=root/'work/gx-recount';code=root/'work/crouch-cleanup-meta/analysis/gx-recount';out=work/'R7CE_full_GX';dest=root/'outputs/gx-recount';stage=dest/'gx_validation';dest.mkdir(exist_ok=True);stage.mkdir(exist_ok=True)
prov=json.loads((out/'provenance.json').read_text());idx=json.loads((work/'BAM_INDEX_METADATA.json').read_text());src=json.loads((work/'BAM_VERIFIED.json').read_text())
assert prov['full_BAM_processed'] is True and prov['counters']['records']==idx['expected_total_alignment_records']
assert prov['author_cohort_cells']==1759 and prov['all_reference_genes']==33694 and not prov['normalized']
assert src['ENA_MD5_verified'] and src['BAM_quickcheck_pass']
files=[(out/'GW18_R7CE_BAM_GX_author_cells_raw_counts.rds','GW18_R7CE_BAM_GX_author_cells_raw_counts.rds'),(work/'inputs/author_50_R7CE_slots.rds','author_50_R7CE_slots.rds'),(code/'compare_GX_STARsolo_RStudio.R','compare_GX_STARsolo_RStudio.R')]
for name in ['final_audit.json','report.txt','provenance.json','source_provenance.json','GX_vs_author_50_cells.tsv','GX_vs_author_summary.tsv','author_counts_ratio_summary.tsv','author_counts_ratio_notes.txt','author_counts_over_GX_per_cell_ratios.tsv','author_counts_over_GX_per_gene_ratios.tsv']:
 files.append((out/name,name))
files.extend([(work/'BAM_VERIFIED.json','BAM_VERIFIED.json'),(work/'BAM_INDEX_METADATA.json','BAM_INDEX_METADATA.json')])
for source,name in files:
 assert source.is_file(),source;shutil.copy2(source,stage/name)
readme='''GW18 R7CE 完整 BAM GX 计数验证包

本包仅对应 GW18 的 R7ctxendo 文库，保留作者的全部 1,759 个细胞、旧参考全部 33,694 个基因。不是全部 GW18，也不是试点或部分染色体结果。
原始 BAM 已通过 ENA MD5、BAM 格式及索引 MD5 校验；扫描记录数与 BAI 独立统计一致。

Windows 操作：
1. 在 Git Bash 解压桌面上的包：
   tar -xzf "$HOME/Desktop/crouch_GW18_R7CE_GX_validation.tar.gz" -C "$HOME/Desktop"
2. 在 RStudio 打开 gx_validation/compare_GX_STARsolo_RStudio.R，点击 Source（如需指定编码请选择 UTF-8）。只需要 Matrix 包。
3. 依次选择：
   本包 GW18_R7CE_BAM_GX_author_cells_raw_counts.rds
   之前通过最终核验的 GW18_author_cells_GRCh38_1.2.0_raw_counts.rds（13,906 细胞、33,694 基因、28,758,442 UMIs）
   本包 author_50_R7CE_slots.rds
4. 输出保存在 GX RDS 旁新建的 GX_STARsolo_validation_日期时间/ 目录。
5. 将 GX_vs_STARsolo_summary.tsv 与 GX_STARsolo_vs_author_summary.tsv 的内容或截图发回来即可，不需要回传大矩阵。

计数：主比对、MAPQ255、唯一 GX、完整有效 CB/UB，按同一物理文库内（完整 CB 含 GEM 后缀，GX，UB）去重。没有额外细胞 QC、基因删除、归一化、log 或 scale。
PCR duplicate 标记本身不用于删去 reads；同一 UMI 的多条支持 reads 最终计一次。跨基因 CB/UB 冲突单独统计，没有假设等同于作者 Cell Ranger 全部原始分子过滤规则。
作者 RNA/counts 是分数表达值，RNA/data 是处理后表达，均不作为原始整数 UMI 真值。先看 GX 与 STARsolo 的精确计数差异，再看相关性。
作者表达比较固定同一 50 个细胞、唯一共同 gene symbols、相同三方表达基因并集。归一化分母使用各自全部原始参考基因。
'''
(stage/'README_先读我.txt').write_text(readme)
manifest=[]
for f in sorted(stage.iterdir()):
 if f.is_file() and f.name!='SHA256SUMS':manifest.append(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.name)
(stage/'SHA256SUMS').write_text('\n'.join(manifest)+'\n')
archive=dest/'crouch_GW18_R7CE_GX_validation.tar.gz'
assert not archive.exists(),'Preserving existing archive'
with tarfile.open(archive,'w:gz') as t:t.add(stage,arcname='gx_validation')
digest=hashlib.sha256(archive.read_bytes()).hexdigest();(dest/(archive.name+'.SHA256SUMS')).write_text(digest+'  '+archive.name+'\n')
print('READY',archive,'bytes',archive.stat().st_size,'SHA256',digest)
