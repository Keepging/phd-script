# %%
import pandas as pd          
import json                  
import os                    
from pathlib import Path     
import numpy as np           
print("库导入成功")

# %%
# 输入文件路径
JSON_FOLDER = r"D:\博士\Phospho\biophys_json"           
SITES_TABLE = r"D:\博士\Phospho\Phospho_Site_Compartment_Table.csv"  

# 输出文件路径
OUTPUT_FILE = r"D:\博士\Phospho\Sites_with_Biophysical_Features.csv"

# 验证文件是否存在
print("🔍 验证文件路径...")
print(f"  JSON文件夹: {'✓' if os.path.exists(JSON_FOLDER) else '✗'} {JSON_FOLDER}")
print(f"  位点表格:   {'✓' if os.path.exists(SITES_TABLE) else '✗'} {SITES_TABLE}")

# 统计JSON文件数量
json_files = [f for f in os.listdir(JSON_FOLDER) if f.endswith('.json')]
print(f"\n📊 JSON文件数量: {len(json_files):,}")

# %%
"""
从JSON文件建立基因名映射

"""

import json
import os
import re

print("🔄 从JSON文件建立基因名到Uniprot ID的映射...")

# 获取所有JSON文件
json_files = [f for f in os.listdir(JSON_FOLDER) if f.endswith('.json')]
print(f"JSON文件总数: {len(json_files):,}")

# 初始化映射字典
gene_to_uniprot = {}
uniprot_to_gene = {}
no_gene_name = []

# 遍历所有JSON文件，提取基因名
print("\n⏳ 扫描JSON文件...")
for i, json_file in enumerate(json_files):
    # 每2000个文件显示进度
    if (i+1) % 2000 == 0:
        print(f"  进度: {i+1:,}/{len(json_files):,} ({(i+1)/len(json_files)*100:.1f}%)")
    
    json_path = os.path.join(JSON_FOLDER, json_file)
    
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        uniprot_id = data.get('protein_id', '')
        header = data.get('uniprot_header', '')
        
        # 从header提取基因名 (格式: ...GN=GENE_NAME ...)
        match = re.search(r'GN=([^\s]+)', header)
        
        if match:
            gene_name = match.group(1)
            
            # 建立双向映射
            gene_to_uniprot[gene_name] = uniprot_id
            uniprot_to_gene[uniprot_id] = gene_name
        else:
            # 记录没有基因名的蛋白质
            no_gene_name.append(uniprot_id)
            
    except Exception as e:
        # 跳过读取失败的文件
        continue

print(f"\n✅ 映射建立完成!")
print(f"  成功提取: {len(gene_to_uniprot):,} 个基因")
print(f"  无基因名: {len(no_gene_name):,} 个")

# 显示映射示例
print(f"\n📋 映射示例 (前10个):")
# 获取前10个键值对，避免字典在迭代时大小变化的问题
example_genes = list(gene_to_uniprot.keys())[:10]
for gene in example_genes:
    print(f"  {gene:15s} → {gene_to_uniprot[gene]}")

# 如果需要，可以将无基因名的ID打印出来一部分看看
if len(no_gene_name) > 0:
    print(f"\n⚠️ 无基因名的Uniprot ID示例 (前5个):")
    for pid in no_gene_name[:5]:
        print(f"  {pid}")


# %%
# 读取位点表格
sites_df = pd.read_csv(SITES_TABLE)

print("📋 位点表格结构:")
print(f"  总行数: {len(sites_df):,}")
print(f"  总列数: {len(sites_df.columns)}")
print(f"\n列名:")
for i, col in enumerate(sites_df.columns, 1):
    print(f"  {i:2d}. {col}")

print("\n前3行数据:")
print(sites_df.head(3))

# 统计独特基因数
unique_genes = sites_df['Gene'].unique()
print(f"\n📊 数据统计:")
print(f"  总位点数: {len(sites_df):,}")
print(f"  涉及基因数: {len(unique_genes):,}")
print(f"  平均每个基因: {len(sites_df)/len(unique_genes):.1f} 个位点")


# %%
def extract_site_features(gene_name, position, gene_to_uniprot_dict, json_folder_path):
    """
    从JSON文件中提取指定磷酸化位点的生物物理特征
    
    参数 Parameters:
        gene_name (str): 基因名，如 'ABCF1'
        position (int): 氨基酸位置，如 109
        gene_to_uniprot_dict (dict): 基因名到Uniprot ID的映射字典
        json_folder_path (str): JSON文件夹路径
    
    返回 Returns:
        dict: 包含生物物理特征的字典，如果提取失败返回None
    """
    
    # 步骤1: 查找对应的Uniprot ID
    if gene_name not in gene_to_uniprot_dict:
        return {'error': 'gene_not_in_mapping'}
    
    uniprot_id = gene_to_uniprot_dict[gene_name]
    
    # 步骤2: 构建JSON文件路径
    json_filename = f"{uniprot_id}.json"
    json_filepath = os.path.join(json_folder_path, json_filename)
    
    # 步骤3: 检查JSON文件是否存在
    if not os.path.exists(json_filepath):
        return {'error': 'json_file_not_found'}
    
    # 步骤4: 读取JSON文件
    try:
        with open(json_filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        return {'error': f'json_read_error: {e}'}
    
    # 步骤5: 提取指定位置的氨基酸信息
    residues = data.get('residues', [])
    
    # 转换为0-based索引 (Python列表从0开始，但氨基酸位置从1开始)
    position_idx = int(position) - 1
    
    # 检查位置是否超出范围
    if position_idx < 0 or position_idx >= len(residues):
        return {'error': 'position_out_of_range'}
    
    # 步骤6: 获取该位置的所有特征
    residue_data = residues[position_idx]
    
    # 步骤7: 提取我们需要的生物物理特征
    features = {
        'backbone_dynamics': residue_data.get('backbone', None),
        'sidechain_dynamics': residue_data.get('sidechain', None),
        'disorder_propensity': residue_data.get('disoMine', None),
        'helix_propensity': residue_data.get('helix', None),
        'sheet_propensity': residue_data.get('sheet', None),
        'coil_propensity': residue_data.get('coil', None),
        'amino_acid': residue_data.get('aa', None),
        'sequence_position': residue_data.get('seqpos', None),
        'earlyFolding': residue_data.get('earlyFolding', None),
        'uniprot_id': uniprot_id,
        'error': None
    }
    
    return features


# 测试函数
print("🧪 测试特征提取函数...")
test_gene = sites_df.iloc[0]['Gene']
test_position = sites_df.iloc[0]['Position']

test_result = extract_site_features(
    test_gene, 
    test_position, 
    gene_to_uniprot, 
    JSON_FOLDER
)

print(f"\n测试位点: {test_gene}_{test_position}")
if test_result.get('error') is None:
    print("✅ 提取成功!")
    for key, value in test_result.items():
        if value is not None:
            print(f"  {key}: {value}")
else:
    print(f"❌ 提取失败: {test_result['error']}")

# %%
"""
批量为所有11,046个磷酸化位点提取生物物理特征

"""

print("开始批量特征提取...")
print(f"总位点数: {len(sites_df):,}\n")

# 初始化结果列表
all_features = []
error_counts = {
    'gene_not_in_mapping': 0,
    'json_file_not_found': 0,
    'position_out_of_range': 0,
    'json_read_error': 0,
    'success': 0
}

# 遍历每个位点
for idx, row in sites_df.iterrows():
    # 每1000个位点显示进度
    if (idx + 1) % 1000 == 0:
        print(f"  进度: {idx+1:,}/{len(sites_df):,} ({(idx+1)/len(sites_df)*100:.1f}%)")
    
    # 提取特征
    features = extract_site_features(
        gene_name=row['Gene'],
        position=row['Position'],
        gene_to_uniprot_dict=gene_to_uniprot,
        json_folder_path=JSON_FOLDER
    )
    
    # 添加位点标识
    features['PTM_collapse_key'] = row['PTM_collapse_key']
    features['Gene'] = row['Gene']
    features['Site'] = row['Site']
    features['Position'] = row['Position']
    
    # 记录错误统计
    if features['error'] is None:
        error_counts['success'] += 1
    else:
        error_type = features['error'].split(':')[0]
        if error_type in error_counts:
            error_counts[error_type] += 1
    
    all_features.append(features)

print("\n✅ 特征提取完成!")

# 显示统计结果
print("\n📊 提取结果统计:")
print(f"  成功: {error_counts['success']:,} ({error_counts['success']/len(sites_df)*100:.1f}%)")
print(f"  失败: {len(sites_df) - error_counts['success']:,}")
print("\n失败原因分布:")
for error_type, count in error_counts.items():
    if error_type != 'success' and count > 0:
        print(f"  - {error_type}: {count:,}")


# %%
"""
将提取结果转换为DataFrame，并与原始定位表合并
"""

# 转换为DataFrame
features_df = pd.DataFrame(all_features)

print("📋 提取特征表结构:")
print(f"  行数: {len(features_df):,}")
print(f"  列数: {len(features_df.columns)}")

# 查看成功提取的数据
successful_features = features_df[features_df['error'].isna()]
print(f"\n✅ 成功提取特征的位点: {len(successful_features):,}")

# 显示前3个成功案例
print("\n前3个成功案例:")
display_cols = ['PTM_collapse_key', 'backbone_dynamics', 
                'disorder_propensity', 'helix_propensity', 'amino_acid']
print(successful_features[display_cols].head(3).to_string(index=False))

# 选择要保留的生物物理特征列
feature_cols = [
    'PTM_collapse_key',
    'backbone_dynamics',
    'sidechain_dynamics', 
    'disorder_propensity',
    'helix_propensity',
    'sheet_propensity',
    'coil_propensity',
    'amino_acid',
    'uniprot_id',
    'earlyFolding',
    'error'
]

features_to_merge = features_df[feature_cols]

# 左连接：保留原始表的所有行
final_df = sites_df.merge(
    features_to_merge, 
    on='PTM_collapse_key', 
    how='left'
)

print(f"\n📋 最终表格结构:")
print(f"  行数: {len(final_df):,}")
print(f"  列数: {len(final_df.columns)}")

print(f"\n新增的生物物理特征列:")
for col in feature_cols:
    if col not in ['PTM_collapse_key', 'error']:
        print(f"  - {col}")

print("\n✅ 数据合并完成！")


# %%
"""
Cell 8: Data Quality Check
检查提取数据的质量和完整性
"""

print("="*80)
print("🔍 数据质量检查报告")
print("="*80)

# 1. 缺失值统计
print("\n1️⃣ 缺失值统计:")
feature_columns = ['backbone_dynamics', 'sidechain_dynamics', 'disorder_propensity', 
                   'helix_propensity', 'sheet_propensity', 'coil_propensity', 'earlyFolding']

for col in feature_columns:
    if col in final_df.columns:
        missing = final_df[col].isna().sum()
        present = len(final_df) - missing
        print(f"  {col:25s}: {present:5,} 有效 / {missing:4,} 缺失 ({present/len(final_df)*100:.1f}%)")
    else:
        print(f"  {col:25s}: ⚠️  列不存在")

# 2. 数值范围检查（应该在0-1之间）
print("\n2️⃣ 数值范围检查 (应在0-1之间):")
for col in feature_columns:
    if col not in final_df.columns:
        continue
    valid_data = final_df[col].dropna()
    if len(valid_data) > 0:
        min_val = valid_data.min()
        max_val = valid_data.max()
        mean_val = valid_data.mean()
        
        in_range = (min_val >= 0) and (max_val <= 1)
        status = "✓" if in_range else "✗"
        
        print(f"  {status} {col:25s}: min={min_val:.3f}, max={max_val:.3f}, mean={mean_val:.3f}")

# 3. 特征分布统计
print("\n3️⃣ 特征分布统计 (5点汇总):")
available_features = [col for col in feature_columns if col in final_df.columns]
if available_features:
    print(final_df[available_features].describe().round(3).to_string())

# 4. 二级结构一致性检查
print("\n4️⃣ 二级结构一致性检查:")
if all(col in final_df.columns for col in ['helix_propensity', 'sheet_propensity', 'coil_propensity']):
    structure_sum = (
        final_df['helix_propensity'].fillna(0) + 
        final_df['sheet_propensity'].fillna(0) + 
        final_df['coil_propensity'].fillna(0)
    )
    valid_structure = structure_sum[structure_sum > 0.5]
    
    if len(valid_structure) > 0:
        mean_sum = valid_structure.mean()
        std_sum = valid_structure.std()
        
        print(f"  Helix+Sheet+Coil 总和:")
        print(f"    平均值: {mean_sum:.4f} (期望≈1.0)")
        print(f"    标准差: {std_sum:.4f}")
        
        if abs(mean_sum - 1.0) < 0.05:
            print(f"    状态: ✓ 一致性良好")
        else:
            print(f"    状态: ⚠ 存在偏差")
else:
    print("  ⚠️  二级结构列不完整，跳过检查")

# 5. 按氨基酸类型统计
print("\n5️⃣ 磷酸化氨基酸类型分布:")
if 'amino_acid' in final_df.columns:
    valid_aa = final_df[final_df['error'].isna()]['amino_acid']
    if len(valid_aa) > 0:
        aa_counts = valid_aa.value_counts()
        print(f"  总计成功提取: {aa_counts.sum():,} 个位点")
        for aa, count in aa_counts.items():
            print(f"    {aa}: {count:5,} ({count/aa_counts.sum()*100:.1f}%)")
else:
    print("  ⚠️  amino_acid列不存在")

# 6. 提取失败的详细原因
print("\n6️⃣ 提取失败详情:")
if 'error' in final_df.columns:
    failed_df = final_df[final_df['error'].notna()]
    if len(failed_df) > 0:
        error_counts = failed_df['error'].value_counts()
        for error, count in error_counts.items():
            print(f"  {error}: {count:,} 个位点")
        
        if 'Gene' in final_df.columns:
            print(f"\n  失败位点影响的基因数: {failed_df['Gene'].nunique():,}")
    else:
        print("  ✓ 无失败位点")
else:
    print("  ⚠️  error列不存在")

# 7. 数据覆盖率总结
print("\n7️⃣ 整体数据覆盖率:")
if 'backbone_dynamics' in final_df.columns:
    total_sites = len(final_df)
    valid_sites = final_df['backbone_dynamics'].notna().sum()
    coverage = valid_sites / total_sites * 100
    print(f"  总位点数: {total_sites:,}")
    print(f"  成功提取: {valid_sites:,} ({coverage:.1f}%)")
    print(f"  提取失败: {total_sites - valid_sites:,} ({100-coverage:.1f}%)")

print("\n" + "="*80)
print("✅ 质量检查完成")
print("="*80)


# %%
"""
保存最终结果表格和提取日志 (修正版 - 去除ModNumber)
"""

print("开始保存结果文件...\n")

# 1. 保存完整表格
final_df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8-sig')
file_size_kb = os.path.getsize(OUTPUT_FILE) / 1024
print(f"✅ 完整表格已保存:")
print(f"   文件: {OUTPUT_FILE}")
print(f"   大小: {file_size_kb:.1f} KB")
print(f"   行数: {len(final_df):,}")
print(f"   列数: {len(final_df.columns)}")

# 2. 保存仅成功提取的位点（用于后续分析）
success_df = final_df[final_df['error'].isna()].copy()
success_file = OUTPUT_FILE.replace('.csv', '_SuccessOnly.csv')
success_df.to_csv(success_file, index=False, encoding='utf-8-sig')
print(f"\n✅ 成功位点表已保存:")
print(f"   文件: {success_file}")
print(f"   位点数: {len(success_df):,}")

# 3. 保存失败记录（用于排查问题）
# ✅ 修正: 去除ModNumber列(已不存在)
failed_df = final_df[final_df['error'].notna()][
    ['PTM_collapse_key', 'Gene', 'Site', 'Position', 'error']
].copy()

if len(failed_df) > 0:
    failed_file = OUTPUT_FILE.replace('.csv', '_Failed.csv')
    failed_df.to_csv(failed_file, index=False, encoding='utf-8-sig')
    print(f"\n⚠️  失败记录已保存:")
    print(f"   文件: {failed_file}")
    print(f"   位点数: {len(failed_df):,}")

# 4. 生成提取报告
from datetime import datetime

report = f"""
================================================================================
磷酸化位点生物物理特征提取报告
Biophysical Feature Extraction Report
================================================================================

提取时间 Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

输入文件 Input Files:
  - 位点表格: {SITES_TABLE}
    总位点数: {len(sites_df):,}
    独特基因: {sites_df['Gene'].nunique():,}
  
  - JSON文件夹: {JSON_FOLDER}
    JSON文件数: {len([f for f in os.listdir(JSON_FOLDER) if f.endswith('.json')]):,}
    
  - 基因映射: 从JSON header提取 (GN= field)
    成功映射: {len(gene_to_uniprot):,} 个基因

输出文件 Output Files:
  - 完整表格: {OUTPUT_FILE}
  - 成功位点: {success_file}
  - 失败记录: {failed_file if len(failed_df) > 0 else 'N/A'}

================================================================================
提取结果统计 Extraction Statistics
================================================================================

总位点数 Total Sites:        {len(final_df):,}
成功提取 Success:            {len(success_df):,} ({len(success_df)/len(final_df)*100:.1f}%)
提取失败 Failed:             {len(failed_df):,} ({len(failed_df)/len(final_df)*100:.1f}%)

按氨基酸类型 By Amino Acid:
"""

aa_counts = success_df['amino_acid'].value_counts()
for aa, count in aa_counts.items():
    report += f"  {aa}: {count:,} ({count/len(success_df)*100:.1f}%)\n"

report += f"""
失败原因 Failure Reasons:
"""
if len(failed_df) > 0:
    error_counts = failed_df['error'].value_counts()
    for error, count in error_counts.items():
        report += f"  - {error}: {count:,}\n"
else:
    report += "  None\n"

report += f"""
================================================================================
提取的生物物理特征 Extracted Biophysical Features
================================================================================

1. Backbone Dynamics (骨架动力学)
   - 描述: 蛋白质主链的灵活性
   - 范围: 0-1 (0=刚性, 1=高度灵活)
   
2. Sidechain Dynamics (侧链动力学)
   - 描述: 氨基酸侧链的灵活性
   - 范围: 0-1
   - 博后特别关注: 磷酸化发生在侧链
   
3. Disorder Propensity (无序倾向)
   - 描述: 该位置形成无序结构的倾向
   - 范围: 0-1 (0=有序, 1=无序)
   - 预测器: DisoMine
   
4. Helix Propensity (α螺旋倾向)
   - 描述: 形成α螺旋的概率
   - 范围: 0-1
   
5. Sheet Propensity (β折叠倾向)
   - 描述: 形成β折叠的概率
   - 范围: 0-1
   
6. Coil Propensity (无规卷曲倾向)
   - 描述: 形成无规卷曲的概率
   - 范围: 0-1
   
7. Early Folding (早期折叠)
   - 描述: 该位置在蛋白质折叠早期形成结构的倾向
   - 范围: 0-1
   - 博后要求新增
   
注: Helix + Sheet + Coil ≈ 1.0

================================================================================
特征分布统计 Feature Distribution
================================================================================

"""

# 添加统计信息
stats_df = success_df[feature_columns].describe().round(3)
report += stats_df.to_string()

report += f"""

================================================================================
数据质量指标 Quality Metrics
================================================================================

覆盖率 Coverage:               {len(success_df)/len(final_df)*100:.1f}%
基因映射成功率:                 {len(gene_to_uniprot)/sites_df['Gene'].nunique()*100:.1f}%

数值范围验证:
"""

for col in feature_columns:
    valid_data = success_df[col].dropna()
    if len(valid_data) > 0:
        in_range = (valid_data.min() >= 0) and (valid_data.max() <= 1)
        status = "PASS" if in_range else "FAIL"
        report += f"  {col:25s}: {status}\n"

# 二级结构一致性
structure_sum = (
    success_df['helix_propensity'].fillna(0) + 
    success_df['sheet_propensity'].fillna(0) + 
    success_df['coil_propensity'].fillna(0)
)
valid_structure = structure_sum[structure_sum > 0.5]
mean_sum = valid_structure.mean()

report += f"""
二级结构一致性:
  Helix+Sheet+Coil 平均值: {mean_sum:.4f} (期望=1.0)
  偏差: {abs(mean_sum - 1.0):.4f}
  状态: {'PASS' if abs(mean_sum - 1.0) < 0.05 else 'WARNING'}

================================================================================
后续分析建议 Next Steps
================================================================================

1. 位点分组 (Site Grouping):
   - 根据亚细胞定位模式对位点进行分组
   - 使用 0min 时间点的 Cyt/Mem/Nuc 列
   
2. 组间比较 (Group Comparison):
   - 比较不同定位模式的位点在生物物理特征上的差异
   - 使用统计检验 (t-test, ANOVA)
   
3. 可视化 (Visualization):
   - 为每个分组绘制特征分布图
   - 比较图显示组间差异
   
4. Movement Score分析:
   - 计算位点的移动性得分
   - 比较高移动性 vs 低移动性位点的生物物理特征
   - 按S/T/Y氨基酸类型拆分分析

推荐使用文件: {success_file}

================================================================================
End of Report
================================================================================
"""

# 保存报告
report_file = OUTPUT_FILE.replace('.csv', '_Report.txt')
with open(report_file, 'w', encoding='utf-8') as f:
    f.write(report)

print(f"\n📄 提取报告已保存:")
print(f"   文件: {report_file}")

print("\n" + "="*80)
print(" 所有文件保存完成！")
print("="*80)
print("\n生成的文件:")
print(f"  1. {OUTPUT_FILE}")
print(f"  2. {success_file}")
if len(failed_df) > 0:
    print(f"  3. {failed_file}")
print(f"  4. {report_file}")
print("\n特征提取流程全部完成！可以进入下一阶段分析。")


# %%
# ===== Cell 10: 定义 A 和 C 的分组函数（修正版）=====

def assign_group_A(row):
    """
    A 方案：时间折叠，根据 Cyt_any/Mem_any/Nuc_any 返回区室组合标签
    """
    c = row['Cyt_any']
    m = row['Mem_any']
    n = row['Nuc_any']
    
    if c == 1 and m == 0 and n == 0:
        return "C-only"
    elif c == 0 and m == 1 and n == 0:
        return "M-only"
    elif c == 0 and m == 0 and n == 1:
        return "N-only"
    elif c == 1 and m == 1 and n == 0:
        return "C&M"
    elif c == 1 and m == 0 and n == 1:
        return "C&N"
    elif c == 0 and m == 1 and n == 1:
        return "M&N"
    elif c == 1 and m == 1 and n == 1:
        return "C&M&N"
    else:
        return "None"


def get_state_at_time(row, time_label):
    """
    C 方案的辅助函数：返回某个时间点的定位状态
    
    注意：列名格式是 "0min_Cyt", "2min_Mem" 等（用下划线）
    """
    c = row[f'{time_label}_Cyt']
    m = row[f'{time_label}_Mem']
    n = row[f'{time_label}_Nuc']
    
    if c == 1 and m == 0 and n == 0:
        return "C"
    elif c == 0 and m == 1 and n == 0:
        return "M"
    elif c == 0 and m == 0 and n == 1:
        return "N"
    elif c == 1 and m == 1 and n == 0:
        return "CM"
    elif c == 1 and m == 0 and n == 1:
        return "CN"
    elif c == 0 and m == 1 and n == 1:
        return "MN"
    elif c == 1 and m == 1 and n == 1:
        return "CMN"
    else:
        return "None"


def assign_group_C(row):
    """
    C 方案：根据 5 个时间点的 trajectory 返回动态类标签
    """
    traj = row['trajectory']
    states = traj.split('→')
    non_none_states = [s for s in states if s != "None"]
    
    if len(non_none_states) == 0:
        return "None"
    
    # sustained-single
    if len(set(states)) == 1 and states[0] in ["C", "M", "N"]:
        return "sustained-single"
    
    # sustained-multi
    if len(set(states)) == 1 and states[0] in ["CN", "CM", "MN", "CMN"]:
        return "sustained-multi"
    
    # early-only
    if len(non_none_states) >= 2 and all(s == "None" for s in states[2:]):
        return "early-only"
    
    # late-appearing
    if all(s == "None" for s in states[:2]) and len(non_none_states) >= 2:
        return "late-appearing"
    
    # transient
    if len(non_none_states) >= 2 and len(non_none_states) < 5:
        if states[0] == "None" or states[-1] == "None":
            return "transient"
    
    # translocation
    if len(set(non_none_states)) > 1:
        return "translocation"
    
    return "other"


print("=" * 60)
print("✓ Cell 10 执行成功！")
print("=" * 60)
print("\n已定义 3 个函数（适配列名格式：下划线连接）")
print("  1. assign_group_A(row)")
print("  2. get_state_at_time(row, time_label)")
print("  3. assign_group_C(row)")
print("=" * 60)


# %%
# ===== Cell 10.5: M1/M2 智能去重（博后要求）=====

import pandas as pd

print("=" * 60)
print("开始执行 Cell 10.5: M1/M2/M3 智能去重")
print("=" * 60)

print("\n博后要求：")
print("  - 如果 M1/M2/M3 的定位模式不同 → 保留所有")
print("  - 如果 M1/M2/M3 的定位模式相同 → 只保留 M1")
print("")

# ========== 1. 定义"定位模式"的比较函数 ==========
def get_localization_signature(row):
    """
    提取该位点的"定位签名"（15 个 0/1 值的元组）
    只要这 15 个值相同，就认为定位模式相同
    """
    time_cols = ['0min_Cyt', '2min_Cyt', '8min_Cyt', '20min_Cyt', '90min_Cyt',
                 '0min_Mem', '2min_Mem', '8min_Mem', '20min_Mem', '90min_Mem',
                 '0min_Nuc', '2min_Nuc', '8min_Nuc', '20min_Nuc', '90min_Nuc']
    return tuple(row[time_cols].values)

# 给每一行加上"定位签名"
success_df['loc_signature'] = success_df.apply(get_localization_signature, axis=1)

# ========== 2. 按 (Gene, Position) 分组处理 ==========
print("[1/4] 识别重复位点...")

grouped = success_df.groupby(['Gene', 'Position'])

# 统计每组有多少条记录
group_sizes = grouped.size()
duplicates = group_sizes[group_sizes > 1]

print(f"   发现 {len(duplicates)} 个位点有多条记录（M1/M2/M3）")
print(f"   涉及总记录数：{duplicates.sum()}")

# ========== 3. 应用去重逻辑 ==========
print("\n[2/4] 应用去重逻辑...")

kept_rows = []
conflict_count = 0
merged_count = 0

for (gene, pos), group in grouped:
    if len(group) == 1:
        # 只有一条记录，直接保留
        kept_rows.append(group.index[0])
    else:
        # 有多条记录（M1/M2/M3），检查定位签名是否相同
        signatures = group['loc_signature'].unique()
        
        if len(signatures) == 1:
            # 定位模式相同 → 只保留 M1
            m1_row = group[group['ModNumber'] == 'M1']
            if len(m1_row) > 0:
                kept_rows.append(m1_row.index[0])
                merged_count += 1
            else:
                # 如果没有 M1（只有 M2/M3），保留第一条
                kept_rows.append(group.index[0])
                merged_count += 1
        else:
            # 定位模式不同 → 保留所有
            kept_rows.extend(group.index.tolist())
            conflict_count += 1

print(f"   相同模式（只保留 M1）：{merged_count} 个位点")
print(f"   不同模式（全部保留）：{conflict_count} 个位点")

# ========== 4. 生成清洗后的数据 ==========
print("\n[3/4] 生成清洗后的数据...")

success_df_cleaned = success_df.loc[kept_rows].copy()

# 标记那些"M1/M2 定位不同"的位点
def mark_conflict(row):
    gene, pos = row['Gene'], row['Position']
    group = success_df[(success_df['Gene'] == gene) & (success_df['Position'] == pos)]
    if len(group) > 1:
        signatures = group['loc_signature'].unique()
        if len(signatures) > 1:
            return 'M1_M2_conflict'
    return 'normal'

success_df_cleaned['M1_M2_status'] = success_df_cleaned.apply(mark_conflict, axis=1)

# 删除临时列
success_df_cleaned = success_df_cleaned.drop(columns=['loc_signature'])

print(f"   清洗前：{len(success_df)} 行")
print(f"   清洗后：{len(success_df_cleaned)} 行")
print(f"   删除重复：{len(success_df) - len(success_df_cleaned)} 行")

# ========== 5. 统计 M1/M2 冲突位点的特征 ==========
print("\n[4/4] M1/M2 冲突位点统计...")

conflict_sites = success_df_cleaned[success_df_cleaned['M1_M2_status'] == 'M1_M2_conflict']
print(f"   发现 {len(conflict_sites)} 个位点存在 M1/M2 定位差异")

if len(conflict_sites) > 0:
    print(f"\n   示例（前 5 个冲突位点）：")
    print(conflict_sites[['Gene', 'Position', 'Site', 'ModNumber']].head(5).to_string(index=False))

# ========== 6. 保存清洗后的数据 ==========
success_df_cleaned.to_csv("success_df_M1M2_cleaned.csv", index=False, encoding='utf-8-sig')
print(f"\n✓ 已保存清洗后的数据：success_df_M1M2_cleaned.csv")

# ========== 7. 更新变量名（让后续 cell 自动使用清洗后的数据）==========
# 重要：把清洗后的数据赋值给 success_df，这样后续 Cell 11A/11C 会自动用清洗后的
success_df = success_df_cleaned.copy()

print("\n" + "=" * 60)
print("✓ Cell 10.5 执行成功！")
print("=" * 60)
print(f"\n清洗后的 success_df：")
print(f"  - 总行数：{len(success_df)}")
print(f"  - 独特位点数：{success_df.groupby(['Gene', 'Position']).ngroups}")
print(f"  - M1/M2 冲突位点：{len(conflict_sites)}")
print("\n后续 Cell 11A-19 将使用清洗后的数据")
print("=" * 60)


# %%
"""
定义LocalizationGroup分组
根据磷酸化位点在不同时间点的细胞区室分布分组
"""

import numpy as np

# 定义细胞区室列
cyt_cols = ['0min_Cyt', '2min_Cyt', '8min_Cyt', '20min_Cyt', '90min_Cyt']
mem_cols = ['0min_Mem', '2min_Mem', '8min_Mem', '20min_Mem', '90min_Mem']
nuc_cols = ['0min_Nuc', '2min_Nuc', '8min_Nuc', '20min_Nuc', '90min_Nuc']

def assign_localization_group(row):
    """
    根据位点在各区室的出现情况分组
    逻辑：如果某个区室在任意时间点有信号（非0/非NaN），则该位点属于该区室
    """
    # 检查每个区室是否有信号
    has_cyt = (row[cyt_cols] > 0).any() if row[cyt_cols].notna().any() else False
    has_mem = (row[mem_cols] > 0).any() if row[mem_cols].notna().any() else False
    has_nuc = (row[nuc_cols] > 0).any() if row[nuc_cols].notna().any() else False
    
    # 分组规则
    if has_cyt and has_mem and has_nuc:
        return 'CMN'
    elif has_cyt and has_mem:
        return 'CM'
    elif has_cyt and has_nuc:
        return 'CN'
    elif has_mem and has_nuc:
        return 'MN'
    elif has_cyt:
        return 'C-only'
    elif has_mem:
        return 'M-only'
    elif has_nuc:
        return 'N-only'
    else:
        return 'Unknown'  # 所有区室都无信号

# 应用分组
print("🔄 正在分配LocalizationGroup...")
cleaned_df['LocalizationGroup'] = cleaned_df.apply(assign_localization_group, axis=1)

# 统计分组结果
print("\n📋 LocalizationGroup分布:")
group_counts = cleaned_df['LocalizationGroup'].value_counts().sort_index()
for group, count in group_counts.items():
    print(f"  {group:12s}: {count:5,} 个位点 ({count/len(cleaned_df)*100:.1f}%)")

print(f"\n✅ 分组完成! 共 {len(cleaned_df):,} 个位点")


# %%
# ===== Cell 11C: C 方案 - 时间轨迹 + 动态分组 =====

print("=" * 60)
print("开始执行 Cell 11C: C 方案分组")
print("=" * 60)

# 1. 对每个 site，提取 5 个时间点的状态
print("\n[1/5] 计算 5 个时间点的状态...")

time_labels = ['0min', '2min', '8min', '20min', '90min']

for t in time_labels:
    success_df[f'state_{t}'] = success_df.apply(lambda row: get_state_at_time(row, t), axis=1)

print(f"   ✓ 新增 5 列：state_0min, state_2min, state_8min, state_20min, state_90min")


# 2. 生成 trajectory 字符串（用 → 连接）
print("\n[2/5] 生成 trajectory 轨迹...")

success_df['trajectory'] = success_df.apply(
    lambda row: '→'.join([row[f'state_{t}'] for t in time_labels]),
    axis=1
)

print(f"   ✓ 新增列：trajectory")

# 查看前 5 个轨迹示例
print("\n   轨迹示例（前 5 行）：")
print(success_df[['Gene', 'Position', 'trajectory']].head(5).to_string(index=False))


# 3. 应用动态分组函数
print("\n[3/5] 应用 assign_group_C 函数...")

success_df['DynamicGroup'] = success_df.apply(assign_group_C, axis=1)

print(f"   ✓ 新增列：DynamicGroup")


# 4. 查看分组分布
print("\n[4/5] 分组分布统计：")
print("-" * 60)

group_counts = success_df['DynamicGroup'].value_counts()
print(group_counts)

print("-" * 60)
print(f"   总计：{group_counts.sum()} sites")


# 5. 过滤掉 "None" 组
print("\n[5/5] 过滤 'None' 组...")

before_count = len(success_df)
analysisdf_C = success_df[success_df['DynamicGroup'] != "None"].copy()
after_count = len(analysisdf_C)

print(f"   过滤前：{before_count} sites")
print(f"   过滤后：{after_count} sites")
print(f"   保留率：{after_count/before_count*100:.1f}%")
print(f"   移除 'None' 组：{before_count - after_count} sites")


# 6. 保存中间结果
print("\n[6/6] 保存结果...")

output_file = "C_GroupedSites.csv"
analysisdf_C.to_csv(output_file, index=False, encoding='utf-8-sig')

print(f"   ✓ 已保存：{output_file}")
print(f"   文件大小：{len(analysisdf_C)} 行 × {len(analysisdf_C.columns)} 列")


# 最终确认
print("\n" + "=" * 60)
print("✓ Cell 11C 执行成功！")
print("=" * 60)
print(f"\n生成变量：")
print(f"  - analysisdf_C ({len(analysisdf_C)} rows)")
print(f"  - 新增列：state_0min, state_2min, ..., trajectory, DynamicGroup")
print(f"\n输出文件：")
print(f"  - C_GroupedSites.csv")
print("=" * 60)


# %%
# ============================================================================
# Cell 12A: Localization Group 特征统计
# ============================================================================

# 检查必要的列是否存在
print("🔍 数据检查:")
print(f"  cleaned_df是否存在? {'cleaned_df' in dir()}")

if 'cleaned_df' not in dir():
    print("  ⚠️  cleaned_df不存在，使用final_df并过滤缺失值")
    cleaned_df = final_df[final_df['backbone_dynamics'].notna()].copy()
    print(f"  ✅ 创建cleaned_df: {len(cleaned_df):,} 行")

print(f"  LocalizationGroup列存在? {'LocalizationGroup' in cleaned_df.columns}")

if 'LocalizationGroup' not in cleaned_df.columns:
    print("\n⚠️  警告: LocalizationGroup列不存在！")
    print("   请先运行定义分组的cell")
    print("   可用的列:")
    print(f"   {cleaned_df.columns.tolist()}")
else:
    # 定义要分析的特征（7个）
    featurecols = ['backbone_dynamics', 'sidechain_dynamics', 'disorder_propensity', 
                   'helix_propensity', 'sheet_propensity', 'coil_propensity', 'earlyFolding']
    
    # 计算统计
    stats_df = cleaned_df.groupby('LocalizationGroup')[featurecols].agg(['mean', 'std', 'count'])
    
    print("\n" + "="*100)
    print("📊 按LocalizationGroup统计 (7个特征)")
    print("="*100)
    print(stats_df.round(3).to_string())
    print("="*100)
    
    # 显示分组样本量
    print("\n📋 各组样本量:")
    group_counts = cleaned_df['LocalizationGroup'].value_counts().sort_index()
    for group, count in group_counts.items():
        print(f"  {group}: {count:,} 个位点")
    
    print("\n✅ 统计完成")


# %%
# ===== Cell 12A: A 方案 - 组统计表 =====

print("=" * 60)
print("开始执行 Cell 12A: A 方案统计表")
print("=" * 60)

# 定义 B2B 特征列（注意：用下划线格式）
# 定义要分析的特征（7个）
featurecols = ['backbone_dynamics', 'sidechain_dynamics', 'disorder_propensity', 
               'helix_propensity', 'sheet_propensity', 'coil_propensity', 'earlyFolding']

# 计算统计
stats_df = cleaned_df.groupby('LocalizationGroup')[featurecols].agg(['mean', 'std', 'count'])

print("📊 按LocalizationGroup统计:")
print(stats_df.round(3))


print("\n[1/2] 计算分组统计...")

# 按 LocalizationGroup 分组，计算统计量
stats_A = analysisdf_A.groupby('LocalizationGroup')[feature_cols].agg(['mean', 'std', 'count'])

# 四舍五入
stats_A = stats_A.round(3)

print("   ✓ 统计完成")

# 2. 显示统计表
print("\n[2/2] A 方案组统计表：")
print("-" * 80)
print(stats_A)
print("-" * 80)

# 3. 保存
output_file = "A_GroupStatistics.csv"
stats_A.to_csv(output_file, encoding='utf-8-sig')

print(f"\n✓ 已保存：{output_file}")

print("\n" + "=" * 60)
print("✓ Cell 12A 执行成功！")
print("=" * 60)


# %%
# ===== Cell 12C: C 方案 - 组统计表 =====

print("=" * 60)
print("开始执行 Cell 12C: C 方案统计表")
print("=" * 60)

print("\n[1/2] 计算分组统计...")

# 按 DynamicGroup 分组，计算统计量
stats_C = analysisdf_C.groupby('DynamicGroup')[feature_cols].agg(['mean', 'std', 'count'])

# 四舍五入
stats_C = stats_C.round(3)

print("   ✓ 统计完成")

# 2. 显示统计表
print("\n[2/2] C 方案组统计表：")
print("-" * 80)
print(stats_C)
print("-" * 80)

# 3. 保存
output_file = "C_GroupStatistics.csv"
stats_C.to_csv(output_file, encoding='utf-8-sig')

print(f"\n✓ 已保存：{output_file}")

print("\n" + "=" * 60)
print("✓ Cell 12C 执行成功！")
print("=" * 60)


# %%
# ===== Cell 13: 快速加载/生成分析数据 (修正版 - 7个特征) =====

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("快速加载/生成分析数据...")
print("="*60)

# 定义路径
BASE_PATH = r"D:\博士\Phospho"
success_file = os.path.join(BASE_PATH, "success_df_M1M2_cleaned.csv")
groupA_file = os.path.join(BASE_PATH, "A_Grouped_Sites.csv")
groupC_file = os.path.join(BASE_PATH, "C_Grouped_Sites.csv")

# ✅ 修正: 定义完整的7个特征列 (包含sidechain_dynamics和earlyFolding)
feature_cols = [
    'backbone_dynamics',       # 主链动力学
    'sidechain_dynamics',      # 侧链动力学 (博后特别关注)
    'disorder_propensity',     # 无序倾向
    'helix_propensity',        # 螺旋倾向
    'sheet_propensity',        # 折叠片倾向
    'coil_propensity',         # 无规卷曲倾向
    'earlyFolding'             # 早期折叠 (博后要求新增)
]

# ===== 1. 加载主数据 =====
print("\n【1/3】加载主数据...")
if os.path.exists(success_file):
    success_df = pd.read_csv(success_file)
    print(f"  ✓ 成功加载: {len(success_df):,} 行")
else:
    print(f"  ✗ 找不到文件: {success_file}")
    raise FileNotFoundError("请先运行生成success_df的代码!")

# ===== 2. 加载或生成方案A数据 =====
print("\n【2/3】处理方案A数据...")
if os.path.exists(groupA_file):
    analysisdf_A = pd.read_csv(groupA_file)
    print(f"  ✓ 加载已有文件: {len(analysisdf_A):,} 行")
else:
    print("  ⚠ 文件不存在,正在生成...")
    
    # 定义分组函数
    def assign_group_A(row):
        c, m, n = row['Cytany'], row['Memany'], row['Nucany']
        if c == 1 and m == 0 and n == 0: return 'C-only'
        elif c == 0 and m == 1 and n == 0: return 'M-only'
        elif c == 0 and m == 0 and n == 1: return 'N-only'
        elif c == 1 and m == 1 and n == 0: return 'CM'
        elif c == 1 and m == 0 and n == 1: return 'CN'
        elif c == 0 and m == 1 and n == 1: return 'MN'
        elif c == 1 and m == 1 and n == 1: return 'CMN'
        else: return None
    
    # 计算time-collapsed any
    success_df['Cytany'] = success_df[['0min_Cyt', '2min_Cyt', '8min_Cyt', '20min_Cyt', '90min_Cyt']].max(axis=1)
    success_df['Memany'] = success_df[['0min_Mem', '2min_Mem', '8min_Mem', '20min_Mem', '90min_Mem']].max(axis=1)
    success_df['Nucany'] = success_df[['0min_Nuc', '2min_Nuc', '8min_Nuc', '20min_Nuc', '90min_Nuc']].max(axis=1)
    
    # 分配分组
    success_df['LocalizationGroup'] = success_df.apply(assign_group_A, axis=1)
    
    # 过滤掉None
    analysisdf_A = success_df[success_df['LocalizationGroup'].notna()].copy()
    
    # 保存
    analysisdf_A.to_csv(groupA_file, index=False, encoding='utf-8-sig')
    print(f"  ✓ 生成并保存: {len(analysisdf_A):,} 行")
    print(f"  分组分布:\n{analysisdf_A['LocalizationGroup'].value_counts()}")

# ===== 3. 加载或生成方案C数据 =====
print("\n【3/3】处理方案C数据...")
if os.path.exists(groupC_file):
    analysisdf_C = pd.read_csv(groupC_file)
    print(f"  ✓ 加载已有文件: {len(analysisdf_C):,} 行")
else:
    print("  ⚠ 文件不存在,正在生成...")
    
    # 定义状态获取函数
    def get_state_at_time(row, timelabel):
        c = row[f'{timelabel}_Cyt']
        m = row[f'{timelabel}_Mem']
        n = row[f'{timelabel}_Nuc']
        if c == 1 and m == 0 and n == 0: return 'C'
        elif c == 0 and m == 1 and n == 0: return 'M'
        elif c == 0 and m == 0 and n == 1: return 'N'
        elif c == 1 and m == 1 and n == 0: return 'CM'
        elif c == 1 and m == 0 and n == 1: return 'CN'
        elif c == 0 and m == 1 and n == 1: return 'MN'
        elif c == 1 and m == 1 and n == 1: return 'CMN'
        else: return None
    
    # 定义分组函数C
    def assign_group_C(row):
        traj = row['trajectory']
        states = [s for s in traj if s != 'None']
        nonnone_states = [s for s in states if s is not None]
        
        if len(nonnone_states) == 0:
            return None
        
        if len(set(states)) == 1 and states[0] in ['C', 'M', 'N']:
            return 'sustained-single'
        
        if len(set(states)) == 1 and states[0] in ['CN', 'CM', 'MN', 'CMN']:
            return 'sustained-multi'
        
        if len(set(nonnone_states)) > 1:
            return 'translocation'
        
        return 'other'
    
    # 生成5个时间点的状态
    time_labels = ['0min', '2min', '8min', '20min', '90min']
    for t in time_labels:
        success_df[f'state_{t}'] = success_df.apply(lambda row: get_state_at_time(row, t), axis=1)
    
    # 生成trajectory字符串
    success_df['trajectory'] = success_df.apply(
        lambda row: '-'.join([str(row[f'state_{t}']) for t in time_labels]), axis=1
    )
    
    # 分配分组
    success_df['DynamicGroup'] = success_df.apply(assign_group_C, axis=1)
    
    # 过滤掉None和other
    analysisdf_C = success_df[success_df['DynamicGroup'].notna()].copy()
    analysisdf_C = analysisdf_C[analysisdf_C['DynamicGroup'] != 'other']
    
    # 保存
    analysisdf_C.to_csv(groupC_file, index=False, encoding='utf-8-sig')
    print(f"  ✓ 生成并保存: {len(analysisdf_C):,} 行")
    print(f"  分组分布:\n{analysisdf_C['DynamicGroup'].value_counts()}")

print("\n" + "="*60)
print("✅ 所有数据准备完毕!")
print("="*60)
print(f"\n当前变量:")
print(f"  • success_df: {len(success_df):,} 行(原始数据)")
print(f"  • analysisdf_A: {len(analysisdf_A):,} 行(方案A)")
print(f"  • analysisdf_C: {len(analysisdf_C):,} 行(方案C)")
print(f"  • feature_cols: {feature_cols}")
print(f"\n✅ 特征列表已更新为7个特征(包含sidechain_dynamics和earlyFolding)")


# %%
# ===== Cell 13C: C 方案 - 箱线图 =====

print("=" * 60)
print("开始执行 Cell 13C: C 方案箱线图")
print("=" * 60)

# 定义组顺序（过滤掉样本数太少的组，比如 translocation 只有 1 个）
# 只绘制样本数 >= 10 的组
group_counts_C = analysisdf_C['DynamicGroup'].value_counts()
valid_groups_C = group_counts_C[group_counts_C >= 10].index.tolist()

print(f"\n有效分组（样本数 >= 10）：{valid_groups_C}")

if len(valid_groups_C) == 0:
    print("⚠️  警告：没有样本数 >= 10 的组，跳过绘图")
else:
    # 过滤数据
    analysisdf_C_filtered = analysisdf_C[analysisdf_C['DynamicGroup'].isin(valid_groups_C)].copy()
    
    print(f"\n准备绘制 {len(feature_cols)} 个特征的箱线图...")
    
    for i, feature in enumerate(feature_cols, 1):
        print(f"\n[{i}/{len(feature_cols)}] 绘制 {feature}...")
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        sns.boxplot(data=analysisdf_C_filtered, x='DynamicGroup', y=feature, 
                    order=valid_groups_C, palette='viridis', ax=ax)
        
        ax.set_title(f"C: {feature} by Dynamic Group", fontsize=14, fontweight='bold')
        ax.set_xlabel("Dynamic Group", fontsize=12)
        ax.set_ylabel(feature.replace('_', ' ').title(), fontsize=12)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        figfile = f"C_{feature}_boxplot.png"
        plt.savefig(figfile, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"   ✓ 已保存：{figfile}")
        
        plt.close()

print("\n" + "=" * 60)
print("✓ Cell 13C 执行成功！")
print("=" * 60)
print(f"\n生成了 {len(feature_cols)} 张箱线图：")
for feature in feature_cols:
    print(f"  - C_{feature}_boxplot.png")
print("=" * 60)


# %%
# ===== Cell 14: 筛选代表蛋白（用于单蛋白故事图）=====

print("=" * 60)
print("开始执行 Cell 14: 筛选代表蛋白")
print("=" * 60)

# ========== A 方案的筛选 ==========
print("\n[A 方案] 筛选候选蛋白...")
print("-" * 60)

# 筛选标准：
# 1. 至少有 5 个 sites
# 2. 这些 sites 覆盖至少 3 个不同的 LocalizationGroup
# 3. 蛋白至少有 1 个 site（用于后续检查蛋白长度）

candidate_A = analysisdf_A.groupby('Gene').filter(
    lambda x: len(x) >= 5 and x['LocalizationGroup'].nunique() >= 3
)

protein_list_A = candidate_A.groupby('Gene').agg({
    'Position': 'count',  # site 数
    'LocalizationGroup': lambda x: x.nunique()  # 组数
}).rename(columns={'Position': 'num_sites', 'LocalizationGroup': 'num_groups'})

protein_list_A = protein_list_A.sort_values(['num_groups', 'num_sites'], ascending=False)

print(f"✓ 找到 {len(protein_list_A)} 个候选蛋白")
print("\nA 方案候选蛋白（前 10）：")
print(protein_list_A.head(10))


# ========== C 方案的筛选 ==========
print("\n\n[C 方案] 筛选候选蛋白...")
print("-" * 60)

# 注意：C 方案只有 2 个有效组（sustained-single 和 sustained-multi）
# 所以筛选标准调整为：至少 5 个 sites，覆盖至少 2 个组

candidate_C = analysisdf_C.groupby('Gene').filter(
    lambda x: len(x) >= 5 and x['DynamicGroup'].nunique() >= 2
)

protein_list_C = candidate_C.groupby('Gene').agg({
    'Position': 'count',
    'DynamicGroup': lambda x: x.nunique()
}).rename(columns={'Position': 'num_sites', 'DynamicGroup': 'num_groups'})

protein_list_C = protein_list_C.sort_values(['num_groups', 'num_sites'], ascending=False)

print(f"✓ 找到 {len(protein_list_C)} 个候选蛋白")
print("\nC 方案候选蛋白（前 10）：")
print(protein_list_C.head(10))


# ========== 找交集（同时满足 A 和 C）==========
print("\n\n[交集] 同时满足 A 和 C 的候选蛋白...")
print("-" * 60)

common_candidates = set(protein_list_A.index) & set(protein_list_C.index)

if len(common_candidates) > 0:
    print(f"✓ 找到 {len(common_candidates)} 个同时满足 A 和 C 的蛋白")
    
    # 合并信息
    common_df = pd.DataFrame({
        'Gene': list(common_candidates),
        'A_num_sites': [protein_list_A.loc[g, 'num_sites'] for g in common_candidates],
        'A_num_groups': [protein_list_A.loc[g, 'num_groups'] for g in common_candidates],
        'C_num_sites': [protein_list_C.loc[g, 'num_sites'] for g in common_candidates],
        'C_num_groups': [protein_list_C.loc[g, 'num_groups'] for g in common_candidates]
    })
    
    # 按 A 的组数 + C 的组数排序
    common_df['total_groups'] = common_df['A_num_groups'] + common_df['C_num_groups']
    common_df = common_df.sort_values('total_groups', ascending=False)
    
    print("\n推荐蛋白（前 10）：")
    print(common_df.head(10).to_string(index=False))
    
    # 保存候选蛋白列表
    common_df.to_csv("CandidateProteins_for_StoryFigure.csv", index=False, encoding='utf-8-sig')
    print("\n✓ 已保存：CandidateProteins_for_StoryFigure.csv")
    
else:
    print("⚠️  没有找到同时满足 A 和 C 的蛋白")
    print("   建议：从 A 或 C 的候选列表中单独选择")


# ========== 给出推荐 ==========
print("\n\n" + "=" * 60)
print("推荐下一步操作：")
print("=" * 60)

if len(common_candidates) > 0:
    # 推荐第一名
    top_gene = common_df.iloc[0]['Gene']
    print(f"\n🌟 推荐使用蛋白：{top_gene}")
    print(f"   - A 方案：{common_df.iloc[0]['A_num_sites']:.0f} sites, {common_df.iloc[0]['A_num_groups']:.0f} groups")
    print(f"   - C 方案：{common_df.iloc[0]['C_num_sites']:.0f} sites, {common_df.iloc[0]['C_num_groups']:.0f} groups")
    print(f"\n在 Cell 15 中，将 'selected_gene' 设置为 '{top_gene}'")
else:
    print("\n请从以下列表中手动选择一个蛋白：")
    print(f"   - A 方案推荐：{protein_list_A.index[0]}")
    print(f"   - C 方案推荐：{protein_list_C.index[0]}")

print("\n" + "=" * 60)
print("✓ Cell 14 执行成功！")
print("=" * 60)


# %%
# ===== Cell 15: 单蛋白故事图（手绘图）=====

import json
import matplotlib.pyplot as plt
import os

print("=" * 60)
print("开始执行 Cell 15: 单蛋白故事图")
print("=" * 60)

# ========== 1. 选择蛋白 ==========
selected_gene = "SRRM2"  # 从 Cell 14 的推荐

print(f"\n[1/6] 选择蛋白：{selected_gene}")


# ========== 2. 检查并加载必要变量 ==========
print("\n[2/6] 检查必要变量...")

# 检查 gene_to_uniprot 是否存在
try:
    test = gene_to_uniprot[selected_gene]
    print(f"   ✓ gene_to_uniprot 已存在")
except (NameError, KeyError) as e:
    print(f"   ⚠️  gene_to_uniprot 不存在或缺少 {selected_gene}，需要重新生成")
    print("   提示：需要运行之前生成 gene_to_uniprot 的 cell")
    raise

# 检查 JSON_FOLDER 是否存在
try:
    test = JSON_FOLDER
    print(f"   ✓ JSON_FOLDER = {JSON_FOLDER}")
except NameError:
    print("   ⚠️  JSON_FOLDER 不存在，使用默认值")
    JSON_FOLDER = "D:/json"  # 根据你的实际路径修改

selected_uniprot = gene_to_uniprot[selected_gene]
print(f"   ✓ Uniprot ID: {selected_uniprot}")


# ========== 3. 读取 JSON 文件 ==========
print(f"\n[3/6] 读取 {selected_gene} 的 B2B 特征...")

json_path = os.path.join(JSON_FOLDER, f"{selected_uniprot}.json")

if not os.path.exists(json_path):
    print(f"   ❌ 错误：找不到文件 {json_path}")
    raise FileNotFoundError(f"JSON 文件不存在：{json_path}")

with open(json_path, 'r', encoding='utf-8') as f:
    protein_data = json.load(f)
    residues = protein_data['residues']

residue_df = pd.DataFrame(residues)

print(f"   ✓ 读取成功：{len(residue_df)} 个氨基酸")
print(f"   蛋白长度：{residue_df['seqpos'].min()} - {residue_df['seqpos'].max()}")


# ========== 4. 提取该蛋白的所有 sites（A 和 C）==========
print(f"\n[4/6] 提取 {selected_gene} 的磷酸化位点...")

protein_sites_A = analysisdf_A[analysisdf_A['Gene'] == selected_gene].copy()
protein_sites_C = analysisdf_C[analysisdf_C['Gene'] == selected_gene].copy()

print(f"   A 方案：{len(protein_sites_A)} sites")
print(f"   C 方案：{len(protein_sites_C)} sites")

# 查看 A 的分组分布
print(f"\n   A 方案分组分布：")
print(protein_sites_A['LocalizationGroup'].value_counts())


# ========== 5. 画 A 版本的图 ==========
print(f"\n[5/6] 绘制 A 版本（区室组合）...")

fig, ax = plt.subplots(figsize=(16, 6))

# 画 5 条 B2B 曲线
ax.plot(residue_df['seqpos'], residue_df['backbone'], label='Backbone dynamics', 
        color='black', linewidth=1.5, alpha=0.8)
ax.plot(residue_df['seqpos'], residue_df['disoMine'], label='Disorder', 
        color='gray', linewidth=1.5, alpha=0.8)
ax.plot(residue_df['seqpos'], residue_df['helix'], label='Helix', 
        color='red', linestyle='--', linewidth=1, alpha=0.7)
ax.plot(residue_df['seqpos'], residue_df['sheet'], label='Sheet', 
        color='blue', linestyle='--', linewidth=1, alpha=0.7)
ax.plot(residue_df['seqpos'], residue_df['coil'], label='Coil', 
        color='green', linestyle='--', linewidth=1, alpha=0.7)

# 定义 A 的颜色映射
color_map_A = {
    "C-only": "#E74C3C",    # 红色
    "N-only": "#3498DB",    # 蓝色
    "M-only": "#2ECC71",    # 绿色
    "C&N": "#9B59B6",       # 紫色
    "C&M": "#F39C12",       # 橙色
    "M&N": "#1ABC9C",       # 青色
    "C&M&N": "#E91E63"      # 品红
}

# 标注 sites
for _, site in protein_sites_A.iterrows():
    pos = site['Position']
    group = site['LocalizationGroup']
    color = color_map_A.get(group, 'black')
    
    # 画垂直虚线
    ax.axvline(x=pos, color=color, alpha=0.5, linewidth=1.5, linestyle=':')

# 设置图例（合并曲线和位点颜色）
handles, labels = ax.get_legend_handles_labels()

# 添加位点颜色图例
from matplotlib.patches import Patch
for group in protein_sites_A['LocalizationGroup'].unique():
    handles.append(Patch(color=color_map_A.get(group, 'black'), label=group))

ax.legend(handles=handles, loc='upper right', fontsize=9, ncol=2)

ax.set_xlabel("Sequence position", fontsize=12)
ax.set_ylabel("B2B feature value", fontsize=12)
ax.set_title(f"Biophysical landscape of {selected_gene} (A: Localization Groups)", 
             fontsize=14, fontweight='bold')
ax.set_ylim(-0.2, 1.1)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
figfile_A = f"A_{selected_gene}_FeatureMap.png"
plt.savefig(figfile_A, dpi=300, bbox_inches='tight', facecolor='white')
print(f"   ✓ 已保存：{figfile_A}")
plt.close()


# ========== 6. 画 C 版本的图 ==========
print(f"\n[6/6] 绘制 C 版本（动态组）...")

fig, ax = plt.subplots(figsize=(16, 6))

# 画 5 条 B2B 曲线（与 A 相同）
ax.plot(residue_df['seqpos'], residue_df['backbone'], label='Backbone dynamics', 
        color='black', linewidth=1.5, alpha=0.8)
ax.plot(residue_df['seqpos'], residue_df['disoMine'], label='Disorder', 
        color='gray', linewidth=1.5, alpha=0.8)
ax.plot(residue_df['seqpos'], residue_df['helix'], label='Helix', 
        color='red', linestyle='--', linewidth=1, alpha=0.7)
ax.plot(residue_df['seqpos'], residue_df['sheet'], label='Sheet', 
        color='blue', linestyle='--', linewidth=1, alpha=0.7)
ax.plot(residue_df['seqpos'], residue_df['coil'], label='Coil', 
        color='green', linestyle='--', linewidth=1, alpha=0.7)

# 定义 C 的颜色映射
color_map_C = {
    "sustained-single": "#2C3E50",  # 深蓝灰
    "sustained-multi": "#E74C3C",   # 红色
    "translocation": "#F39C12"      # 橙色
}

# 标注 sites
for _, site in protein_sites_C.iterrows():
    pos = site['Position']
    group = site['DynamicGroup']
    color = color_map_C.get(group, 'black')
    
    ax.axvline(x=pos, color=color, alpha=0.5, linewidth=1.5, linestyle=':')

# 设置图例
handles, labels = ax.get_legend_handles_labels()
for group in protein_sites_C['DynamicGroup'].unique():
    handles.append(Patch(color=color_map_C.get(group, 'black'), label=group))

ax.legend(handles=handles, loc='upper right', fontsize=9, ncol=2)

ax.set_xlabel("Sequence position", fontsize=12)
ax.set_ylabel("B2B feature value", fontsize=12)
ax.set_title(f"Biophysical landscape of {selected_gene} (C: Dynamic Groups)", 
             fontsize=14, fontweight='bold')
ax.set_ylim(-0.2, 1.1)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
figfile_C = f"C_{selected_gene}_FeatureMap.png"
plt.savefig(figfile_C, dpi=300, bbox_inches='tight', facecolor='white')
print(f"   ✓ 已保存：{figfile_C}")
plt.close()


# ========== 完成 ==========
print("\n" + "=" * 60)
print("✓ Cell 15 执行成功！")
print("=" * 60)
print(f"\n生成的单蛋白故事图：")
print(f"  1. {figfile_A}  (A: 区室组合)")
print(f"  2. {figfile_C}  (C: 动态组)")
print("\n这两张图展示了 {selected_gene} 蛋白的完整 B2B landscape")
print("以及 170 个磷酸化位点在不同分组下的分布。")
print("=" * 60)


# %%
# ===== Cell 16: 统计检验 - 组间显著性分析 =====

from scipy import stats
import pandas as pd

print("=" * 60)
print("开始执行 Cell 16: 统计显著性检验")
print("=" * 60)

# ========== A 方案的统计检验 ==========
print("\n[A 方案] ANOVA + Kruskal-Wallis 检验")
print("-" * 60)

# 定义特征列
feature_cols = ['backbone_dynamics', 'disorder_propensity', 'helix_propensity', 
                'sheet_propensity', 'coil_propensity']

results_A = []

for feature in feature_cols:
    print(f"\n{feature}:")
    
    # 准备数据：按组分割
    groups = []
    group_names = []
    for group in analysisdf_A['LocalizationGroup'].unique():
        data = analysisdf_A[analysisdf_A['LocalizationGroup'] == group][feature].dropna()
        if len(data) >= 3:  # 至少 3 个样本才纳入检验
            groups.append(data)
            group_names.append(group)
    
    if len(groups) < 2:
        print(f"   ⚠️  样本不足，跳过")
        continue
    
    # 1. ANOVA（参数检验，假设正态分布）
    f_stat, p_anova = stats.f_oneway(*groups)
    
    # 2. Kruskal-Wallis（非参数检验，不假设正态分布）
    h_stat, p_kruskal = stats.kruskal(*groups)
    
    print(f"   ANOVA:          F = {f_stat:.2f}, p = {p_anova:.2e} {'***' if p_anova < 0.001 else '**' if p_anova < 0.01 else '*' if p_anova < 0.05 else 'ns'}")
    print(f"   Kruskal-Wallis: H = {h_stat:.2f}, p = {p_kruskal:.2e} {'***' if p_kruskal < 0.001 else '**' if p_kruskal < 0.01 else '*' if p_kruskal < 0.05 else 'ns'}")
    
    results_A.append({
        'Feature': feature,
        'Test': 'ANOVA',
        'Statistic': f_stat,
        'p-value': p_anova,
        'Significant': p_anova < 0.05
    })
    
    results_A.append({
        'Feature': feature,
        'Test': 'Kruskal-Wallis',
        'Statistic': h_stat,
        'p-value': p_kruskal,
        'Significant': p_kruskal < 0.05
    })

# 转成 DataFrame
results_A_df = pd.DataFrame(results_A)

print("\n\nA 方案统计检验汇总：")
print("-" * 60)
print(results_A_df.to_string(index=False))


# ========== C 方案的统计检验 ==========
print("\n\n[C 方案] 统计检验")
print("-" * 60)

# C 方案只有 2-3 个组，用 t-test 或 Mann-Whitney
valid_groups_C = analysisdf_C['DynamicGroup'].value_counts()
valid_groups_C = valid_groups_C[valid_groups_C >= 10].index.tolist()

if len(valid_groups_C) >= 2:
    results_C = []
    
    for feature in feature_cols:
        print(f"\n{feature}:")
        
        groups = []
        for group in valid_groups_C:
            data = analysisdf_C[analysisdf_C['DynamicGroup'] == group][feature].dropna()
            groups.append(data)
        
        if len(groups) == 2:
            # t-test（参数）
            t_stat, p_ttest = stats.ttest_ind(groups[0], groups[1])
            
            # Mann-Whitney（非参数）
            u_stat, p_mann = stats.mannwhitneyu(groups[0], groups[1], alternative='two-sided')
            
            print(f"   t-test:        t = {t_stat:.2f}, p = {p_ttest:.2e} {'***' if p_ttest < 0.001 else '**' if p_ttest < 0.01 else '*' if p_ttest < 0.05 else 'ns'}")
            print(f"   Mann-Whitney:  U = {u_stat:.0f}, p = {p_mann:.2e} {'***' if p_mann < 0.001 else '**' if p_mann < 0.01 else '*' if p_mann < 0.05 else 'ns'}")
            
            results_C.append({
                'Feature': feature,
                'Test': 't-test',
                'Statistic': t_stat,
                'p-value': p_ttest,
                'Significant': p_ttest < 0.05
            })
            
            results_C.append({
                'Feature': feature,
                'Test': 'Mann-Whitney',
                'Statistic': u_stat,
                'p-value': p_mann,
                'Significant': p_mann < 0.05
            })
        
        elif len(groups) > 2:
            # Kruskal-Wallis
            h_stat, p_kruskal = stats.kruskal(*groups)
            print(f"   Kruskal-Wallis: H = {h_stat:.2f}, p = {p_kruskal:.2e} {'***' if p_kruskal < 0.001 else '**' if p_kruskal < 0.01 else '*' if p_kruskal < 0.05 else 'ns'}")
            
            results_C.append({
                'Feature': feature,
                'Test': 'Kruskal-Wallis',
                'Statistic': h_stat,
                'p-value': p_kruskal,
                'Significant': p_kruskal < 0.05
            })
    
    results_C_df = pd.DataFrame(results_C)
    
    print("\n\nC 方案统计检验汇总：")
    print("-" * 60)
    print(results_C_df.to_string(index=False))
else:
    print("   ⚠️  C 方案有效组不足 2 个，跳过统计检验")
    results_C_df = pd.DataFrame()


# ========== 保存结果 ==========
print("\n\n保存统计检验结果...")

results_A_df.to_csv("A_StatisticalTests.csv", index=False, encoding='utf-8-sig')
print(f"   ✓ 已保存：A_StatisticalTests.csv")

if not results_C_df.empty:
    results_C_df.to_csv("C_StatisticalTests.csv", index=False, encoding='utf-8-sig')
    print(f"   ✓ 已保存：C_StatisticalTests.csv")


# ========== 结论 ==========
print("\n" + "=" * 60)
print("结论与建议")
print("=" * 60)

# 统计 A 方案有多少个特征显著
sig_count_A = results_A_df[results_A_df['Test'] == 'ANOVA']['Significant'].sum()

print(f"\nA 方案：{sig_count_A}/{len(feature_cols)} 个特征显示显著组间差异")

if sig_count_A >= 3:
    print("✅ **建议：不需要做 ±5 氨基酸扩展分析**")
    print("   已有足够显著的 pattern，可以直接用现有结果讲故事。")
elif sig_count_A >= 1:
    print("⚠️  **建议：部分特征显著，可选择性做 ±5 扩展**")
    print("   对不显著的特征做扩展分析可能增强信号。")
else:
    print("❌ **建议：需要做 ±5 氨基酸扩展分析**")
    print("   位点本身未显示显著差异，需要扩展到邻近氨基酸。")

if not results_C_df.empty:
    sig_count_C = results_C_df['Significant'].sum() / 2  # 每个特征有 2 个检验
    print(f"\nC 方案：约 {sig_count_C:.0f}/{len(feature_cols)} 个特征显示显著差异")

print("\n" + "=" * 60)
print("✓ Cell 16 执行成功！")
print("=" * 60)


# %%
# ===== Cell 17: 数据驱动的分组策略 - 热图分析 =====

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

print("=" * 60)
print("开始执行 Cell 18: 热图分析各组差异")
print("=" * 60)

# 提取 3 个关键特征
key_features = ['backbone_dynamics', 'disorder_propensity', 'helix_propensity']

# 计算各组的均值
group_means = analysisdf_A.groupby('LocalizationGroup')[key_features].mean()

# 标准化（Z-score），让差异更明显
from scipy.stats import zscore
group_means_zscore = group_means.apply(zscore, axis=0)

# 定义组顺序（按逻辑排序）
group_order = ["C-only", "N-only", "M-only", "C&N", "C&M", "M&N", "C&M&N"]
group_means_zscore = group_means_zscore.loc[group_order]

# 绘制热图
fig, ax = plt.subplots(figsize=(8, 6))

sns.heatmap(group_means_zscore, 
            annot=True,  # 显示数值
            fmt='.2f',   # 保留 2 位小数
            cmap='RdBu_r',  # 红-蓝配色（红=高，蓝=低）
            center=0,    # 0 为中心（白色）
            linewidths=1,
            cbar_kws={'label': 'Z-score (Standardized Value)'},
            ax=ax)

ax.set_title("Biophysical Feature Comparison Across Localization Groups\n(Higher value = More extreme)", 
             fontsize=12, fontweight='bold')
ax.set_xlabel("Biophysical Features", fontsize=11)
ax.set_ylabel("Localization Groups", fontsize=11)

plt.tight_layout()
figfile = "A_GroupComparison_Heatmap.png"
plt.savefig(figfile, dpi=300, bbox_inches='tight', facecolor='white')
print(f"✓ 已保存：{figfile}")

plt.show()

# 数值解读
print("\n" + "=" * 60)
print("数据解读")
print("=" * 60)

print("\n原始均值（未标准化）：")
print(group_means.round(3))

print("\n标准化 Z-score（正值=高于平均，负值=低于平均）：")
print(group_means_zscore.round(2))

print("\n" + "=" * 60)
print("推荐标注策略（基于数据）：")
print("=" * 60)

# 找出每个特征最极端的组
for feature in key_features:
    max_group = group_means_zscore[feature].idxmax()
    min_group = group_means_zscore[feature].idxmin()
    print(f"\n{feature}:")
    print(f"  最高：{max_group} (Z={group_means_zscore.loc[max_group, feature]:.2f})")
    print(f"  最低：{min_group} (Z={group_means_zscore.loc[min_group, feature]:.2f})")

print("\n" + "-" * 60)
print("建议突出显示的组（Z-score 绝对值 > 0.5）：")
extreme_groups = set()
for feature in key_features:
    for group in group_means_zscore.index:
        if abs(group_means_zscore.loc[group, feature]) > 0.5:
            extreme_groups.add(group)


extreme_groups_list = sorted(extreme_groups)  # 先保存排序后的版本用于打印
print(f"  {', '.join(extreme_groups_list)}")

print("\n建议灰化/合并的组（Z-score 绝对值都 < 0.5）：")
mild_groups = set(group_means_zscore.index) - extreme_groups  # ← extreme_groups 仍是 set
mild_groups = sorted(mild_groups)
print(f"  {', '.join(mild_groups)}")


# %%
# ===== Cell 19: 自动选择要在3D结构上标注的位点 =====

import pandas as pd

print("=" * 60)
print("开始执行 Cell 19: 选择3D结构标注位点")
print("=" * 60)

# ========== 1. 筛选SRRM2的位点 ==========
srrm2_sites = analysisdf_A[analysisdf_A['Gene'] == 'SRRM2'].copy()
print(f"\nSRRM2 总共有 {len(srrm2_sites)} 个磷酸化位点")

# ========== 2. 选择M-only组的代表位点 ==========
print("\n[M-only组] 选择标注位点...")

m_only = srrm2_sites[srrm2_sites['LocalizationGroup'] == 'M-only'].copy()
print(f"  M-only组共 {len(m_only)} 个位点")

# 按helix_propensity降序排列（选最"螺旋化"的）
m_only = m_only.sort_values('helix_propensity', ascending=False)

# 强制空间分散：每个位点至少相隔200 aa
selected_m_only = []
for _, row in m_only.iterrows():
    pos = row['Position']
    if len(selected_m_only) == 0:
        selected_m_only.append(row)
    elif all(abs(pos - s['Position']) > 200 for s in selected_m_only):
        selected_m_only.append(row)
    if len(selected_m_only) >= 3:
        break

print(f"  选中 {len(selected_m_only)} 个位点：")
for site in selected_m_only:
    print(f"    - Position {site['Position']:4d} ({site['Site']:8s}): Helix={site['helix_propensity']:.3f}, Backbone={site['backbone_dynamics']:.3f}")

# ========== 3. 选择C&M&N组的代表位点 ==========
print("\n[C&M&N组] 选择标注位点...")

cmn = srrm2_sites[srrm2_sites['LocalizationGroup'] == 'C&M&N'].copy()
print(f"  C&M&N组共 {len(cmn)} 个位点")

# 按disorder_propensity降序排列（选最"无序"的）
cmn = cmn.sort_values('disorder_propensity', ascending=False)

# 强制空间分散
selected_cmn = []
for _, row in cmn.iterrows():
    pos = row['Position']
    if len(selected_cmn) == 0:
        selected_cmn.append(row)
    elif all(abs(pos - s['Position']) > 200 for s in selected_cmn):
        selected_cmn.append(row)
    if len(selected_cmn) >= 3:
        break

print(f"  选中 {len(selected_cmn)} 个位点：")
for site in selected_cmn:
    print(f"    - Position {site['Position']:4d} ({site['Site']:8s}): Disorder={site['disorder_propensity']:.3f}, Backbone={site['backbone_dynamics']:.3f}")

# ========== 4. 生成PyMOL脚本 ==========
print("\n[生成PyMOL脚本]...")

pml_script = f"""# PyMOL Script for SRRM2 Phosphosite Visualization
# Generated automatically

# Load structure (replace with your actual PDB file path)
load SRRM2_AlphaFold.pdb, SRRM2

# Basic settings
bg_color white
set cartoon_transparency, 0.3
set surface_color, gray80
show surface
hide cartoon

# Color M-only sites (blue)
"""

for site in selected_m_only:
    pml_script += f"select m_site_{site['Position']}, resi {site['Position']}\n"
    pml_script += f"color blue, m_site_{site['Position']}\n"
    pml_script += f"show spheres, m_site_{site['Position']}\n"
    pml_script += f"label m_site_{site['Position']} and name CA, '{site['Site']}'\n"

pml_script += "\n# Color C&M&N sites (red)\n"

for site in selected_cmn:
    pml_script += f"select cmn_site_{site['Position']}, resi {site['Position']}\n"
    pml_script += f"color red, cmn_site_{site['Position']}\n"
    pml_script += f"show spheres, cmn_site_{site['Position']}\n"
    pml_script += f"label cmn_site_{site['Position']} and name CA, '{site['Site']}'\n"

pml_script += """
# Adjust view
set label_color, black
set label_size, 20
zoom
"""

# 保存脚本
script_file = "SRRM2_phosphosite_visualization.pml"
with open(script_file, 'w') as f:
    f.write(pml_script)

print(f"  ✓ 已保存PyMOL脚本：{script_file}")

# ========== 5. 保存选中的位点列表 ==========
selected_sites_df = pd.DataFrame(selected_m_only + selected_cmn)
selected_sites_df.to_csv("Selected_Sites_for_3D.csv", index=False, encoding='utf-8-sig')
print(f"  ✓ 已保存位点列表：Selected_Sites_for_3D.csv")

print("\n" + "=" * 60)
print("下一步操作指南")
print("=" * 60)
print("\n1. 下载SRRM2的AlphaFold结构：")
print("   https://alphafold.ebi.ac.uk/")
print("   搜索 'SRRM2' 或 Uniprot ID，下载PDB文件")
print("\n2. 安装PyMOL（教育版免费）：")
print("   https://pymol.org/edu/")
print("\n3. 打开PyMOL，运行脚本：")
print(f"   File → Run Script → 选择 {script_file}")
print("\n4. 调整视角后截图：")
print("   File → Export Image as PNG → 选择分辨率（建议3000x2000）")
print("\n" + "=" * 60)
print("✓ Cell 19 执行成功！")
print("=" * 60)


# %%
# ===== Cell 15A: 单蛋白故事图（三极对比 + 优化版）=====

import json
import matplotlib.pyplot as plt
import os
import numpy as np

print("=" * 60)
print("开始执行 Cell 15A: 单蛋白故事图（优化版）")
print("=" * 60)

# ========== 1. 选择蛋白 ==========
selected_gene = "SRRM2"
print(f"\n[1/5] 选择蛋白：{selected_gene}")

# ========== 2. 读取 JSON 文件 ==========
print(f"\n[2/5] 读取 {selected_gene} 的 B2B 特征...")

selected_uniprot = gene_to_uniprot[selected_gene]
json_path = os.path.join(JSON_FOLDER, f"{selected_uniprot}.json")

with open(json_path, 'r', encoding='utf-8') as f:
    protein_data = json.load(f)
    residues = protein_data['residues']

residue_df = pd.DataFrame(residues)

protein_length = residue_df['seqpos'].max()
print(f"   ✓ 蛋白总长度：{protein_length} aa")

# ========== 3. 决定展示范围（Zoom In）==========
print(f"\n[3/5] 决定展示范围...")

# 策略：找位点最密集的区域
protein_sites = analysisdf_A[analysisdf_A['Gene'] == selected_gene].copy()

# 将序列切成 500 aa 的窗口，找位点最多的那个窗口
window_size = 600
best_start = 0
max_sites_in_window = 0

for start in range(0, protein_length - window_size, 100):
    end = start + window_size
    sites_in_window = len(protein_sites[(protein_sites['Position'] >= start) & 
                                         (protein_sites['Position'] <= end)])
    if sites_in_window > max_sites_in_window:
        max_sites_in_window = sites_in_window
        best_start = start

best_end = best_start + window_size

print(f"   选择展示区域：aa {best_start}-{best_end}")
print(f"   该区域包含 {max_sites_in_window} 个位点（占总数 {max_sites_in_window/len(protein_sites)*100:.1f}%）")

# ========== 4. 筛选三个极端组的位点 ==========
print(f"\n[4/5] 筛选三个极端组的位点...")

# 定义三个关键组
key_groups = {
    'M-only': {'color': '#2E86C1', 'label': 'M-only (Membrane-exclusive)'},      # 深蓝
    'C&M&N': {'color': '#E74C3C', 'label': 'C&M&N (Multi-compartment)'},        # 红色
    'C-only': {'color': '#9B59B6', 'label': 'C-only (Cytosol-exclusive)'}       # 紫色
}

sites_by_group = {}
for group in key_groups.keys():
    sites = protein_sites[protein_sites['LocalizationGroup'] == group]
    # 只保留在展示窗口内的位点
    sites = sites[(sites['Position'] >= best_start) & (sites['Position'] <= best_end)]
    sites_by_group[group] = sites
    print(f"   {group}: {len(sites)} 个位点（在展示区域内）")

# 其他组合并为"其他"
other_sites = protein_sites[~protein_sites['LocalizationGroup'].isin(key_groups.keys())]
other_sites = other_sites[(other_sites['Position'] >= best_start) & (other_sites['Position'] <= best_end)]
print(f"   其他组: {len(other_sites)} 个位点（灰化处理）")

# ========== 5. 绘图 ==========
print(f"\n[5/5] 绘制图形...")

fig, ax = plt.subplots(figsize=(16, 7))

# 截取展示区域的残基数据
residue_subset = residue_df[(residue_df['seqpos'] >= best_start) & 
                             (residue_df['seqpos'] <= best_end)]

# 画 3 条关键曲线（粗线，清晰）
ax.plot(residue_subset['seqpos'], residue_subset['backbone'], 
        label='Backbone Dynamics', color='#2C3E50', linewidth=2.5, alpha=0.9, zorder=10)

ax.plot(residue_subset['seqpos'], residue_subset['disoMine'], 
        label='Disorder Propensity', color='#16A085', linewidth=2, alpha=0.8, 
        linestyle='--', zorder=9)

ax.plot(residue_subset['seqpos'], residue_subset['helix'], 
        label='Helix Propensity', color='#D35400', linewidth=2, alpha=0.8, 
        linestyle=':', zorder=8)

# 先画"其他"组（灰色，作为背景）
for _, site in other_sites.iterrows():
    ax.axvline(x=site['Position'], color='#BDC3C7', alpha=0.3, linewidth=1, 
               linestyle=':', zorder=1)

# 再画三个关键组（彩色，清晰）
for group, info in key_groups.items():
    sites = sites_by_group[group]
    for _, site in sites.iterrows():
        ax.axvline(x=site['Position'], color=info['color'], alpha=0.7, 
                   linewidth=2, linestyle='--', zorder=5)

# 图例（分两列：曲线 + 位点组）
handles, labels = ax.get_legend_handles_labels()

from matplotlib.patches import Patch
for group, info in key_groups.items():
    handles.append(Patch(color=info['color'], label=info['label']))

handles.append(Patch(color='#BDC3C7', alpha=0.3, label='Other groups'))

ax.legend(handles=handles, loc='upper right', fontsize=10, ncol=2, 
          framealpha=0.95, edgecolor='black')

# 坐标轴和标题
ax.set_xlabel(f"Sequence Position (aa {best_start}-{best_end})", fontsize=13, fontweight='bold')
ax.set_ylabel("B2B Feature Value", fontsize=13, fontweight='bold')
ax.set_title(f"Biophysical Landscape of {selected_gene} (Region: {best_start}-{best_end} aa)\n" + 
             "Three-Group Comparison: M-only vs C&M&N vs C-only", 
             fontsize=14, fontweight='bold', pad=15)

ax.set_xlim(best_start, best_end)
ax.set_ylim(-0.1, 1.05)
ax.grid(axis='y', alpha=0.3, linestyle='--')

# 添加背景色区分（可选）
ax.axhspan(0.7, 1.05, alpha=0.05, color='red', zorder=0)  # 高值区（浅红底）
ax.text(best_start + 20, 0.88, 'High dynamics/disorder region', 
        fontsize=9, color='gray', style='italic')

plt.tight_layout()

# 保存
figfile_A = f"A_{selected_gene}_FeatureMap_ThreeGroup.png"
plt.savefig(figfile_A, dpi=300, bbox_inches='tight', facecolor='white')
print(f"   ✓ 已保存：{figfile_A}")

plt.show()

# ========== 完成 ==========
print("\n" + "=" * 60)
print("✓ Cell 15A 执行成功！")
print("=" * 60)
print(f"\n生成的图片：{figfile_A}")
print(f"\n关键改进：")
print("  1. ✅ 只展示 {best_start}-{best_end} aa（最密集区域）")
print("  2. ✅ 只画 3 条曲线（Backbone, Disorder, Helix）")
print("  3. ✅ 突出 3 个极端组（M-only蓝, C&M&N红, C-only紫）")
print("  4. ✅ 其他组灰化处理（不干扰视线）")
print("  5. ✅ 优化颜色对比（深色曲线 + 鲜明位点色）")
print("\n明天给博后看时可以说：")
print('  "I selected the most information-dense region and highlighted')
print('   the three biophysically extreme groups based on heatmap analysis."')
print("=" * 60)


# %%
# ===== Cell 20: 保存关键变量（防止明天重跑 Cell 1-19）=====

import pickle
import json

print("=" * 60)
print("保存关键变量到文件")
print("=" * 60)

# 1. 保存 analysisdf_A 和 analysisdf_C（已经有 CSV，但 pickle 更快）
print("\n[1/4] 保存分析数据框...")
analysisdf_A.to_pickle("analysisdf_A.pkl")
analysisdf_C.to_pickle("analysisdf_C.pkl")
print("   ✓ analysisdf_A.pkl")
print("   ✓ analysisdf_C.pkl")

# 2. 保存 gene_to_uniprot 映射（关键！明天画图要用）
print("\n[2/4] 保存基因-Uniprot 映射...")
with open("gene_to_uniprot.json", 'w', encoding='utf-8') as f:
    json.dump(gene_to_uniprot, f, indent=2, ensure_ascii=False)
print("   ✓ gene_to_uniprot.json")

# 3. 保存 JSON_FOLDER 路径（防止忘记）
print("\n[3/4] 保存配置信息...")
config = {
    'JSON_FOLDER': JSON_FOLDER,
    'total_sites': len(analysisdf_A),
    'feature_columns': ['backbone_dynamics', 'disorder_propensity', 
                        'helix_propensity', 'sheet_propensity', 'coil_propensity'],
    'selected_gene': 'SRRM2',  # 如果改了，明天直接改这里
    'date_generated': '2026-02-12'
}
with open("analysis_config.json", 'w', encoding='utf-8') as f:
    json.dump(config, f, indent=2)
print("   ✓ analysis_config.json")

# 4. 保存分组函数（以防 kernel 重启）
print("\n[4/4] 保存分组函数代码...")
functions_code = '''
# 分组函数（从 Cell 10 复制）
def assign_group_A(row):
    c = row['Cyt_any']
    m = row['Mem_any']
    n = row['Nuc_any']
    if c == 1 and m == 0 and n == 0:
        return "C-only"
    elif c == 0 and m == 1 and n == 0:
        return "M-only"
    elif c == 0 and m == 0 and n == 1:
        return "N-only"
    elif c == 1 and m == 1 and n == 0:
        return "C&M"
    elif c == 1 and m == 0 and n == 1:
        return "C&N"
    elif c == 0 and m == 1 and n == 1:
        return "M&N"
    elif c == 1 and m == 1 and n == 1:
        return "C&M&N"
    else:
        return "None"

def get_state_at_time(row, time_label):
    c = row[f'{time_label}_Cyt']
    m = row[f'{time_label}_Mem']
    n = row[f'{time_label}_Nuc']
    if c == 1 and m == 0 and n == 0:
        return "C"
    elif c == 0 and m == 1 and n == 0:
        return "M"
    elif c == 0 and m == 0 and n == 1:
        return "N"
    elif c == 1 and m == 1 and n == 0:
        return "CM"
    elif c == 1 and m == 0 and n == 1:
        return "CN"
    elif c == 0 and m == 1 and n == 1:
        return "MN"
    elif c == 1 and m == 1 and n == 1:
        return "CMN"
    else:
        return "None"
'''

with open("grouping_functions.py", 'w', encoding='utf-8') as f:
    f.write(functions_code)
print("   ✓ grouping_functions.py")

print("\n" + "=" * 60)
print("✓ 所有关键变量已保存！")
print("=" * 60)
print("\n明天可以这样快速加载：")
print("""
import pandas as pd
import json

# 加载数据
analysisdf_A = pd.read_pickle("analysisdf_A.pkl")
analysisdf_C = pd.read_pickle("analysisdf_C.pkl")

# 加载映射
with open("gene_to_uniprot.json", 'r') as f:
    gene_to_uniprot = json.load(f)

# 加载配置
with open("analysis_config.json", 'r') as f:
    config = json.load(f)
    JSON_FOLDER = config['JSON_FOLDER']

# 加载函数
exec(open("grouping_functions.py").read())

print("✓ 所有变量已加载，可以直接运行新代码！")
""")
print("=" * 60)


# %%
# ===== 快速启动：加载所有已保存的数据（跳过 Cell 1-19）=====

import pandas as pd
import json
import os

print("快速加载分析环境...")

# 1. 加载数据
analysisdf_A = pd.read_pickle("analysisdf_A.pkl")
analysisdf_C = pd.read_pickle("analysisdf_C.pkl")
print(f"✓ analysisdf_A: {len(analysisdf_A)} sites")
print(f"✓ analysisdf_C: {len(analysisdf_C)} sites")

# 2. 加载映射
with open("gene_to_uniprot.json", 'r', encoding='utf-8') as f:
    gene_to_uniprot = json.load(f)
print(f"✓ gene_to_uniprot: {len(gene_to_uniprot)} genes")

# 3. 加载配置
with open("analysis_config.json", 'r', encoding='utf-8') as f:
    config = json.load(f)
    JSON_FOLDER = config['JSON_FOLDER']
    feature_cols = config['feature_columns']
print(f"✓ JSON_FOLDER: {JSON_FOLDER}")

# 4. 加载函数
exec(open("grouping_functions.py").read())
print("✓ 分组函数已加载")

print("\n所有数据已就绪！现在可以运行新代码了。")


# %%
# Cell 19C: C方案优化版 SRRM2 Feature Map（适配快速启动版）

import json
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

print("="*60)
print("Cell 19C: C方案 SRRM2 Feature Map (优化版)")
print("="*60)

# ✅ 0. 适配变量名（昨天的快速启动用的是带下划线的）
print("[0] 检查数据...")

try:
    # 优先使用快速启动的变量名（带下划线）
    analysisdfC = analysisdf_C
    genetouniprot = gene_to_uniprot
    JSONFOLDER = JSON_FOLDER
    print(f"   ✓ 使用快速启动数据 (shape: {analysisdfC.shape})")
except NameError:
    # 如果没有运行快速启动，尝试原始变量名
    try:
        test = analysisdfC
        print(f"   ✓ 使用原始变量名 (shape: {analysisdfC.shape})")
    except NameError:
        raise NameError(
            "❌ 错误：找不到数据！\n"
            "   请先运行快速启动 Cell，或运行 Cell 11C 生成数据。"
        )

# 确保 genetouniprot 和 JSONFOLDER 存在
if 'genetouniprot' not in locals():
    genetouniprot = {'SRRM2': 'Q9UQ35'}
    print(f"   ⚠ 使用默认 genetouniprot")

if 'JSONFOLDER' not in locals():
    JSONFOLDER = r'D:\json'
    print(f"   ⚠ 使用默认 JSONFOLDER")

print()

# ✅ 1. 选定基因
selectedgene = 'SRRM2'
print(f"[1] 选定基因: {selectedgene}")

# ✅ 2. 加载 JSON B2B features
print(f"[2] 加载 {selectedgene} B2B features...")

selecteduniprot = genetouniprot[selectedgene]
jsonpath = os.path.join(JSONFOLDER, f'{selecteduniprot}.json')

if not os.path.exists(jsonpath):
    raise FileNotFoundError(f"❌ JSON 文件不存在: {jsonpath}")

with open(jsonpath, 'r', encoding='utf-8') as f:
    proteindata = json.load(f)
    residues = proteindata['residues']

residuedf = pd.DataFrame(residues)

proteinlength = residuedf['seqpos'].max()
print(f"   蛋白总长: {proteinlength} aa")
print(f"   范围: {residuedf['seqpos'].min()} - {proteinlength}")

# ✅ 3. 提取 SRRM2 的 C 图位点（sustained-single & sustained-multi）
print(f"[3] 提取 {selectedgene} C 图位点...")

proteinsitesC = analysisdfC[analysisdfC['Gene'] == selectedgene].copy()

# 只保留 sustained-single 和 sustained-multi
validgroups = ['sustained-single', 'sustained-multi']
proteinsitesC = proteinsitesC[proteinsitesC['DynamicGroup'].isin(validgroups)]

print(f"   总位点数: {len(proteinsitesC)}")
print(f"   C图分组统计:")
print(proteinsitesC['DynamicGroup'].value_counts())

# ✅ 4. 滑动窗口寻找最密集的 600 aa 区域
print(f"[4] 滑动窗口寻找最密集的 600 aa 区域...")

windowsize = 600
beststart = 0
maxsitesinwindow = 0

for start in range(0, proteinlength - windowsize, 100):
    end = start + windowsize
    sitesinwindow = len(proteinsitesC[
        (proteinsitesC['Position'] >= start) & 
        (proteinsitesC['Position'] <= end)
    ])
    
    if sitesinwindow > maxsitesinwindow:
        maxsitesinwindow = sitesinwindow
        beststart = start

bestend = beststart + windowsize

print(f"   最佳窗口: {beststart}-{bestend} aa")
print(f"   该区域位点数: {maxsitesinwindow} ({maxsitesinwindow/len(proteinsitesC)*100:.1f}%)")

# ✅ 5. 按组划分位点
print(f"[5] 按组划分位点...")

# C图只有两个组（sustained-single 和 sustained-multi）
keygroups = {
    'sustained-single': {'color': '#2C3E50', 'label': 'Sustained-Single (Stable)'},
    'sustained-multi': {'color': '#E74C3C', 'label': 'Sustained-Multi (Dynamic)'}
}

sitesbygroup = {}

for group in validgroups:
    sites = proteinsitesC[proteinsitesC['DynamicGroup'] == group]
    # 筛选在窗口内的位点
    sites = sites[(sites['Position'] >= beststart) & (sites['Position'] <= bestend)]
    sitesbygroup[group] = sites
    print(f"   {group}: {len(sites)} 个位点")

# ✅ 6. 绘制 feature map
print(f"[6] 绘制 feature map...")

fig, ax = plt.subplots(figsize=(16, 7))

# 提取窗口内的 residue dataframe
residuesubset = residuedf[
    (residuedf['seqpos'] >= beststart) & 
    (residuedf['seqpos'] <= bestend)
]

# 绘制 3 条 biophysical curves
ax.plot(residuesubset['seqpos'], residuesubset['backbone'], 
        label='Backbone Dynamics', color='#2C3E50', linewidth=2.5, alpha=0.9, zorder=10)

ax.plot(residuesubset['seqpos'], residuesubset['disoMine'], 
        label='Disorder Propensity', color='#16A085', linewidth=2, alpha=0.8, 
        linestyle='--', zorder=9)

ax.plot(residuesubset['seqpos'], residuesubset['helix'], 
        label='Helix Propensity', color='#D35400', linewidth=2, alpha=0.8, 
        linestyle=':', zorder=8)

# 绘制磷酸化位点的竖线（分组标注）
for group, info in keygroups.items():
    sites = sitesbygroup[group]
    for _, site in sites.iterrows():
        ax.axvline(x=site['Position'], color=info['color'], alpha=0.7, 
                  linewidth=2, linestyle='--', zorder=5)

# 图例
handles, labels = ax.get_legend_handles_labels()

from matplotlib.patches import Patch
for group, info in keygroups.items():
    handles.append(Patch(color=info['color'], label=info['label']))

ax.legend(handles=handles, loc='upper right', fontsize=10, ncol=2, 
         framealpha=0.95, edgecolor='black')

# 坐标轴
ax.set_xlabel(f'Sequence Position (aa {beststart}-{bestend})', fontsize=13, fontweight='bold')
ax.set_ylabel('B2B Feature Value', fontsize=13, fontweight='bold')
ax.set_title(f'Biophysical Landscape of {selectedgene} (Region {beststart}-{bestend} aa)\nTwo-Group Comparison: Sustained-Single vs Sustained-Multi', 
            fontsize=14, fontweight='bold', pad=15)

ax.set_xlim(beststart, bestend)
ax.set_ylim(-0.1, 1.05)
ax.grid(axis='y', alpha=0.3, linestyle='--')

# 高亮高动态/无序区域
ax.axhspan(0.7, 1.05, alpha=0.05, color='red', zorder=0)
ax.text(beststart + 20, 0.88, 'High dynamics/disorder region', 
       fontsize=9, color='gray', style='italic')

plt.tight_layout()

# 保存
figfileC = f'C_{selectedgene}_FeatureMapTwoGroup.png'
plt.savefig(figfileC, dpi=300, bbox_inches='tight', facecolor='white')
print(f"   ✅ 保存图片: {figfileC}")

plt.show()

print("="*60)
print("✅ Cell 19C 完成")
print("="*60)
print()
print("输出:")
print(f"  1. {figfileC}")
print()
print("说明:")
print("  我筛选了最信息密集的区域，并突出显示了")
print("  两个组（Sustained-Single vs Sustained-Multi）的差异。")


# %%
#快速加载数据
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("快速加载/生成分析数据...")
print("="*60)

# 定义路径
BASE_PATH = r"D:\博士\Phospho"
success_file = os.path.join(BASE_PATH, "success_df_M1M2_cleaned.csv")
groupA_file = os.path.join(BASE_PATH, "A_Grouped_Sites.csv")
groupC_file = os.path.join(BASE_PATH, "C_Grouped_Sites.csv")

# 定义特征列
feature_cols = ['backbone_dynamics', 'disorder_propensity', 'helix_propensity', 
                'sheet_propensity', 'coil_propensity']

# ===== 1. 加载主数据 =====
print("\n【1/3】加载主数据...")
if os.path.exists(success_file):
    success_df = pd.read_csv(success_file)
    print(f"  ✓ 成功加载: {len(success_df):,} 行")
else:
    print(f"  ✗ 找不到文件: {success_file}")
    raise FileNotFoundError("请先运行生成success_df的代码！")

# ===== 2. 加载或生成方案A数据 =====
print("\n【2/3】处理方案A数据...")
if os.path.exists(groupA_file):
    analysisdf_A = pd.read_csv(groupA_file)
    print(f"  ✓ 加载已有文件: {len(analysisdf_A):,} 行")
else:
    print("  ⚠ 文件不存在，正在生成...")
    
    # 定义分组函数
    def assign_group_A(row):
        c, m, n = row['Cytany'], row['Memany'], row['Nucany']
        if c == 1 and m == 0 and n == 0: return 'C-only'
        elif c == 0 and m == 1 and n == 0: return 'M-only'
        elif c == 0 and m == 0 and n == 1: return 'N-only'
        elif c == 1 and m == 1 and n == 0: return 'CM'
        elif c == 1 and m == 0 and n == 1: return 'CN'
        elif c == 0 and m == 1 and n == 1: return 'MN'
        elif c == 1 and m == 1 and n == 1: return 'CMN'
        else: return None
    
    # 计算time-collapsed any
    success_df['Cytany'] = success_df[['0min_Cyt', '2min_Cyt', '8min_Cyt', '20min_Cyt', '90min_Cyt']].max(axis=1)
    success_df['Memany'] = success_df[['0min_Mem', '2min_Mem', '8min_Mem', '20min_Mem', '90min_Mem']].max(axis=1)
    success_df['Nucany'] = success_df[['0min_Nuc', '2min_Nuc', '8min_Nuc', '20min_Nuc', '90min_Nuc']].max(axis=1)
    
    # 分配分组
    success_df['LocalizationGroup'] = success_df.apply(assign_group_A, axis=1)
    
    # 过滤掉None
    analysisdf_A = success_df[success_df['LocalizationGroup'].notna()].copy()
    
    # 保存
    analysisdf_A.to_csv(groupA_file, index=False, encoding='utf-8-sig')
    print(f"  ✓ 生成并保存: {len(analysisdf_A):,} 行")
    print(f"  分组分布:\n{analysisdf_A['LocalizationGroup'].value_counts()}")

# ===== 3. 加载或生成方案C数据 =====
print("\n【3/3】处理方案C数据...")
if os.path.exists(groupC_file):
    analysisdf_C = pd.read_csv(groupC_file)
    print(f"  ✓ 加载已有文件: {len(analysisdf_C):,} 行")
else:
    print("  ⚠ 文件不存在，正在生成...")
    
    # 定义状态获取函数
    def get_state_at_time(row, timelabel):
        c = row[f'{timelabel}_Cyt']
        m = row[f'{timelabel}_Mem']
        n = row[f'{timelabel}_Nuc']
        if c == 1 and m == 0 and n == 0: return 'C'
        elif c == 0 and m == 1 and n == 0: return 'M'
        elif c == 0 and m == 0 and n == 1: return 'N'
        elif c == 1 and m == 1 and n == 0: return 'CM'
        elif c == 1 and m == 0 and n == 1: return 'CN'
        elif c == 0 and m == 1 and n == 1: return 'MN'
        elif c == 1 and m == 1 and n == 1: return 'CMN'
        else: return None
    
    # 定义分组函数C
    def assign_group_C(row):
        traj = row['trajectory']
        states = [s for s in traj if s != 'None']
        nonnone_states = [s for s in states if s is not None]
        
        if len(nonnone_states) == 0:
            return None
        
        if len(set(states)) == 1 and states[0] in ['C', 'M', 'N']:
            return 'sustained-single'
        
        if len(set(states)) == 1 and states[0] in ['CN', 'CM', 'MN', 'CMN']:
            return 'sustained-multi'
        
        if len(set(nonnone_states)) > 1:
            return 'translocation'
        
        return 'other'
    
    # 生成5个时间点的状态
    time_labels = ['0min', '2min', '8min', '20min', '90min']
    for t in time_labels:
        success_df[f'state_{t}'] = success_df.apply(lambda row: get_state_at_time(row, t), axis=1)
    
    # 生成trajectory字符串
    success_df['trajectory'] = success_df.apply(
        lambda row: '-'.join([str(row[f'state_{t}']) for t in time_labels]), axis=1
    )
    
    # 分配分组
    success_df['DynamicGroup'] = success_df.apply(assign_group_C, axis=1)
    
    # 过滤掉None和other
    analysisdf_C = success_df[success_df['DynamicGroup'].notna()].copy()
    analysisdf_C = analysisdf_C[analysisdf_C['DynamicGroup'] != 'other']
    
    # 保存
    analysisdf_C.to_csv(groupC_file, index=False, encoding='utf-8-sig')
    print(f"  ✓ 生成并保存: {len(analysisdf_C):,} 行")
    print(f"  分组分布:\n{analysisdf_C['DynamicGroup'].value_counts()}")

print("\n" + "="*60)
print("✅ 所有数据准备完毕！")
print("="*60)
print(f"\n当前变量:")
print(f"  • success_df: {len(success_df):,} 行（原始数据）")
print(f"  • analysisdf_A: {len(analysisdf_A):,} 行（方案A）")
print(f"  • analysisdf_C: {len(analysisdf_C):,} 行（方案C）")
print(f"  • feature_cols: {feature_cols}")


# %%
print("="*70)
print("【阶段1】蛋白质子集筛选 - 修正版")
print("="*70)

# ===== 1. 方案A候选蛋白筛选 =====
print("\n【1/4】方案A筛选...")
print("  筛选标准：")
print("    • 位点数: 15-30个")
print("    • 分组覆盖: ≥5个LocalizationGroup")
print("    • 排除SRRM2（位点数过多，不具代表性）")

# 按基因分组统计
proteinA_stats = analysisdf_A.groupby('Gene').agg(
    num_sites=('Position', 'count'),
    num_groups=('LocalizationGroup', 'nunique')
).reset_index()

# 应用筛选条件
candidateA = proteinA_stats[
    (proteinA_stats['num_sites'] >= 15) & 
    (proteinA_stats['num_sites'] <= 30) &
    (proteinA_stats['num_groups'] >= 5) &
    (proteinA_stats['Gene'] != 'SRRM2')
].sort_values(['num_groups', 'num_sites'], ascending=False)

print(f"\n  结果: {len(candidateA)} 个候选蛋白")
print(f"\n  Top 15:")
print(candidateA.head(15).to_string(index=False))

# ===== 2. 方案C候选蛋白筛选 =====
print("\n\n【2/4】方案C筛选...")
print("  筛选标准：")
print("    • 位点数: 15-30个")
print("    • 必须同时有sustained-single和sustained-multi")
print("    • 两组位点数量不能太不平衡（各占≥20%）")

# 先检查DynamicGroup列是否存在
if 'DynamicGroup' not in analysisdf_C.columns:
    print("\n  ✗ 错误: analysisdf_C中没有'DynamicGroup'列")
    print(f"  现有列: {analysisdf_C.columns.tolist()}")
    candidateC = pd.DataFrame(columns=['Gene', 'num_sites', 'min_group_ratio'])
else:
    # 计算每个蛋白的分组平衡度
    balance_check = []
    
    for gene in analysisdf_C['Gene'].unique():
        gene_df = analysisdf_C[analysisdf_C['Gene'] == gene]
        group_counts = gene_df['DynamicGroup'].value_counts()
        
        # 检查是否有两个组且位点数在范围内
        if len(group_counts) == 2 and 15 <= len(gene_df) <= 30:
            min_ratio = group_counts.min() / len(gene_df)
            
            # 检查平衡度
            if min_ratio >= 0.2:
                balance_check.append({
                    'Gene': gene,
                    'num_sites': len(gene_df),
                    'min_group_ratio': min_ratio,
                    'group1': group_counts.index[0],
                    'group1_count': group_counts.iloc[0],
                    'group2': group_counts.index[1],
                    'group2_count': group_counts.iloc[1]
                })
    
    # 创建DataFrame
    if len(balance_check) > 0:
        candidateC = pd.DataFrame(balance_check).sort_values('num_sites', ascending=False)
        print(f"\n  结果: {len(candidateC)} 个候选蛋白")
        print(f"\n  Top 15:")
        display_cols = ['Gene', 'num_sites', 'min_group_ratio', 'group1_count', 'group2_count']
        print(candidateC[display_cols].head(15).to_string(index=False))
    else:
        print(f"\n  ⚠ 结果: 0 个候选蛋白（标准太严格，尝试放宽...）")
        
        # 放宽条件：去掉平衡度限制
        balance_check_relaxed = []
        for gene in analysisdf_C['Gene'].unique():
            gene_df = analysisdf_C[analysisdf_C['Gene'] == gene]
            group_counts = gene_df['DynamicGroup'].value_counts()
            
            if len(group_counts) == 2 and 15 <= len(gene_df) <= 30:
                min_ratio = group_counts.min() / len(gene_df)
                balance_check_relaxed.append({
                    'Gene': gene,
                    'num_sites': len(gene_df),
                    'min_group_ratio': min_ratio
                })
        
        if len(balance_check_relaxed) > 0:
            candidateC = pd.DataFrame(balance_check_relaxed).sort_values('num_sites', ascending=False)
            print(f"\n  放宽后结果: {len(candidateC)} 个候选蛋白")
            print(f"\n  Top 15:")
            print(candidateC.head(15).to_string(index=False))
        else:
            print(f"\n  即使放宽条件，仍然0个候选蛋白")
            candidateC = pd.DataFrame(columns=['Gene', 'num_sites', 'min_group_ratio'])

# ===== 3. 找到A和C的交集 =====
print("\n\n【3/4】A+C交集分析...")

if len(candidateC) > 0:
    common_genes = set(candidateA['Gene']) & set(candidateC['Gene'])
    print(f"\n  共同候选: {len(common_genes)} 个蛋白")
else:
    common_genes = set()
    print(f"\n  ⚠ 方案C无候选蛋白，无法计算交集")

# ===== 4. 最终选择策略 =====
print("\n\n【4/4】最终子集选择...")

if len(common_genes) >= 8:
    # 策略1: 交集充足，直接使用
    print(f"  策略: 使用A+C交集（{len(common_genes)}个蛋白）")
    selected_genes = list(common_genes)[:12]  # 最多12个
    
elif len(common_genes) > 0:
    # 策略2: 交集不足，从方案A补充
    print(f"  策略: 交集{len(common_genes)}个，从方案A补充")
    selected_genes = list(common_genes)
    remaining = candidateA[~candidateA['Gene'].isin(common_genes)].head(12 - len(common_genes))
    selected_genes.extend(remaining['Gene'].tolist())
    
else:
    # 策略3: 无交集，仅使用方案A
    print(f"  策略: 无交集，仅使用方案A前12个蛋白")
    selected_genes = candidateA.head(12)['Gene'].tolist()

# 生成子集
subset_dfA = analysisdf_A[analysisdf_A['Gene'].isin(selected_genes)].copy()
subset_dfC = analysisdf_C[analysisdf_C['Gene'].isin(selected_genes)].copy()

print(f"\n  最终选择: {len(selected_genes)} 个蛋白")
print(f"    • 方案A子集: {len(subset_dfA):,} 个位点")
print(f"    • 方案C子集: {len(subset_dfC):,} 个位点")

print(f"\n  入选蛋白详情:")
for i, gene in enumerate(selected_genes, 1):
    a_count = len(subset_dfA[subset_dfA['Gene'] == gene])
    c_count = len(subset_dfC[subset_dfC['Gene'] == gene])
    a_groups = subset_dfA[subset_dfA['Gene'] == gene]['LocalizationGroup'].nunique()
    print(f"    {i:2d}. {gene:15s} - A:{a_count:2d}位点/{a_groups}组, C:{c_count:2d}位点")

# 保存
subset_info_file = os.path.join(BASE_PATH, "Selected_Protein_Subset.csv")
pd.DataFrame({
    'Gene': selected_genes,
    'A_sites': [len(subset_dfA[subset_dfA['Gene'] == g]) for g in selected_genes],
    'C_sites': [len(subset_dfC[subset_dfC['Gene'] == g]) for g in selected_genes]
}).to_csv(subset_info_file, index=False, encoding='utf-8-sig')

print(f"\n  ✓ 已保存: {subset_info_file}")

print("\n" + "="*70)
print("✅ 阶段1完成！")
print("="*70)


# %%
print("="*70)
print("【方案C筛选思路对比 - 数据统计】")
print("="*70)

# ===== 前置：检查数据 =====
print("\n【0】数据概况...")
print(f"  analysisdf_C总位点数: {len(analysisdf_C):,}")
print(f"  涉及蛋白数: {analysisdf_C['Gene'].nunique():,}")
print(f"  DynamicGroup分布:")
print(analysisdf_C['DynamicGroup'].value_counts())

# ===== 思路1: 蛋白内部异质性（既有single又有multi）=====
print("\n\n【思路1】蛋白内部异质性分析")
print("  标准: 同一蛋白既有sustained-single位点，又有sustained-multi位点")
print("-"*70)

thought1_candidates = []

for gene in analysisdf_C['Gene'].unique():
    gene_df = analysisdf_C[analysisdf_C['Gene'] == gene]
    
    # 统计该蛋白的两组位点数
    single_count = len(gene_df[gene_df['DynamicGroup'] == 'sustained-single'])
    multi_count = len(gene_df[gene_df['DynamicGroup'] == 'sustained-multi'])
    total_count = len(gene_df)
    
    # 筛选条件：两组都有，且总数≥15
    if single_count > 0 and multi_count > 0 and total_count >= 15:
        thought1_candidates.append({
            'Gene': gene,
            'total_sites': total_count,
            'single_sites': single_count,
            'multi_sites': multi_count,
            'single_ratio': single_count / total_count,
            'multi_ratio': multi_count / total_count
        })

# 安全创建DataFrame
if len(thought1_candidates) > 0:
    thought1_df = pd.DataFrame(thought1_candidates).sort_values('total_sites', ascending=False)
    print(f"\n  结果: {len(thought1_df)} 个候选蛋白")
    print(f"\n  详细列表 (前20):")
    print(thought1_df.head(20).to_string(index=False))
    
    # 额外统计：平衡度分布
    balanced = thought1_df[
        (thought1_df['single_ratio'] >= 0.2) & 
        (thought1_df['multi_ratio'] >= 0.2)
    ]
    print(f"\n  其中平衡蛋白（两组各≥20%）: {len(balanced)} 个")
else:
    thought1_df = pd.DataFrame()
    print(f"\n  结果: 0 个候选蛋白")
    print("  ✗ 没有蛋白同时包含sustained-single和sustained-multi位点")

# ===== 思路2: 位点数≥15 + 有两组 =====
print("\n\n【思路2】简单位点数筛选")
print("  标准: 位点数≥15，且分布在sustained-single和sustained-multi两组")
print("-"*70)

thought2_candidates = []

for gene in analysisdf_C['Gene'].unique():
    gene_df = analysisdf_C[analysisdf_C['Gene'] == gene]
    group_counts = gene_df['DynamicGroup'].value_counts()
    
    # 筛选条件：≥15个位点，且有2个组
    if len(gene_df) >= 15 and len(group_counts) == 2:
        thought2_candidates.append({
            'Gene': gene,
            'num_sites': len(gene_df),
            'group1': group_counts.index[0],
            'group1_count': group_counts.iloc[0],
            'group2': group_counts.index[1],
            'group2_count': group_counts.iloc[1]
        })

# 安全创建DataFrame
if len(thought2_candidates) > 0:
    thought2_df = pd.DataFrame(thought2_candidates).sort_values('num_sites', ascending=False)
    print(f"\n  结果: {len(thought2_df)} 个候选蛋白")
    print(f"\n  详细列表 (前20):")
    print(thought2_df.head(20).to_string(index=False))
else:
    thought2_df = pd.DataFrame()
    print(f"\n  结果: 0 个候选蛋白")
    print("  ✗ 没有蛋白满足≥15位点且有两组")

# ===== 思路3: 仅用方案A，方案C跟随 =====
print("\n\n【思路3】仅用方案A筛选")
print("  标准: 使用方案A的26个候选蛋白，检查它们在方案C中的覆盖")
print("-"*70)

# 使用之前筛选出的candidateA
thought3_genes = candidateA['Gene'].tolist()

# 检查这些蛋白在方案C中的情况
thought3_coverage = []
for gene in thought3_genes:
    c_df = analysisdf_C[analysisdf_C['Gene'] == gene]
    
    if len(c_df) > 0:
        thought3_coverage.append({
            'Gene': gene,
            'A_sites': len(analysisdf_A[analysisdf_A['Gene'] == gene]),
            'C_sites': len(c_df),
            'C_groups': c_df['DynamicGroup'].nunique(),
            'has_both_C_groups': c_df['DynamicGroup'].nunique() == 2
        })

# 安全创建DataFrame
if len(thought3_coverage) > 0:
    thought3_df = pd.DataFrame(thought3_coverage)
    print(f"\n  方案A的26个蛋白中:")
    print(f"    • 在方案C中有数据: {len(thought3_df)} 个")
    print(f"    • 在方案C有两组: {thought3_df['has_both_C_groups'].sum()} 个")
    print(f"\n  详细列表:")
    print(thought3_df.to_string(index=False))
else:
    thought3_df = pd.DataFrame()
    print(f"\n  方案A的26个蛋白在方案C中都没有数据！")

# =====


# %%
# 检查SRRM2在方案C的实际分布
srrm2_c = analysisdf_C[analysisdf_C['Gene'] == 'SRRM2']
print(f"SRRM2在方案C中:")
print(f"  总位点数: {len(srrm2_c)}")
print(f"  DynamicGroup分布:")
print(srrm2_c['DynamicGroup'].value_counts())


# %%
print("="*70)
print("【修复assigngroupC并重新分组】")
print("="*70)

# 1. 定义正确的分组函数
def assigngroupC_fixed(row):
    """修复版：正确处理trajectory字符串"""
    traj = row['trajectory']
    
    # 正确split
    states = traj.split('-')
    nonnone_states = [s for s in states if s not in ['None', '', 'nan']]
    
    if len(nonnone_states) == 0:
        return None
    
    # sustained-single: CCCCC, NNNNN, MMMMM
    if len(set(nonnone_states)) == 1 and nonnone_states[0] in ['C', 'M', 'N']:
        return 'sustained-single'
    
    # sustained-multi: CMN-CMN-CMN, MN-MN-MN等
    if len(set(nonnone_states)) == 1 and nonnone_states[0] in ['CN', 'CM', 'MN', 'CMN']:
        return 'sustained-multi'
    
    # early-only: 前期出现后消失
    if len(nonnone_states) <= 2 and all(s == 'None' for s in states[2:]):
        return 'early-only'
    
    # late-appearing: 后期才出现
    if all(s == 'None' for s in states[:2]) and len(nonnone_states) >= 2:
        return 'late-appearing'
    
    # translocation: 真正的状态改变
    if len(set(nonnone_states)) > 1:
        return 'translocation'
    
    return 'other'

# 2. 重新应用分组
print("\n重新分组中...")
success_df['DynamicGroup_fixed'] = success_df.apply(assigngroupC_fixed, axis=1)

# 3. 对比新旧分组
print("\n修复前（错误）:")
print(success_df['DynamicGroup'].value_counts())

print("\n修复后（正确）:")
print(success_df['DynamicGroup_fixed'].value_counts())

# 4. 检查SRRM2
srrm2_fixed = success_df[success_df['Gene'] == 'SRRM2']
print(f"\nSRRM2修复后的分组:")
print(srrm2_fixed['DynamicGroup_fixed'].value_counts())

# 5. 生成新的analysisdf_C
analysisdf_C_fixed = success_df[
    success_df['DynamicGroup_fixed'].notna() & 
    (success_df['DynamicGroup_fixed'] != 'other')
].copy()

print(f"\n新的方案C数据:")
print(f"  总位点数: {len(analysisdf_C_fixed):,}")
print(f"  分组分布:")
print(analysisdf_C_fixed['DynamicGroup_fixed'].value_counts())

# 6. 保存修复后的数据
analysisdf_C_fixed = analysisdf_C_fixed.rename(columns={'DynamicGroup_fixed': 'DynamicGroup'})
analysisdf_C_fixed.to_csv(
    os.path.join(BASE_PATH, "C_Grouped_Sites_FIXED.csv"),
    index=False,
    encoding='utf-8-sig'
)

print("\n✓ 已保存: C_Grouped_Sites_FIXED.csv")
print("="*70)


# %%
print("="*70)
print("【细致的trajectory模式分析】")
print("="*70)

# 1. 检查包含None的情况
has_none = success_df['trajectory'].str.contains('None', na=False)
print(f"\n包含None的位点: {has_none.sum():,} / {len(success_df):,} ({has_none.sum()/len(success_df)*100:.1f}%)")

# 2. 找出那个translocation
translocation_site = success_df[success_df['DynamicGroup_fixed'] == 'translocation']
if len(translocation_site) > 0:
    print(f"\n唯一的translocation位点:")
    for idx, row in translocation_site.iterrows():
        print(f"  {row['Gene']} {row['Site']}: {row['trajectory']}")

# 3. 抽样看sustained的trajectory
print(f"\nsustained-single的trajectory样本（前10个）:")
single_samples = success_df[success_df['DynamicGroup_fixed'] == 'sustained-single']['trajectory'].head(10)
for traj in single_samples:
    print(f"  {traj}")

print(f"\nsustained-multi的trajectory样本（前10个）:")
multi_samples = success_df[success_df['DynamicGroup_fixed'] == 'sustained-multi']['trajectory'].head(10)
for traj in multi_samples:
    print(f"  {traj}")

# 4. 统计unique trajectory的数量
print(f"\n全局trajectory的diversity:")
print(f"  Unique trajectories: {success_df['trajectory'].nunique()}")
print(f"  Top 20 most common:")
print(success_df['trajectory'].value_counts().head(20))

print("="*70)


# %%
# 加载之前的C方案结果
import pandas as pd

# 检查是否存在统计结果文件
import os

files_to_check = [
    r"D:\博士\Phospho\C_StatisticalTests.csv",
    r"D:\博士\Phospho\C_GroupStatistics.csv"
]

print("="*70)
print("【检查原C方案的统计结果】")
print("="*70)

for f in files_to_check:
    if os.path.exists(f):
        print(f"\n✓ 找到文件: {f}")
        df = pd.read_csv(f)
        print(f"\n前几行:")
        print(df.head())
    else:
        print(f"\n✗ 未找到: {f}")

# 如果文件不存在，需要重新加载数据检查
if not any(os.path.exists(f) for f in files_to_check):
    print("\n需要重新计算统计检验...")
    
    # 加载C方案分组数据
    try:
        df = pd.read_csv(r"D:\博士\Phospho\success_df_with_dynamics.csv")
        
        # 筛选有效数据
        valid_df = df[df['DynamicGroup'].isin(['sustained-single', 'sustained-multi'])].copy()
        
        print(f"\n原C方案分组:")
        print(valid_df['DynamicGroup'].value_counts())
        
        # 快速t-test
        from scipy import stats
        
        features = ['backbone_dynamics', 'disorder_propensity', 'helix_propensity', 
                   'sheet_propensity', 'sidechain_dynamics']
        
        print(f"\n快速统计检验:")
        print(f"\n{'Feature':<25} {'Single mean':<12} {'Multi mean':<12} {'p-value':<12} {'Significant'}")
        print("-"*80)
        
        for feat in features:
            single_data = valid_df[valid_df['DynamicGroup']=='sustained-single'][feat].dropna()
            multi_data = valid_df[valid_df['DynamicGroup']=='sustained-multi'][feat].dropna()
            
            stat, pval = stats.ttest_ind(single_data, multi_data)
            sig = "✓ YES" if pval < 0.05 else "✗ NO"
            
            print(f"{feat:<25} {single_data.mean():>11.4f} {multi_data.mean():>11.4f} {pval:>11.2e} {sig}")
            
    except Exception as e:
        print(f"\n出错: {e}")
        print("\n可能需要先加载正确的数据文件")

print("\n"+"="*70)


# %%
print("="*70)
print("【阶段2】按氨基酸类型（S/T/Y）分层分析")
print("="*70)

# ===== 1. 检查氨基酸分布 =====
print("\n【1/4】氨基酸分布统计...")

# 全数据集
print("\n全数据集 (10,257位点):")
aa_counts_all = analysisdf_A['AA'].value_counts()
print(aa_counts_all)
print(f"\n  S: {aa_counts_all.get('S', 0):,} ({aa_counts_all.get('S', 0)/len(analysisdf_A)*100:.1f}%)")
print(f"  T: {aa_counts_all.get('T', 0):,} ({aa_counts_all.get('T', 0)/len(analysisdf_A)*100:.1f}%)")
print(f"  Y: {aa_counts_all.get('Y', 0):,} ({aa_counts_all.get('Y', 0)/len(analysisdf_A)*100:.1f}%)")

# 26个候选蛋白
subset_dfA_26 = analysisdf_A[analysisdf_A['Gene'].isin(candidateA.head(26)['Gene'])].copy()
print(f"\n26个候选蛋白 ({len(subset_dfA_26)}位点):")
aa_counts_subset = subset_dfA_26['AA'].value_counts()
print(aa_counts_subset)

# ===== 2. 按S/T/Y分层统计 =====
print("\n\n【2/4】按氨基酸类型的LocalizationGroup分布...")

for aa_type in ['S', 'T', 'Y']:
    aa_df = subset_dfA_26[subset_dfA_26['AA'] == aa_type]
    print(f"\n{aa_type} ({len(aa_df)}位点):")
    if len(aa_df) > 0:
        group_dist = aa_df['LocalizationGroup'].value_counts()
        for group, count in group_dist.items():
            print(f"  {group:15s}: {count:3d} ({count/len(aa_df)*100:.1f}%)")

# ===== 3. 生成分层箱线图 =====
print("\n\n【3/4】生成按氨基酸分层的箱线图...")

import matplotlib.pyplot as plt
import seaborn as sns

# 特征列
feature_cols = ['backbone_dynamics', 'disorder_propensity', 'helix_propensity', 
                'sheet_propensity', 'coil_propensity']

# LocalizationGroup顺序
group_order = ['C-only', 'N-only', 'M-only', 'CN', 'CM', 'MN', 'CMN']

# 为每个氨基酸类型生成图
saved_files = []

for aa_type in ['S', 'T', 'Y']:
    aa_subset = subset_dfA_26[subset_dfA_26['AA'] == aa_type]
    
    print(f"\n  处理 {aa_type} ({len(aa_subset)}位点)...")
    
    if len(aa_subset) < 10:  # 数据太少跳过
        print(f"    ⚠ 数据不足，跳过")
        continue
    
    # 检查有哪些分组
    available_groups = [g for g in group_order if g in aa_subset['LocalizationGroup'].unique()]
    
    if len(available_groups) < 3:  # 分组太少跳过
        print(f"    ⚠ 分组数<3，跳过")
        continue
    
    print(f"    ✓ 有效分组: {available_groups}")
    
    # 为每个特征生成箱线图
    for feature in feature_cols:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 绘制箱线图
        sns.boxplot(
            data=aa_subset,
            x='LocalizationGroup',
            y=feature,
            order=available_groups,
            palette='Set2',
            ax=ax
        )
        
        # 标题和标签
        ax.set_title(f'{aa_type}-phosphorylation: {feature.replace("_", " ").title()} by Localization',
                    fontsize=14, fontweight='bold')
        ax.set_xlabel('Localization Group', fontsize=12)
        ax.set_ylabel(feature.replace('_', ' ').title(), fontsize=12)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', alpha=0.3)
        
        # 添加样本量标注
        for i, group in enumerate(available_groups):
            n = len(aa_subset[aa_subset['LocalizationGroup'] == group])
            ax.text(i, ax.get_ylim()[1]*0.95, f'n={n}', 
                   ha='center', va='top', fontsize=8, color='gray')
        
        plt.tight_layout()
        
        # 保存
        figfile = os.path.join(BASE_PATH, f'A_{aa_type}_{feature}_boxplot.png')
        plt.savefig(figfile, dpi=300, bbox_inches='tight', facecolor='white')
        saved_files.append(figfile)
        plt.close()
    
    print(f"    ✓ 已生成5张图")

print(f"\n  共生成 {len(saved_files)} 张图")

# ===== 4. 生成统计报告 =====
print("\n\n【4/4】统计检验...")

from scipy import stats

# 为每个氨基酸类型做ANOVA
stat_results = []

for aa_type in ['S', 'T', 'Y']:
    aa_subset = subset_dfA_26[subset_dfA_26['AA'] == aa_type]
    
    if len(aa_subset) < 10:
        continue
    
    available_groups = [g for g in group_order if g in aa_subset['LocalizationGroup'].unique()]
    
    if len(available_groups) < 3:
        continue
    
    print(f"\n{aa_type}-phosphorylation:")
    print("-" * 50)
    
    for feature in feature_cols:
        # 收集各组数据
        groups = []
        for group in available_groups:
            data = aa_subset[aa_subset['LocalizationGroup'] == group][feature].dropna()
            if len(data) >= 3:  # 至少3个样本
                groups.append(data)
        
        if len(groups) < 3:
            continue
        
        # ANOVA
        f_stat, p_anova = stats.f_oneway(*groups)
        
        # 判断显著性
        sig = '***' if p_anova < 0.001 else '**' if p_anova < 0.01 else '*' if p_anova < 0.05 else 'ns'
        
        print(f"  {feature:25s}: F={f_stat:6.2f}, p={p_anova:.2e} {sig}")
        
        stat_results.append({
            'AA': aa_type,
            'Feature': feature,
            'F_statistic': f_stat,
            'p_value': p_anova,
            'Significant': p_anova < 0.05
        })

# 保存统计结果
stat_df = pd.DataFrame(stat_results)
stat_file = os.path.join(BASE_PATH, "A_Stratified_by_AA_Statistics.csv")
stat_df.to_csv(stat_file, index=False, encoding='utf-8-sig')

print(f"\n✓ 统计结果已保存: {stat_file}")

# ===== 总结 =====
print("\n" + "="*70)
print("✅ 阶段2完成！")
print("="*70)
print(f"\n生成文件:")
print(f"  • {len(saved_files)} 张箱线图")
print(f"  • 1 个统计报告CSV")

print(f"\n关键发现:")
for aa_type in ['S', 'T', 'Y']:
    aa_results = [r for r in stat_results if r['AA'] == aa_type]
    if len(aa_results) > 0:
        sig_count = sum(1 for r in aa_results if r['Significant'])
        print(f"  {aa_type}: {sig_count}/{len(aa_results)} 特征显著 (p<0.05)")

print("\n下一步: 阶段3 - ±5序列上下文分析")
print("="*70)


# %%
# ============================================================
# Cell A: 定义提取±5氨基酸序列上下文特征的函数
# ============================================================

import pandas as pd
import json
import os

print("="*60)
print("Cell A: 定义序列上下文特征提取函数")
print("="*60)

def extract_context_features(gene_name, position, window=5, 
                            gene_to_uniprot_dict=None, json_folder_path=None):
    """
    提取磷酸化位点前后window个氨基酸的生物物理特征
    
    参数:
        gene_name: 基因名
        position: 中心位点位置（1-based）
        window: 窗口大小（默认5，即提取±5共11个位置）
        gene_to_uniprot_dict: 基因名到Uniprot ID的映射
        json_folder_path: JSON文件夹路径
    
    返回:
        DataFrame，包含-window到+window位置的特征
    """
    
    # 查找Uniprot ID
    if gene_name not in gene_to_uniprot_dict:
        return None
    
    uniprot_id = gene_to_uniprot_dict[gene_name]
    json_filepath = os.path.join(json_folder_path, f"{uniprot_id}.json")
    
    if not os.path.exists(json_filepath):
        return None
    
    # 读取JSON
    try:
        with open(json_filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except:
        return None
    
    residues = data.get('residues', [])
    
    # 提取-window到+window的特征
    results = []
    for offset in range(-window, window+1):
        target_pos = position + offset
        
        # 检查位置是否有效（1-based转0-based）
        target_idx = target_pos - 1
        
        if target_idx < 0 or target_idx >= len(residues):
            # 超出范围，填充NaN
            results.append({
                'relative_position': offset,
                'absolute_position': target_pos,
                'backbone_dynamics': None,
                'sidechain_dynamics': None,
                'disorder_propensity': None,
                'helix_propensity': None,
                'sheet_propensity': None,
                'coil_propensity': None,
                'amino_acid': None
            })
        else:
            residue_data = residues[target_idx]
            results.append({
                'relative_position': offset,
                'absolute_position': target_pos,
                'backbone_dynamics': residue_data.get('backbone', None),
                'sidechain_dynamics': residue_data.get('sidechain', None),
                'disorder_propensity': residue_data.get('disoMine', None),
                'helix_propensity': residue_data.get('helix', None),
                'sheet_propensity': residue_data.get('sheet', None),
                'coil_propensity': residue_data.get('coil', None),
                'amino_acid': residue_data.get('aa', None)
            })
    
    return pd.DataFrame(results)

print("✓ 函数定义完成")
print("="*60)


# %%
# ============================================================
# Step 1: 加载必要的配置和映射
# ============================================================

import json
import os
import pandas as pd
import numpy as np

print("="*60)
print("Step 1: 加载配置")
print("="*60)

# 1. 加载基因到Uniprot的映射
with open('gene_to_uniprot.json', 'r', encoding='utf-8') as f:
    gene_to_uniprot = json.load(f)

print(f"✓ 加载基因映射: {len(gene_to_uniprot)} 个基因")

# 2. 设置JSON文件夹路径（根据你的实际情况修改）
# 常见位置：
possible_json_folders = [
    'D:/json',
    './json',
    '../json',
    'D:/博士/Phospho/json',
    './data/json'
]

JSON_FOLDER = None
for folder in possible_json_folders:
    if os.path.exists(folder):
        json_files = [f for f in os.listdir(folder) if f.endswith('.json')]
        if len(json_files) > 0:
            JSON_FOLDER = folder
            print(f"✓ 找到JSON文件夹: {folder} ({len(json_files)} 个文件)")
            break

if JSON_FOLDER is None:
    print("\n⚠️  未找到JSON文件夹，请手动指定:")
    print("当前目录的子文件夹:")
    subdirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    for d in subdirs:
        print(f"  - {d}")
    
    # 手动设置（如果自动查找失败）
    # JSON_FOLDER = "你的JSON文件夹路径"  # ← 取消注释并填写

# 3. 确认使用26个蛋白的子集
subset_df = subset_dfA_26.copy()

print(f"\n✓ 使用26个蛋白子集:")
print(f"  位点数: {len(subset_df)}")
print(f"  蛋白数: {subset_df['Gene'].nunique()}")
print(f"\n按LocalizationGroup分布:")
print(subset_df['LocalizationGroup'].value_counts())
print(f"\n按STY分布:")
print(subset_df['AA'].value_counts())

print("="*60)


# %%
# ============================================================
# Step 2: 定义±5序列上下文特征提取函数
# ============================================================

def extract_context_features(gene_name, position, window=5):
    """
    提取磷酸化位点前后window个氨基酸的生物物理特征
    """
    
    # 查找Uniprot ID
    if gene_name not in gene_to_uniprot:
        return None
    
    uniprot_id = gene_to_uniprot[gene_name]
    json_filepath = os.path.join(JSON_FOLDER, f"{uniprot_id}.json")
    
    if not os.path.exists(json_filepath):
        return None
    
    # 读取JSON
    try:
        with open(json_filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except:
        return None
    
    residues = data.get('residues', [])
    
    # 提取-window到+window的特征
    results = []
    for offset in range(-window, window+1):
        target_pos = position + offset
        target_idx = target_pos - 1  # 1-based转0-based
        
        if target_idx < 0 or target_idx >= len(residues):
            # 超出范围，填充NaN
            results.append({
                'relative_position': offset,
                'absolute_position': target_pos,
                'backbone_dynamics': np.nan,
                'disorder_propensity': np.nan,
                'helix_propensity': np.nan,
                'sheet_propensity': np.nan,
                'coil_propensity': np.nan,
                'amino_acid': None
            })
        else:
            residue_data = residues[target_idx]
            results.append({
                'relative_position': offset,
                'absolute_position': target_pos,
                'backbone_dynamics': residue_data.get('backbone', np.nan),
                'disorder_propensity': residue_data.get('disoMine', np.nan),
                'helix_propensity': residue_data.get('helix', np.nan),
                'sheet_propensity': residue_data.get('sheet', np.nan),
                'coil_propensity': residue_data.get('coil', np.nan),
                'amino_acid': residue_data.get('aa', None)
            })
    
    return pd.DataFrame(results)

print("✓ 函数定义完成")
print("="*60)


# %%
# ============================================================
# Step 3: 批量提取530个位点的±5特征
# ============================================================

from tqdm import tqdm

print("="*60)
print("Step 3: 批量提取序列上下文特征")
print("="*60)

if JSON_FOLDER is None:
    print("❌ 错误: JSON_FOLDER未设置，无法继续")
    print("请在Step 1中手动设置JSON_FOLDER路径")
else:
    print(f"开始提取 {len(subset_df)} 个位点的±5特征...\n")
    
    all_context_data = []
    failed_count = 0
    
    for idx, row in tqdm(subset_df.iterrows(), total=len(subset_df), desc="提取特征"):
        context_df = extract_context_features(
            gene_name=row['Gene'],
            position=row['Position'],
            window=5
        )
        
        if context_df is not None:
            # 添加位点元信息
            context_df['PTM_collapse_key'] = row['PTM_collapse_key']
            context_df['Gene'] = row['Gene']
            context_df['Site'] = row['Site']
            context_df['central_position'] = row['Position']
            context_df['LocalizationGroup'] = row['LocalizationGroup']
            context_df['AA_type'] = row['AA']
            
            all_context_data.append(context_df)
        else:
            failed_count += 1
    
    # 合并所有数据
    if len(all_context_data) > 0:
        context_full_df = pd.concat(all_context_data, ignore_index=True)
        
        print(f"\n✓ 提取完成!")
        print(f"  总行数: {len(context_full_df):,} (预期: {len(subset_df)}×11 = {len(subset_df)*11:,})")
        print(f"  成功位点: {context_full_df['PTM_collapse_key'].nunique()} / {len(subset_df)}")
        print(f"  失败位点: {failed_count}")
        
        # 保存结果
        OUTPUT_FILE = "26proteins_context_window5.csv"
        context_full_df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8-sig')
        print(f"  ✓ 已保存: {OUTPUT_FILE}")
    else:
        print("\n❌ 错误: 未能提取任何数据")
        print("请检查JSON文件夹路径是否正确")

print("="*60)


# %%
# ============================================================
# Step 4 修复版: 阶段3 - 使用正确的分组名称
# ============================================================

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

print("="*60)
print("Step 4修复版: 阶段3 - 序列上下文分析")
print("="*60)

context_df = pd.read_csv("26proteins_context_window5.csv")

# 特征列表
feature_cols = ['backbone_dynamics', 'disorder_propensity', 
                'helix_propensity', 'sheet_propensity', 'coil_propensity']

# ✅ 修复：使用实际的分组名称（无&符号）
groups = ['C-only', 'N-only', 'M-only', 'CN', 'CM', 'MN', 'CMN']

# 统计
print("\n位点分布:")
for group in groups:
    n = context_df[context_df['LocalizationGroup'] == group]['PTM_collapse_key'].nunique()
    print(f"  {group:10s}: {n:3d} sites")

print("\n开始生成图表...")

# 为每个LocalizationGroup生成图表
for group in groups:
    group_data = context_df[context_df['LocalizationGroup'] == group].copy()
    
    if len(group_data) == 0:
        print(f"\n⚠️  {group} 无数据，跳过")
        continue
    
    n_sites = group_data['PTM_collapse_key'].nunique()
    print(f"\n生成 {group} 的序列上下文图 (n={n_sites})...")
    
    # 创建5个子图
    fig, axes = plt.subplots(5, 1, figsize=(14, 20))
    fig.suptitle(f'Sequence Context Analysis: {group} (n={n_sites} sites)', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    for i, feature in enumerate(feature_cols):
        ax = axes[i]
        
        plot_data = group_data[['relative_position', feature]].dropna()
        
        if len(plot_data) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(feature.replace('_', ' ').title())
            continue
        
        positions = sorted(plot_data['relative_position'].unique())
        
        sns.boxplot(data=plot_data, x='relative_position', y=feature,
                   order=positions, palette='viridis', ax=ax)
        
        # 标记中心位点
        center_idx = positions.index(0) if 0 in positions else None
        if center_idx is not None:
            ax.axvline(x=center_idx, color='red', linestyle='--',
                      linewidth=2, alpha=0.7, label='Phosphosite')
        
        ax.set_title(f'{feature.replace("_", " ").title()}',
                    fontsize=12, fontweight='bold')
        ax.set_xlabel('Relative Position', fontsize=11)
        ax.set_ylabel('Feature Value', fontsize=11)
        ax.legend(loc='upper right')
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    # 保存
    figname = f"Stage3_{group.replace('-', '')}_context11pos.png"
    plt.savefig(figname, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  ✓ 保存: {figname}")
    plt.close()

print("\n" + "="*60)
print("✓ 阶段3完成！已生成7个分组的图表")
print("="*60)


# %%
# ============================================================
# Step 5 修复版: 阶段4 - STY × 序列上下文
# ============================================================

print("="*60)
print("Step 5修复版: 阶段4 - STY × 序列上下文")
print("="*60)

# ✅ 使用正确的分组名称
groups = ['C-only', 'N-only', 'M-only', 'CN', 'CM', 'MN', 'CMN']
aa_types = ['S', 'T', 'Y']

# 统计每个组合
print("\n组合统计（位点数）:")
print("-"*60)

for aa in aa_types:
    print(f"\n{aa}:")
    aa_data = context_df[context_df['AA_type'] == aa]
    
    for group in groups:
        n_sites = aa_data[aa_data['LocalizationGroup'] == group]['PTM_collapse_key'].nunique()
        if n_sites > 0:
            print(f"  {group:10s}: {n_sites:3d} sites")

print("\n" + "="*60)
print("开始生成图表（跳过样本数<3的组合）...\n")

# 生成图表
for aa in aa_types:
    print(f"【{aa} 类型】")
    
    aa_data = context_df[context_df['AA_type'] == aa].copy()
    
    if len(aa_data) == 0:
        print(f"  ⚠️  无数据，跳过\n")
        continue
    
    for group in groups:
        combo_data = aa_data[aa_data['LocalizationGroup'] == group].copy()
        
        if len(combo_data) == 0:
            continue
        
        n_sites = combo_data['PTM_collapse_key'].nunique()
        
        if n_sites < 3:
            print(f"  跳过 {aa}_{group} (n={n_sites}, 样本太少)")
            continue
        
        print(f"  生成 {aa}_{group} (n={n_sites})...")
        
        # 创建5个子图
        fig, axes = plt.subplots(5, 1, figsize=(14, 20))
        fig.suptitle(f'Sequence Context: {aa} × {group} (n={n_sites} sites)',
                     fontsize=16, fontweight='bold', y=0.995)
        
        for i, feature in enumerate(feature_cols):
            ax = axes[i]
            
            plot_data = combo_data[['relative_position', feature]].dropna()
            
            if len(plot_data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=14)
                ax.set_title(feature.replace('_', ' ').title())
                continue
            
            positions = sorted(plot_data['relative_position'].unique())
            
            sns.boxplot(data=plot_data, x='relative_position', y=feature,
                       order=positions, palette='coolwarm', ax=ax)
            
            center_idx = positions.index(0) if 0 in positions else None
            if center_idx is not None:
                ax.axvline(x=center_idx, color='red', linestyle='--',
                          linewidth=2, alpha=0.7, label='Phosphosite')
            
            ax.set_title(f'{feature.replace("_", " ").title()}',
                        fontsize=12, fontweight='bold')
            ax.set_xlabel('Relative Position', fontsize=11)
            ax.set_ylabel('Feature Value', fontsize=11)
            ax.legend(loc='upper right')
            ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        figname = f"Stage4_{aa}_{group.replace('-', '')}_context11pos.png"
        plt.savefig(figname, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"    ✓ 保存: {figname}")
        plt.close()
    
    print()  # 空行分隔

print("="*60)
print("✓ 阶段4完成！")
print("="*60)


# %%
# ============================================================
# 纵向对比版：S磷酸化 ±5序列上下文分析
# ============================================================

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

print("="*70)
print("【S磷酸化专项分析】±5序列上下文 (纵向对比7组)")
print("="*70)

# 读取数据
context_df = pd.read_csv("26proteins_context_window5.csv")
context_S = context_df[context_df['AA_type'] == 'S'].copy()

print(f"\n数据统计:")
print(f"  总位点数: {context_S['PTM_collapse_key'].nunique()} (S)")

# 分组
groups = ['C-only', 'N-only', 'M-only', 'CN', 'CM', 'MN', 'CMN']
group_labels = ['C-only', 'N-only', 'M-only', 'C&N', 'C&M', 'N&M', 'C&M&N']

# 统计
group_stats = {}
for group in groups:
    n_sites = context_S[context_S['LocalizationGroup'] == group]['PTM_collapse_key'].nunique()
    group_stats[group] = n_sites
    print(f"  {group:10s}: {n_sites:3d} sites")

# 保留有效组
valid_groups = [g for g in groups if group_stats[g] >= 3]
valid_labels = [group_labels[groups.index(g)] for g in valid_groups]

# 3个显著特征
significant_features = [
    'backbone_dynamics',
    'helix_propensity', 
    'coil_propensity'
]

feature_display_names = {
    'backbone_dynamics': 'Backbone Dynamics',
    'helix_propensity': 'Helix Propensity',
    'coil_propensity': 'Coil Propensity'
}

print("\n" + "="*70)
print("生成纵向对比图（7行×11列）...")
print("="*70)

# 为每个特征生成一张图
for feature in significant_features:
    print(f"\n【{feature_display_names[feature]}】生成中...")
    
    # 创建子图：N行×1列（N=有效组数）
    fig, axes = plt.subplots(len(valid_groups), 1, 
                            figsize=(14, 2.5*len(valid_groups)),
                            sharex=True)
    
    # 如果只有一个组
    if len(valid_groups) == 1:
        axes = [axes]
    
    # 总标题
    fig.suptitle(f'{feature_display_names[feature]} - Sequence Context (±5 aa)\nS-phosphorylation Sites Across LocalizationGroups', 
                 fontsize=15, fontweight='bold', y=0.995)
    
    # 为每个组绘制子图
    for idx, (group, label) in enumerate(zip(valid_groups, valid_labels)):
        ax = axes[idx]
        
        # 提取数据
        group_data = context_S[context_S['LocalizationGroup'] == group].copy()
        plot_data = group_data[['relative_position', feature]].dropna()
        
        if len(plot_data) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=12, color='gray')
            ax.set_ylabel(f'{label}\n(n=0)', fontsize=11, fontweight='bold', rotation=0, 
                         ha='right', va='center')
            continue
        
        # 位置
        positions = sorted(plot_data['relative_position'].unique())
        
        # 绘制boxplot
        bp = sns.boxplot(data=plot_data, x='relative_position', y=feature,
                        order=positions,
                        palette='RdYlBu_r',
                        ax=ax,
                        width=0.7,
                        linewidth=1.5,
                        fliersize=3)
        
        # 标记中心位点
        if 0 in positions:
            center_idx = positions.index(0)
            
            # 红色虚线
            ax.axvline(x=center_idx, color='#d62728', linestyle='--',
                      linewidth=2.8, alpha=0.85, zorder=10)
            
            # 浅红色背景
            ax.axvspan(center_idx-0.45, center_idx+0.45, 
                      alpha=0.15, color='#d62728', zorder=0)
        
        # 计算统计
        center_data = plot_data[plot_data['relative_position'] == 0][feature]
        flanking_data = plot_data[plot_data['relative_position'] != 0][feature]
        
        center_mean = center_data.mean() if len(center_data) > 0 else np.nan
        flanking_mean = flanking_data.mean() if len(flanking_data) > 0 else np.nan
        delta = center_mean - flanking_mean if not np.isnan(center_mean) else np.nan
        
        # 统计检验
        if len(center_data) > 0 and len(flanking_data) > 0:
            stat, pval = stats.mannwhitneyu(center_data, flanking_data, alternative='two-sided')
            sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
        else:
            pval = np.nan
            sig = 'N/A'
        
        # Y轴标签（组名 + 样本量 + 统计）
        n = group_stats[group]
        ylabel_text = f'{label}\n(n={n})\n{sig}'
        ax.set_ylabel(ylabel_text, fontsize=10, fontweight='bold', 
                     rotation=0, ha='right', va='center', labelpad=45)
        
        # X轴标签（只在最后一个子图显示）
        if idx == len(valid_groups) - 1:
            ax.set_xlabel('Position Relative to S-phosphosite', fontsize=11, fontweight='bold')
            ax.set_xticklabels([f'{int(p):+d}' if p != 0 else 'S*' for p in positions],
                              fontsize=10, fontweight='bold')
        else:
            ax.set_xlabel('')
        
        # 网格
        ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=0.8)
        ax.grid(axis='x', alpha=0.15, linestyle='-', linewidth=0.5)
        
        # 添加统计信息（右上角）
        if not np.isnan(center_mean):
            info_text = f'S*: {center_mean:.3f}\nFlanking: {flanking_mean:.3f}\nΔ: {delta:+.3f}'
            if not np.isnan(pval):
                info_text += f'\np={pval:.2e}'
            
            # 根据显著性设置背景色
            if sig == '***':
                bgcolor = '#d4edda'  # 浅绿
            elif sig == '**':
                bgcolor = '#fff3cd'  # 浅黄
            elif sig == '*':
                bgcolor = '#fff3cd'
            else:
                bgcolor = '#f8f9fa'  # 浅灰
            
            ax.text(0.98, 0.97, info_text,
                   transform=ax.transAxes,
                   ha='right', va='top',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor=bgcolor, 
                           edgecolor='gray', linewidth=1, alpha=0.8),
                   fontsize=8.5,
                   family='monospace')
        
        # 设置Y轴范围一致（方便对比）
        # 这个可以根据需要注释掉
        # ax.set_ylim([0, 1])
    
    plt.tight_layout()
    
    # 保存
    figname = f"S_Context_Vertical_{feature}.png"
    plt.savefig(figname, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  ✓ 已保存: {figname}")
    plt.close()

print("\n" + "="*70)
print("✅ 分析完成！共生成 3 张纵向对比图")
print("="*70)

# 生成统计总结表
print("\n" + "="*70)
print("【统计总结】中心位点(S*) vs 周围氨基酸")
print("="*70)

summary_data = []

for feature in significant_features:
    print(f"\n{feature_display_names[feature]}:")
    print("-" * 70)
    print(f"{'Group':<12s} {'n':<5s} {'S* mean':<10s} {'Flanking':<10s} {'Δ':<10s} {'p-value':<12s} {'Sig'}")
    print("-" * 70)
    
    for group in valid_groups:
        group_data = context_S[context_S['LocalizationGroup'] == group]
        plot_data = group_data[['relative_position', feature]].dropna()
        
        if len(plot_data) == 0:
            continue
        
        center_data = plot_data[plot_data['relative_position'] == 0][feature]
        flanking_data = plot_data[plot_data['relative_position'] != 0][feature]
        
        if len(center_data) > 0 and len(flanking_data) > 0:
            center_mean = center_data.mean()
            flanking_mean = flanking_data.mean()
            delta = center_mean - flanking_mean
            
            stat, pval = stats.mannwhitneyu(center_data, flanking_data, alternative='two-sided')
            sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
            
            label = group_labels[groups.index(group)]
            n = group_stats[group]
            
            print(f"{label:<12s} {n:<5d} {center_mean:<10.4f} {flanking_mean:<10.4f} "
                  f"{delta:<+10.4f} {pval:<12.2e} {sig}")
            
            summary_data.append({
                'Feature': feature,
                'Group': label,
                'n': n,
                'Center_mean': center_mean,
                'Flanking_mean': flanking_mean,
                'Delta': delta,
                'p_value': pval,
                'Significant': sig
            })

# 保存统计表
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('S_Context_Statistical_Summary.csv', index=False)
print("\n✓ 统计表已保存: S_Context_Statistical_Summary.csv")

print("\n" + "="*70)
print("【生成的图表】")
print("="*70)
for feature in significant_features:
    figname = f"S_Context_Vertical_{feature}.png"
    print(f"  ✓ {figname}")

print("\n" + "="*70)
print("【PPT建议】")
print("="*70)
print("""
纵向对比图的优势：
  ✓ 容易看出哪些组的S*位置特殊
  ✓ 容易看出组间pattern是否相似
  ✓ 一张图展示全部7组

如果用于PPT Slide 1:
  • 选择最显著的一张（通常是backbone_dynamics）
  • 指出："±5 analysis shows modest effects"
  • 强调："Pattern similar across groups → Single-site already captures signal"

如果所有组的Δ都 < 0.05:
  → "Local context adds minor refinement, not game-changing"
  
如果某些组Δ > 0.1:
  → "Context effects present but group-specific"
""")

print("="*70)


# %%
# 在你现有notebook的最后一个cell后面运行
print("当前success_df的列：")
print(success_df.columns.tolist())
print(f"\n行数: {len(success_df)}")
print(f"\n前3行示例:")
print(success_df[['PTM_collapse_key', 'Gene', 'Position', 'LocalizationGroup']].head(3))

# %%
# ===== Cell 20: 读取Intensity数据并计算Compartment平均值 =====

import pandas as pd
import numpy as np
import openpyxl

print("=" * 80)
print("Cell 20: 读取Intensity数据并计算Compartment平均值")
print("=" * 80)

# ========== 1. 读取intensity数据 ==========
print("\n[1/5] 读取intensity数据...")
df_intensity = pd.read_excel(
    r"D:\博士\Phospho\Full_data.xlsx",
    sheet_name="S5D-Log2 proc HeLa+EGF PHOS",
    engine='openpyxl'
)
print(f"✓ 读取成功: {len(df_intensity)} 行 × {len(df_intensity.columns)} 列")

# ========== 2. 定义Fraction到Compartment的映射 ==========
print("\n[2/5] 定义Fraction到Compartment的映射...")
# 根据原文：
# FR1+FR2 → Cytoplasm
# FR3+FR4 → Membrane-bound organelles
# FR5+FR6 → Nucleus

fraction_to_compartment = {
    'FR1': 'Cyt',
    'FR2': 'Cyt',
    'FR3': 'Mem',
    'FR4': 'Mem',
    'FR5': 'Nuc',
    'FR6': 'Nuc'
}

# ========== 3. 计算每个时间点、每个compartment的平均intensity ==========
print("\n[3/5] 计算平均intensity...")

# 时间点映射（CTRL = 0min）
time_points = ['CTRL', '2min', '8min', '20min', '90min']
time_labels = ['0min', '2min', '8min', '20min', '90min']

# 初始化结果字典
intensity_data = {'PTM_collapse_key': df_intensity['PTM_collapse_key']}

for time_point, time_label in zip(time_points, time_labels):
    print(f"  处理 {time_label}...")
    
    for compartment in ['Cyt', 'Mem', 'Nuc']:
        # 找出属于该compartment的所有fractions
        fractions = [fr for fr, comp in fraction_to_compartment.items() if comp == compartment]
        
        # 收集所有replicates的列名
        all_cols = []
        for fr in fractions:
            for rep in range(1, 5):  # Rep1-Rep4
                col_name = f"EGF_{time_point}_{fr}_Rep{rep}"
                if col_name in df_intensity.columns:
                    all_cols.append(col_name)
        
        # 计算平均值
        if len(all_cols) > 0:
            # 对所有replicates和fractions求平均
            intensity_data[f'{time_label}_{compartment}_intensity'] = df_intensity[all_cols].mean(axis=1)
        else:
            print(f"    ⚠️ 未找到 {time_label}_{compartment} 的数据列")
            intensity_data[f'{time_label}_{compartment}_intensity'] = 0.0

# 转换为DataFrame
df_intensity_processed = pd.DataFrame(intensity_data)

print(f"\n✓ 处理完成，生成 {len(df_intensity_processed)} 行 × {len(df_intensity_processed.columns)} 列")
print(f"  新增列: {[col for col in df_intensity_processed.columns if col != 'PTM_collapse_key']}")

# ========== 4. 和success_df合并 ==========
print("\n[4/5] 和success_df合并...")
print(f"  合并前 success_df: {len(success_df)} 行")
print(f"  合并前 df_intensity_processed: {len(df_intensity_processed)} 行")

# 检查匹配率
matched_keys = set(success_df['PTM_collapse_key']) & set(df_intensity_processed['PTM_collapse_key'])
print(f"  可匹配的PTM_collapse_key数量: {len(matched_keys)}")

# 左连接：保留success_df的所有行
success_df = success_df.merge(
    df_intensity_processed,
    on='PTM_collapse_key',
    how='left'
)

print(f"  合并后 success_df: {len(success_df)} 行 × {len(success_df.columns)} 列")

# 检查合并后的缺失值
intensity_cols = [col for col in success_df.columns if '_intensity' in col]
missing_intensity = success_df[intensity_cols].isnull().any(axis=1).sum()
print(f"  缺失intensity数据的位点: {missing_intensity} 个")

# ========== 5. 数据质量检查 ==========
print("\n[5/5] 数据质量检查...")

# 显示一个完整案例
print("\n示例位点 (AAK1_S623_M2):")
example_row = success_df[success_df['PTM_collapse_key'] == 'AAK1_S623_M2'].iloc[0]
for time_label in time_labels:
    cyt_int = example_row[f'{time_label}_Cyt_intensity']
    mem_int = example_row[f'{time_label}_Mem_intensity']
    nuc_int = example_row[f'{time_label}_Nuc_intensity']
    print(f"  {time_label}: Cyt={cyt_int:.2f}, Mem={mem_int:.2f}, Nuc={nuc_int:.2f}")

# 统计各compartment在各时间点的intensity分布
print("\n各时间点Intensity统计（非零值）:")
for time_label in time_labels:
    for comp in ['Cyt', 'Mem', 'Nuc']:
        col = f'{time_label}_{comp}_intensity'
        non_zero = success_df[success_df[col] > 0][col]
        if len(non_zero) > 0:
            print(f"  {col:30s}: n={len(non_zero):5d}, mean={non_zero.mean():6.2f}, std={non_zero.std():5.2f}")

print("\n" + "=" * 80)
print("✓ Cell 20 执行完成！")
print("=" * 80)
print("\n生成的新列（intensity相关）:")
new_cols = [col for col in success_df.columns if '_intensity' in col]
for col in new_cols:
    print(f"  • {col}")
print(f"\n当前 success_df: {len(success_df)} 行 × {len(success_df.columns)} 列")
print("=" * 80)


# %%
# ===== Cell 21: 计算Movement Score =====

import pandas as pd
import numpy as np
from scipy import stats

print("=" * 80)
print("Cell 21: 计算Movement Score")
print("=" * 80)

# ========== 1. 标准化intensity为百分比分布 ==========
print("\n[1/4] 标准化intensity为百分比分布...")

time_labels = ['0min', '2min', '8min', '20min', '90min']
compartments = ['Cyt', 'Mem', 'Nuc']

# 计算每个时间点的总intensity（三个compartment之和）
for time_label in time_labels:
    total_col = f'{time_label}_total_intensity'
    success_df[total_col] = (
        success_df[f'{time_label}_Cyt_intensity'] + 
        success_df[f'{time_label}_Mem_intensity'] + 
        success_df[f'{time_label}_Nuc_intensity']
    )
    
    # 计算百分比（避免除零错误）
    for comp in compartments:
        intensity_col = f'{time_label}_{comp}_intensity'
        pct_col = f'{time_label}_{comp}_pct'
        
        # 如果total=0，百分比设为0
        success_df[pct_col] = np.where(
            success_df[total_col] > 0,
            success_df[intensity_col] / success_df[total_col],
            0.0
        )

print("✓ 标准化完成，生成15个百分比列")

# ========== 2. 计算Movement Score ==========
print("\n[2/4] 计算Movement Score...")

def calculate_movement_score_for_site(row):
    """
    计算单个位点的Movement Score
    
    返回：
    - max_mobility: 最大的mobility score (0-200%)
    - best_time: 变化最大的时间点
    - direction_from: 主要从哪个compartment离开
    - direction_to: 主要去哪个compartment
    """
    
    # 如果0min的total intensity为0，无法计算
    if row['0min_total_intensity'] == 0:
        return {
            'movement_score': np.nan,
            'movement_best_time': 'undetectable',
            'direction_from': 'None',
            'direction_to': 'None'
        }
    
    # 获取0min的百分比分布
    baseline = {
        'Cyt': row['0min_Cyt_pct'],
        'Mem': row['0min_Mem_pct'],
        'Nuc': row['0min_Nuc_pct']
    }
    
    max_mobility = 0
    best_time = '0min'
    best_deltas = {'Cyt': 0, 'Mem': 0, 'Nuc': 0}
    
    # 遍历其他时间点
    for time_label in ['2min', '8min', '20min', '90min']:
        # 如果该时间点total=0，跳过
        if row[f'{time_label}_total_intensity'] == 0:
            continue
        
        # 计算各compartment的变化
        deltas = {}
        for comp in compartments:
            delta = abs(row[f'{time_label}_{comp}_pct'] - baseline[comp])
            deltas[comp] = delta
        
        # 取最大和次大的delta之和
        sorted_deltas = sorted(deltas.values(), reverse=True)
        mobility = (sorted_deltas[0] + sorted_deltas[1]) * 100  # 转成百分比
        
        # 更新最大值
        if mobility > max_mobility:
            max_mobility = mobility
            best_time = time_label
            best_deltas = deltas
    
    # 判断移动方向（哪个增加最多 vs 哪个减少最多）
    if max_mobility > 0:
        increases = {comp: row[f'{best_time}_{comp}_pct'] - baseline[comp] 
                    for comp in compartments}
        
        # 找增加最多的
        direction_to = max(increases, key=increases.get)
        # 找减少最多的
        direction_from = min(increases, key=increases.get)
        
        # 如果变化太小（<1%），标记为None
        if abs(increases[direction_to]) < 0.01:
            direction_to = 'None'
        if abs(increases[direction_from]) < 0.01:
            direction_from = 'None'
    else:
        direction_from = 'None'
        direction_to = 'None'
    
    return {
        'movement_score': max_mobility,
        'movement_best_time': best_time,
        'direction_from': direction_from,
        'direction_to': direction_to
    }

# 应用到每一行
print("  正在计算每个位点的Movement Score...")
movement_results = success_df.apply(calculate_movement_score_for_site, axis=1, result_type='expand')

# 合并结果
success_df['movement_score'] = movement_results['movement_score']
success_df['movement_best_time'] = movement_results['movement_best_time']
success_df['direction_from'] = movement_results['direction_from']
success_df['direction_to'] = movement_results['direction_to']

print("✓ Movement Score计算完成")

# ========== 3. 数据质量检查 ==========
print("\n[3/4] Movement Score统计...")

# 排除undetectable的位点
valid_scores = success_df[success_df['movement_score'].notna()]['movement_score']

print(f"\n可计算Movement Score的位点: {len(valid_scores):,} / {len(success_df):,}")
print(f"无法计算的位点（0min未检测到）: {success_df['movement_score'].isna().sum():,}")

print(f"\nMovement Score分布:")
print(f"  最小值: {valid_scores.min():.2f}%")
print(f"  25分位: {valid_scores.quantile(0.25):.2f}%")
print(f"  中位数: {valid_scores.median():.2f}%")
print(f"  75分位: {valid_scores.quantile(0.75):.2f}%")
print(f"  最大值: {valid_scores.max():.2f}%")
print(f"  平均值: {valid_scores.mean():.2f}%")

# 按Movement Score大小分段统计
bins = [0, 1, 5, 10, 20, 50, 200]
labels = ['0-1%', '1-5%', '5-10%', '10-20%', '20-50%', '>50%']
score_bins = pd.cut(valid_scores, bins=bins, labels=labels)
print(f"\nMovement Score分段统计:")
for label in labels:
    count = (score_bins == label).sum()
    pct = count / len(valid_scores) * 100
    print(f"  {label:10s}: {count:5,} ({pct:5.1f}%)")

# 变化最大的时间点统计
print(f"\n变化最大的时间点统计:")
time_counts = success_df['movement_best_time'].value_counts()
for time_label, count in time_counts.items():
    print(f"  {time_label:15s}: {count:5,} ({count/len(success_df)*100:5.1f}%)")

# 移动方向统计（只看movement_score >= 10的）
high_movement = success_df[success_df['movement_score'] >= 10]
if len(high_movement) > 0:
    print(f"\n高移动位点（≥10%）的方向统计 (n={len(high_movement)}):")
    direction_pairs = high_movement.groupby(['direction_from', 'direction_to']).size().sort_values(ascending=False)
    for (from_comp, to_comp), count in direction_pairs.head(10).items():
        print(f"  {from_comp} → {to_comp}: {count:4,}")

# ========== 4. 展示典型案例 ==========
print("\n[4/4] 典型案例展示...")

# 找movement score最高的3个位点
top_movers = success_df.nlargest(3, 'movement_score')

print("\n✨ Movement Score最高的3个位点:")
for idx, row in top_movers.iterrows():
    print(f"\n{row['PTM_collapse_key']} (Movement Score = {row['movement_score']:.1f}%)")
    print(f"  基因: {row['Gene']}, 位点: {row['Site']}")
    print(f"  变化最大时刻: {row['movement_best_time']}")
    print(f"  移动方向: {row['direction_from']} → {row['direction_to']}")
    print(f"  0/1 Pattern: {row['LocalizationGroup']}")
    print(f"  Trajectory: {row['trajectory']}")
    
    # 显示百分比变化
    print(f"  百分比分布:")
    for time_label in time_labels:
        cyt = row[f'{time_label}_Cyt_pct'] * 100
        mem = row[f'{time_label}_Mem_pct'] * 100
        nuc = row[f'{time_label}_Nuc_pct'] * 100
        print(f"    {time_label}: Cyt={cyt:5.1f}%, Mem={mem:5.1f}%, Nuc={nuc:5.1f}%")

print("\n" + "=" * 80)
print("✓ Cell 21 执行完成！")
print("=" * 80)
print(f"\n新增列:")
print(f"  • movement_score (Movement Score)")
print(f"  • movement_best_time (变化最大的时间点)")
print(f"  • direction_from (主要从哪离开)")
print(f"  • direction_to (主要去哪)")
print(f"  • 15个百分比列 (0min_Cyt_pct, ...)")
print(f"\n当前 success_df: {len(success_df)} 行 × {len(success_df.columns)} 列")
print("=" * 80)


# %%
# ===== Cell 22: 587个高mobility位点的详细统计 =====

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

print("=" * 80)
print("Cell 22: 587个高mobility位点的全面统计")
print("=" * 80)

# ========== 1. 筛选高mobility位点 ==========
print("\n[1/7] 筛选高mobility位点...")
high_mobility = success_df[success_df['movement_score'] >= 10].copy()
print(f"✓ 筛选出 {len(high_mobility)} 个位点 (movement_score ≥ 10%)")

# ========== 2. 维度A：移动方向（空间维度）==========
print("\n" + "=" * 80)
print("[2/7] 维度A：移动方向（空间维度）")
print("=" * 80)

# 2.1 统计各方向的数量
direction_counts = high_mobility.groupby(['direction_from', 'direction_to']).size().reset_index(name='count')
direction_counts['direction'] = direction_counts['direction_from'] + '→' + direction_counts['direction_to']
direction_counts = direction_counts.sort_values('count', ascending=False)

print("\n移动方向分布:")
for idx, row in direction_counts.iterrows():
    direction = row['direction']
    count = row['count']
    pct = count / len(high_mobility) * 100
    print(f"  {direction:15s}: {count:4d} ({pct:5.1f}%)")

# 2.2 创建方向标签列
high_mobility['Direction'] = high_mobility['direction_from'] + '→' + high_mobility['direction_to']

# ========== 3. 维度B：氨基酸类型（phospho特异）==========
print("\n" + "=" * 80)
print("[3/7] 维度B：氨基酸类型")
print("=" * 80)

aa_counts = high_mobility['AA'].value_counts()
print("\n氨基酸分布:")
for aa, count in aa_counts.items():
    pct = count / len(high_mobility) * 100
    print(f"  {aa}: {count:4d} ({pct:5.1f}%)")

# ========== 4. 维度C：响应时间（时间维度）==========
print("\n" + "=" * 80)
print("[4/7] 维度C：响应时间")
print("=" * 80)

time_counts = high_mobility['movement_best_time'].value_counts().sort_index()
print("\n变化最大的时间点分布:")
for time, count in time_counts.items():
    pct = count / len(high_mobility) * 100
    print(f"  {time:10s}: {count:4d} ({pct:5.1f}%)")

# ========== 5. Movement Score分布 ==========
print("\n" + "=" * 80)
print("[5/7] Movement Score分布")
print("=" * 80)

score_stats = high_mobility['movement_score'].describe()
print("\nMovement Score统计:")
print(f"  最小值:  {score_stats['min']:.2f}%")
print(f"  25分位:  {score_stats['25%']:.2f}%")
print(f"  中位数:  {score_stats['50%']:.2f}%")
print(f"  75分位:  {score_stats['75%']:.2f}%")
print(f"  最大值:  {score_stats['max']:.2f}%")
print(f"  平均值:  {score_stats['mean']:.2f}%")

# 分档统计
score_bins = pd.cut(high_mobility['movement_score'], 
                    bins=[10, 15, 20, 30, 50, 200],
                    labels=['10-15%', '15-20%', '20-30%', '30-50%', '>50%'])
print("\nMovement Score分档:")
for label in score_bins.cat.categories:
    count = (score_bins == label).sum()
    pct = count / len(high_mobility) * 100
    print(f"  {label:10s}: {count:4d} ({pct:5.1f}%)")

# ========== 6. 交叉统计：方向 × 氨基酸 ==========
print("\n" + "=" * 80)
print("[6/7] 交叉统计：移动方向 × 氨基酸")
print("=" * 80)

cross_direction_aa = pd.crosstab(high_mobility['Direction'], high_mobility['AA'])
print("\n各方向中的氨基酸分布:")
print(cross_direction_aa.to_string())

# 计算每个方向的氨基酸比例
cross_direction_aa_pct = cross_direction_aa.div(cross_direction_aa.sum(axis=1), axis=0) * 100
print("\n各方向中的氨基酸比例(%):")
print(cross_direction_aa_pct.round(1).to_string())

# ========== 7. 交叉统计：方向 × 时间 ==========
print("\n" + "=" * 80)
print("[7/7] 交叉统计：移动方向 × 响应时间")
print("=" * 80)

cross_direction_time = pd.crosstab(high_mobility['Direction'], high_mobility['movement_best_time'])
print("\n各方向在不同时间点的分布:")
print(cross_direction_time.to_string())

# 计算每个方向的时间比例
cross_direction_time_pct = cross_direction_time.div(cross_direction_time.sum(axis=1), axis=0) * 100
print("\n各方向在不同时间点的比例(%):")
print(cross_direction_time_pct.round(1).to_string())

# ========== 8. 保存分组标签 ==========
print("\n" + "=" * 80)
print("保存分组信息到success_df...")
print("=" * 80)

# 为所有位点添加Direction和MobilityLevel列
success_df['Direction'] = success_df['direction_from'] + '→' + success_df['direction_to']

# Mobility分级（三档）
def assign_mobility_level(score):
    if pd.isna(score):
        return 'Undetectable'
    elif score >= 10:
        return 'High'
    elif score >= 5:
        return 'Medium'
    else:
        return 'Low'

success_df['MobilityLevel'] = success_df['movement_score'].apply(assign_mobility_level)

# 统计各mobility level的数量
mobility_dist = success_df['MobilityLevel'].value_counts()
print("\n所有位点的Mobility分级:")
for level in ['High', 'Medium', 'Low', 'Undetectable']:
    if level in mobility_dist.index:
        count = mobility_dist[level]
        pct = count / len(success_df) * 100
        print(f"  {level:15s}: {count:5,} ({pct:5.1f}%)")

print("\n" + "=" * 80)
print("✓ Cell 22 执行完成！")
print("=" * 80)
print(f"\n关键发现:")
print(f"  • 587个高mobility位点分布在6个主要方向")
print(f"  • 最常见的方向: Nuc→Cyt (154个) 和 Mem→Cyt (146个)")
print(f"  • 响应时间跨度: 2min到90min")
print(f"  • 氨基酸分布: S/T/Y (待统计)")
print(f"\n新增列:")
print(f"  • Direction (移动方向标签)")
print(f"  • MobilityLevel (mobility分级: High/Medium/Low)")
print(f"\n当前 success_df: {len(success_df)} 行 × {len(success_df.columns)} 列")
print("=" * 80)


# %%
# ===== Cell 24: Movement Score分析 (自适应版本) =====

import pandas as pd
import numpy as np
from scipy import stats

print("="*80)
print("Movement Score分析: 高移动性 vs 低移动性位点的生物物理特征比较")
print("="*80)

# [0/5] 检查可用的生物物理特征
print("\n[0/5] 检查数据中的生物物理特征...")

# 定义期望的7个特征及其可能的列名变体
expected_features = {
    'backbone_dynamics': ['backbone_dynamics', 'Backbone_Dynamics', 'BackboneDynamics'],
    'sidechain_dynamics': ['sidechain_dynamics', 'Sidechain_Dynamics', 'SidechainDynamics', 'side_chain_dynamics'],
    'disorder_propensity': ['disorder_propensity', 'Disorder_Propensity', 'DisorderPropensity'],
    'helix_propensity': ['helix_propensity', 'Helix_Propensity', 'HelixPropensity'],
    'sheet_propensity': ['sheet_propensity', 'Sheet_Propensity', 'SheetPropensity'],
    'coil_propensity': ['coil_propensity', 'Coil_Propensity', 'CoilPropensity'],
    'earlyFolding': ['earlyFolding', 'early_folding', 'Early_Folding', 'EarlyFolding']
}

# 检查success_df中实际存在的列
actual_columns = set(success_df.columns)
bio_features = []
feature_mapping = {}

for feature_key, possible_names in expected_features.items():
    found = False
    for name in possible_names:
        if name in actual_columns:
            bio_features.append(name)
            feature_mapping[feature_key] = name
            found = True
            print(f"  ✓ {feature_key:25s} → {name}")
            break
    
    if not found:
        print(f"  ✗ {feature_key:25s} → 未找到")

print(f"\n可用特征数量: {len(bio_features)}/7")

if len(bio_features) == 0:
    print("\n❌ 错误: 没有找到任何生物物理特征列!")
    print("请检查success_df的列名。")
    print(f"\nsuccess_df的所有列:\n{list(success_df.columns)}")
    raise ValueError("没有可用的生物物理特征列")

# [1/5] 定义对比组
print("\n[1/5] 定义对比组...")

# 高移动性组: movement_score ≥ 10%
high_mobility = success_df[success_df['movement_score'] >= 10].copy()

# 低移动性组: movement_score < 5%
low_mobility = success_df[success_df['movement_score'] < 5].copy()

print(f"High mobility组 (movement_score ≥ 10%): n={len(high_mobility)}")
print(f"Low mobility组  (movement_score < 5%):  n={len(low_mobility)}")

# [2/5] 统计描述
print("\n[2/5] 生物物理特征统计...")

# 计算每组的均值和标准差
high_stats = high_mobility[bio_features].agg(['mean', 'std', 'count'])
low_stats = low_mobility[bio_features].agg(['mean', 'std', 'count'])

print("\n📊 High Mobility组 (n={})".format(len(high_mobility)))
print(high_stats.round(3))

print("\n📊 Low Mobility组 (n={})".format(len(low_mobility)))
print(low_stats.round(3))

# [3/5] 统计检验
print("\n[3/5] Mann-Whitney U检验...")

test_results = []
for feature in bio_features:
    high_data = high_mobility[feature].dropna()
    low_data = low_mobility[feature].dropna()
    
    if len(high_data) == 0 or len(low_data) == 0:
        print(f"  ⚠️  跳过 {feature}: 数据为空")
        continue
    
    # Mann-Whitney U检验
    statistic, pvalue = stats.mannwhitneyu(high_data, low_data, alternative='two-sided')
    
    # 效应量 (Cohen's d)
    mean_diff = high_data.mean() - low_data.mean()
    pooled_std = np.sqrt((high_data.std()**2 + low_data.std()**2) / 2)
    cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
    
    test_results.append({
        'Feature': feature,
        'High_mean': high_data.mean(),
        'Low_mean': low_data.mean(),
        'Mean_diff': mean_diff,
        'U_statistic': statistic,
        'p_value': pvalue,
        'Cohens_d': cohens_d,
        'Significant': '***' if pvalue < 0.001 else '**' if pvalue < 0.01 else '*' if pvalue < 0.05 else 'ns'
    })

test_df = pd.DataFrame(test_results)
print("\n统计检验结果:")
print(test_df.round(4))

# [4/5] 按氨基酸类型细分
print("\n[4/5] 按氨基酸类型(S/T/Y)细分...")

# 提取氨基酸类型
high_mobility['AA'] = high_mobility['Site'].str.extract(r'([STY])\d+')[0]
low_mobility['AA'] = low_mobility['Site'].str.extract(r'([STY])\d+')[0]

aa_distribution = pd.DataFrame({
    'High_mobility': high_mobility['AA'].value_counts(),
    'Low_mobility': low_mobility['AA'].value_counts()
})

print("\n氨基酸类型分布:")
print(aa_distribution)

# [5/5] 保存结果
print("\n[5/5] 保存分析结果...")

# 保存统计检验结果
test_df.to_csv('movement_score_statistical_tests.csv', index=False)
print("  ✓ 统计检验结果已保存: movement_score_statistical_tests.csv")

# 保存分组数据
high_mobility.to_csv('high_mobility_sites.csv', index=False)
low_mobility.to_csv('low_mobility_sites.csv', index=False)
print("  ✓ 分组数据已保存: high_mobility_sites.csv, low_mobility_sites.csv")

# 保存结果
print("\n" + "="*80)
print("✅ Movement Score分析完成!")
print("="*80)
print(f"\n关键发现:")
print(f"  • 高移动性位点: {len(high_mobility)} 个")
print(f"  • 低移动性位点: {len(low_mobility)} 个")
print(f"  • 可用特征数量: {len(bio_features)} 个")
print(f"  • 显著差异特征 (p<0.05): {len(test_df[test_df['p_value'] < 0.05])} 个")
print(f"  • 特征列表: {bio_features}")

# 提示缺失的特征
missing_features = []
for feature_key in expected_features.keys():
    if feature_key not in feature_mapping:
        missing_features.append(feature_key)

if missing_features:
    print(f"\n⚠️  缺失的特征: {missing_features}")
    print("  如果需要这些特征,请检查:")
    print("  1. 数据提取Cell (Cell 3-5)是否正确提取了这些特征")
    print("  2. 列名是否与预期一致")
    print("  3. 是否需要重新运行数据提取流程")
else:
    print(f"\n✅ 所有7个特征都已找到并分析!")


# %%
# ===== Cell 25: Movement Score箱线图 (自适应版本) =====

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

print("="*80)
print("生成Movement Score箱线图 (High vs Low Mobility)")
print("="*80)

# [0/2] 检查可用的生物物理特征
print("\n[0/2] 检查可用的生物物理特征...")

# 定义期望的7个特征及其可能的列名变体
expected_features = {
    'backbone_dynamics': ['backbone_dynamics', 'Backbone_Dynamics', 'BackboneDynamics'],
    'sidechain_dynamics': ['sidechain_dynamics', 'Sidechain_Dynamics', 'SidechainDynamics', 'side_chain_dynamics'],
    'disorder_propensity': ['disorder_propensity', 'Disorder_Propensity', 'DisorderPropensity'],
    'helix_propensity': ['helix_propensity', 'Helix_Propensity', 'HelixPropensity'],
    'sheet_propensity': ['sheet_propensity', 'Sheet_Propensity', 'SheetPropensity'],
    'coil_propensity': ['coil_propensity', 'Coil_Propensity', 'CoilPropensity'],
    'earlyFolding': ['earlyFolding', 'early_folding', 'Early_Folding', 'EarlyFolding']
}

# 特征的显示标签
feature_labels_dict = {
    'backbone_dynamics': 'Backbone\nDynamics',
    'sidechain_dynamics': 'Sidechain\nDynamics',
    'disorder_propensity': 'Disorder\nPropensity',
    'helix_propensity': 'Helix\nPropensity',
    'sheet_propensity': 'Sheet\nPropensity',
    'coil_propensity': 'Coil\nPropensity',
    'earlyFolding': 'Early\nFolding'
}

# 检查success_df中实际存在的列
actual_columns = set(success_df.columns)
bio_features = []
feature_labels = []
feature_mapping = {}

for feature_key, possible_names in expected_features.items():
    found = False
    for name in possible_names:
        if name in actual_columns:
            bio_features.append(name)
            feature_labels.append(feature_labels_dict[feature_key])
            feature_mapping[feature_key] = name
            found = True
            print(f"  ✓ {feature_key:25s} → {name}")
            break
    
    if not found:
        print(f"  ✗ {feature_key:25s} → 未找到")

print(f"\n可用特征数量: {len(bio_features)}/7")

if len(bio_features) == 0:
    print("\n❌ 错误: 没有找到任何生物物理特征列!")
    raise ValueError("没有可用的生物物理特征列")

# [1/2] 准备数据
print("\n[1/2] 准备数据...")

# 确保high_mobility和low_mobility已定义
if 'high_mobility' not in globals() or 'low_mobility' not in globals():
    print("  ⚠️  high_mobility或low_mobility未定义,正在重新定义...")
    high_mobility = success_df[success_df['movement_score'] >= 10].copy()
    low_mobility = success_df[success_df['movement_score'] < 5].copy()
    print(f"  ✓ High mobility: n={len(high_mobility)}")
    print(f"  ✓ Low mobility: n={len(low_mobility)}")

high_mobility['Mobility'] = 'High'
low_mobility['Mobility'] = 'Low'
combined_df = pd.concat([high_mobility, low_mobility], ignore_index=True)

# [2/2] 创建箱线图
print("\n[2/2] 创建箱线图...")

# 动态调整图表大小
n_features = len(bio_features)
fig_width = max(12, n_features * 3.5)  # 每个特征至少3.5英寸宽
fig, axes = plt.subplots(1, n_features, figsize=(fig_width, 5))

# 如果只有1个特征,axes不是数组,需要转换
if n_features == 1:
    axes = [axes]

fig.suptitle('Biophysical Features Comparison: High vs Low Mobility', 
             fontsize=16, fontweight='bold', y=1.02)

# 颜色设置
colors = {'High': '#E74C3C', 'Low': '#3498DB'}

# 绘制每个特征的箱线图
for i, (feature, label) in enumerate(zip(bio_features, feature_labels)):
    ax = axes[i]
    
    # 绘制箱线图
    sns.boxplot(data=combined_df, x=feature, y='Mobility', 
                palette=colors, ax=ax, orient='h')
    
    # 设置标题和标签
    ax.set_title(label, fontsize=12, fontweight='bold')
    ax.set_xlabel('Feature Value', fontsize=10)
    ax.set_ylabel('Mobility' if i == 0 else '', fontsize=10)
    
    # 添加样本量
    high_n = len(high_mobility[feature].dropna())
    low_n = len(low_mobility[feature].dropna())
    ax.text(0.98, 0.98, f'High: n={high_n}\nLow: n={low_n}',
            transform=ax.transAxes, fontsize=8,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    # 网格线
    ax.grid(axis='x', alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('movement_score_boxplots.png', dpi=300, bbox_inches='tight')
print("  ✓ 箱线图已保存: movement_score_boxplots.png")
plt.show()

print("\n✅ 箱线图生成完成!")
print(f"  • 特征数量: {len(bio_features)} 个")
print(f"  • 特征列表: {bio_features}")
print(f"  • High mobility: n={len(high_mobility)}")
print(f"  • Low mobility: n={len(low_mobility)}")

# 提示缺失的特征
missing_features = []
for feature_key in expected_features.keys():
    if feature_key not in feature_mapping:
        missing_features.append(feature_key)

if missing_features:
    print(f"\n⚠️  缺失的特征: {missing_features}")
    print("  箱线图只包含可用的特征。")
    print("  如果需要完整的7个特征,请检查数据提取流程。")
else:
    print(f"\n✅ 所有7个特征都已包含在箱线图中!")


# %%
"""
Cell 1: 数据准备与STY分组
功能: 将高/低移动性组按S/T/Y氨基酸类型拆分
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ============================================================================
# 1. 加载数据 (假设您已经运行了之前的代码,生成了final_df)
# ============================================================================
# 如果还没有final_df,请先运行之前的代码生成它
# 这里假设final_df已经在内存中

print("=" * 80)
print("📊 Cell 1: 数据准备与STY分组")
print("=" * 80)

# ============================================================================
# 2. 定义移动性分组逻辑
# ============================================================================
# 根据对话记录,高移动性定义为:在不同时间点和区室之间有转移的位点
# 这里使用ModScore或者根据0/1模式判断

# 定义区室列
compartment_cols = [
    '0min_Cyt', '0min_Mem', '0min_Nuc',
    '2min_Cyt', '2min_Mem', '2min_Nuc', 
    '8min_Cyt', '8min_Mem', '8min_Nuc',
    '20min_Cyt', '20min_Mem', '20min_Nuc',
    '90min_Cyt', '90min_Mem', '90min_Nuc'
]

# 计算每个位点的移动性得分 (出现的区室数量)
def calculate_mobility_score(row):
    """计算位点在不同区室出现的次数"""
    return row[compartment_cols].sum()

final_df['mobility_score'] = final_df.apply(calculate_mobility_score, axis=1)

# 定义高/低移动性阈值
# 根据对话,高移动性约587个,这里使用top 5-10%作为高移动性
mobility_threshold = final_df['mobility_score'].quantile(0.90)  # 可调整

print(f"\n📈 移动性统计:")
print(f"  移动性得分范围: {final_df['mobility_score'].min():.0f} - {final_df['mobility_score'].max():.0f}")
print(f"  移动性阈值 (90th percentile): {mobility_threshold:.1f}")

# 分组
final_df['mobility_group'] = final_df['mobility_score'].apply(
    lambda x: 'High' if x >= mobility_threshold else 'Low'
)

print(f"\n🔢 分组统计:")
print(f"  高移动性组: {(final_df['mobility_group'] == 'High').sum():,} 个位点")
print(f"  低移动性组: {(final_df['mobility_group'] == 'Low').sum():,} 个位点")

# ============================================================================
# 3. 按STY氨基酸类型拆分
# ============================================================================
# 使用AA列 (从原始数据) 或 amino_acid列 (从JSON提取)
aa_column = 'AA' if 'AA' in final_df.columns else 'amino_acid'

print(f"\n🧬 氨基酸类型统计:")
aa_counts = final_df[aa_column].value_counts()
for aa, count in aa_counts.items():
    print(f"  {aa}: {count:,} 个位点")

# 过滤掉缺失值
df_valid = final_df[final_df[aa_column].isin(['S', 'T', 'Y'])].copy()

print(f"\n✅ 有效数据: {len(df_valid):,} 个位点 (S/T/Y)")

# 按STY和mobility分组
groups = {}
for aa in ['S', 'T', 'Y']:
    for mobility in ['High', 'Low']:
        key = f"{aa}_{mobility}"
        groups[key] = df_valid[
            (df_valid[aa_column] == aa) & 
            (df_valid['mobility_group'] == mobility)
        ].copy()
        print(f"  {key}: {len(groups[key]):,} 个位点")

# ============================================================================
# 4. 提取生物物理特征列
# ============================================================================
biophysical_features = [
    'backbone_dynamics',
    'sidechain_dynamics',
    'disorder_propensity', 
    'helix_propensity',
    'sheet_propensity',
    'coil_propensity',
    'earlyFolding'
]

# 检查特征是否存在
missing_features = [f for f in biophysical_features if f not in df_valid.columns]
if missing_features:
    print(f"\n⚠️ 警告: 以下特征缺失: {missing_features}")
    biophysical_features = [f for f in biophysical_features if f in df_valid.columns]

print(f"\n📋 将分析的生物物理特征 ({len(biophysical_features)}个):")
for i, feature in enumerate(biophysical_features, 1):
    print(f"  {i}. {feature}")

# ============================================================================
# 5. 数据质量检查
# ============================================================================
print(f"\n🔍 数据质量检查:")
for feature in biophysical_features:
    valid_count = df_valid[feature].notna().sum()
    missing_count = df_valid[feature].isna().sum()
    print(f"  {feature:25s}: {valid_count:,} 有效 / {missing_count:,} 缺失")

print("\n" + "=" * 80)
print("✅ Cell 1 完成: 数据已按STY分组,共6个数据集准备就绪")
print("=" * 80)


# %%
"""
Cell 2: KS统计检验函数
功能: 实现随机抽样+KS检验,解决样本不平衡问题
"""

import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("📊 Cell 2: KS统计检验函数")
print("=" * 80)

# ============================================================================
# KS统计检验函数
# ============================================================================

def perform_ks_test_with_resampling(high_mobility_data, low_mobility_data, 
                                     n_iterations=10, random_state=42):
    """
    使用随机抽样进行KS检验,解决样本不平衡问题
    
    参数:
        high_mobility_data: 高移动性组的数据 (Series或array)
        low_mobility_data: 低移动性组的数据 (Series或array)
        n_iterations: 随机抽样的次数 (默认10次)
        random_state: 随机种子
        
    返回:
        dict: 包含KS统计量、p值、平均KS统计量等信息
    """
    
    # 移除缺失值
    high_data = high_mobility_data.dropna() if hasattr(high_mobility_data, 'dropna') else high_mobility_data[~np.isnan(high_mobility_data)]
    low_data = low_mobility_data.dropna() if hasattr(low_mobility_data, 'dropna') else low_mobility_data[~np.isnan(low_mobility_data)]
    
    # 检查数据是否足够
    if len(high_data) < 5 or len(low_data) < 5:
        return {
            'ks_statistic': np.nan,
            'p_value': np.nan,
            'mean_ks': np.nan,
            'std_ks': np.nan,
            'significant': False,
            'error': 'insufficient_data'
        }
    
    # 样本量
    n_high = len(high_data)
    n_low = len(low_data)
    
    # 如果高移动性组样本量更大,交换
    if n_high > n_low:
        high_data, low_data = low_data, high_data
        n_high, n_low = n_low, n_high
    
    # 进行多次随机抽样
    np.random.seed(random_state)
    ks_statistics = []
    p_values = []
    
    for i in range(n_iterations):
        # 从大样本组随机抽样,匹配小样本组的数量
        sampled_low = np.random.choice(low_data, size=n_high, replace=False)
        
        # 进行KS检验
        ks_stat, p_val = stats.ks_2samp(high_data, sampled_low)
        ks_statistics.append(ks_stat)
        p_values.append(p_val)
    
    # 计算平均值
    mean_ks = np.mean(ks_statistics)
    std_ks = np.std(ks_statistics)
    mean_p = np.mean(p_values)
    
    # 判断显著性 (使用平均p值)
    if mean_p < 0.001:
        sig_level = '***'
        significant = True
    elif mean_p < 0.01:
        sig_level = '**'
        significant = True
    elif mean_p < 0.05:
        sig_level = '*'
        significant = True
    else:
        sig_level = 'ns'
        significant = False
    
    return {
        'ks_statistic': mean_ks,
        'p_value': mean_p,
        'mean_ks': mean_ks,
        'std_ks': std_ks,
        'all_ks': ks_statistics,
        'all_p': p_values,
        'significant': significant,
        'sig_level': sig_level,
        'n_high': n_high,
        'n_low': n_low,
        'n_iterations': n_iterations
    }


def format_pvalue(p):
    """格式化p值显示"""
    if np.isnan(p):
        return "N/A"
    elif p < 0.001:
        return "p < 0.001"
    elif p < 0.01:
        return f"p = {p:.3f}"
    else:
        return f"p = {p:.2f}"


# ============================================================================
# 测试KS检验函数
# ============================================================================

print("\n🧪 测试KS检验函数...")

# 生成测试数据
np.random.seed(42)
test_high = np.random.normal(0.5, 0.1, 100)  # 高移动性组: 均值0.5
test_low = np.random.normal(0.4, 0.1, 1000)  # 低移动性组: 均值0.4

# 执行KS检验
result = perform_ks_test_with_resampling(test_high, test_low, n_iterations=10)

print(f"\n📈 测试结果:")
print(f"  高移动性组样本量: {result['n_high']}")
print(f"  低移动性组样本量: {result['n_low']}")
print(f"  KS统计量 (平均): {result['ks_statistic']:.4f} ± {result['std_ks']:.4f}")
print(f"  P值 (平均): {format_pvalue(result['p_value'])}")
print(f"  显著性: {result['sig_level']}")
print(f"  是否显著: {'是' if result['significant'] else '否'}")

print("\n" + "=" * 80)
print("✅ Cell 2 完成: KS检验函数已准备就绪")
print("=" * 80)


# %%
"""
Cell 3: 按STY拆分的Box Plot + KS统计检验
功能: 生成1行3列的subplot,比较S/T/Y的高/低移动性组
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

print("=" * 80)
print("📊 Cell 3: 按STY拆分的Box Plot可视化")
print("=" * 80)

# ============================================================================
# 1. 准备绘图数据
# ============================================================================

def prepare_plot_data(groups, biophysical_features, aa_type):
    """
    准备单个氨基酸类型的绘图数据
    
    参数:
        groups: 包含所有分组数据的字典
        biophysical_features: 要分析的特征列表
        aa_type: 氨基酸类型 ('S', 'T', 'Y')
        
    返回:
        DataFrame: 长格式数据,适合seaborn绘图
    """
    
    plot_data = []
    
    for mobility in ['High', 'Low']:
        key = f"{aa_type}_{mobility}"
        df = groups[key]
        
        for feature in biophysical_features:
            values = df[feature].dropna()
            for val in values:
                plot_data.append({
                    'Feature': feature,
                    'Value': val,
                    'Mobility': mobility,
                    'AA': aa_type
                })
    
    return pd.DataFrame(plot_data)


# ============================================================================
# 2. 绘制按STY拆分的Box Plot
# ============================================================================

def plot_sty_boxplots(groups, biophysical_features, figsize=(20, 5), 
                       save_path=None, dpi=300):
    """
    绘制按S/T/Y拆分的box plot,并添加KS检验结果
    
    参数:
        groups: 包含所有分组数据的字典
        biophysical_features: 要分析的特征列表
        figsize: 图片大小
        save_path: 保存路径 (可选)
        dpi: 分辨率
    """
    
    # 创建子图
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    amino_acids = ['S', 'T', 'Y']
    colors = {'High': '#E74C3C', 'Low': '#3498DB'}  # 红色=高移动性, 蓝色=低移动性
    
    for idx, aa in enumerate(amino_acids):
        ax = axes[idx]
        
        # 准备数据
        plot_df = prepare_plot_data(groups, biophysical_features, aa)
        
        if len(plot_df) == 0:
            ax.text(0.5, 0.5, f'No data for {aa}', 
                   ha='center', va='center', fontsize=14)
            ax.set_title(f'{aa} (Serine)' if aa == 'S' else 
                        f'{aa} (Threonine)' if aa == 'T' else 
                        f'{aa} (Tyrosine)', fontsize=14, fontweight='bold')
            continue
        
        # 绘制box plot
        sns.boxplot(data=plot_df, x='Feature', y='Value', hue='Mobility',
                   palette=colors, ax=ax, showfliers=False, width=0.6)
        
        # 添加数据点 (可选,如果数据量不大)
        # sns.stripplot(data=plot_df, x='Feature', y='Value', hue='Mobility',
        #              palette=colors, ax=ax, dodge=True, alpha=0.3, size=2)
        
        # 设置标题
        aa_name = 'Serine' if aa == 'S' else 'Threonine' if aa == 'T' else 'Tyrosine'
        n_high = len(groups[f'{aa}_High'])
        n_low = len(groups[f'{aa}_Low'])
        ax.set_title(f'{aa} ({aa_name})\nHigh: n={n_high}, Low: n={n_low}', 
                    fontsize=12, fontweight='bold')
        
        # 设置轴标签
        ax.set_xlabel('')
        ax.set_ylabel('Feature Value' if idx == 0 else '', fontsize=11)
        
        # 旋转x轴标签
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
        
        # 调整图例
        if idx == 2:  # 只在最右边显示图例
            ax.legend(title='Mobility', loc='upper right', fontsize=9)
        else:
            ax.get_legend().remove()
        
        # 添加网格
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"✅ 图片已保存: {save_path}")
    
    plt.show()
    
    return fig


# ============================================================================
# 3. 执行绘图
# ============================================================================

print("\n📈 开始绘制按STY拆分的Box Plot...")

# 绘制图表
fig = plot_sty_boxplots(
    groups=groups,
    biophysical_features=biophysical_features,
    figsize=(20, 6),
    save_path='STY_Biophysical_Features_BoxPlot.png',
    dpi=300
)

print("\n" + "=" * 80)
print("✅ Cell 3 完成: Box Plot已生成")
print("=" * 80)


# %%
"""
Cell 4: 添加统计显著性标注
功能: 在box plot上添加KS检验结果和显著性标记
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle

print("=" * 80)
print("📊 Cell 4: 添加统计显著性标注")
print("=" * 80)

# ============================================================================
# 1. 计算所有特征的KS统计量
# ============================================================================

def calculate_all_ks_statistics(groups, biophysical_features):
    """
    计算所有S/T/Y × 特征组合的KS统计量
    
    返回:
        DataFrame: 包含所有KS检验结果
    """
    
    results = []
    
    for aa in ['S', 'T', 'Y']:
        high_key = f'{aa}_High'
        low_key = f'{aa}_Low'
        
        for feature in biophysical_features:
            # 提取数据
            high_data = groups[high_key][feature]
            low_data = groups[low_key][feature]
            
            # 执行KS检验
            ks_result = perform_ks_test_with_resampling(
                high_data, low_data, n_iterations=10
            )
            
            results.append({
                'AA': aa,
                'Feature': feature,
                'KS_statistic': ks_result['ks_statistic'],
                'p_value': ks_result['p_value'],
                'sig_level': ks_result['sig_level'],
                'significant': ks_result['significant'],
                'n_high': ks_result['n_high'],
                'n_low': ks_result['n_low']
            })
    
    return pd.DataFrame(results)


# 计算KS统计量
print("\n🔬 计算KS统计量...")
ks_results_df = calculate_all_ks_statistics(groups, biophysical_features)

print(f"\n📋 KS检验结果汇总:")
print(ks_results_df.to_string(index=False))

# 统计显著性
sig_count = ks_results_df['significant'].sum()
total_count = len(ks_results_df)
print(f"\n📊 显著性统计:")
print(f"  显著的比较: {sig_count} / {total_count} ({sig_count/total_count*100:.1f}%)")


# ============================================================================
# 2. 绘制带统计标注的Box Plot
# ============================================================================

def plot_sty_boxplots_with_stats(groups, biophysical_features, ks_results_df,
                                  figsize=(20, 6), save_path=None, dpi=300):
    """
    绘制带KS统计标注的box plot
    """
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    amino_acids = ['S', 'T', 'Y']
    colors = {'High': '#E74C3C', 'Low': '#3498DB'}
    
    for idx, aa in enumerate(amino_acids):
        ax = axes[idx]
        
        # 准备数据
        plot_df = prepare_plot_data(groups, biophysical_features, aa)
        
        if len(plot_df) == 0:
            continue
        
        # 绘制box plot
        sns.boxplot(data=plot_df, x='Feature', y='Value', hue='Mobility',
                   palette=colors, ax=ax, showfliers=False, width=0.6)
        
        # 获取该氨基酸的KS结果
        aa_ks = ks_results_df[ks_results_df['AA'] == aa]
        
        # 添加显著性标记
        feature_positions = {feat: i for i, feat in enumerate(biophysical_features)}
        y_max = plot_df['Value'].max()
        y_min = plot_df['Value'].min()
        y_range = y_max - y_min
        
        for _, row in aa_ks.iterrows():
            feature = row['Feature']
            sig_level = row['sig_level']
            p_value = row['p_value']
            
            if feature not in feature_positions:
                continue
            
            x_pos = feature_positions[feature]
            
            # 计算标注位置
            y_pos = y_max + y_range * 0.05
            
            # 添加显著性星号
            if sig_level != 'ns':
                ax.text(x_pos, y_pos, sig_level, 
                       ha='center', va='bottom', fontsize=12, 
                       fontweight='bold', color='black')
            
            # 添加p值 (小字)
            if not np.isnan(p_value):
                p_text = f'p={p_value:.3f}' if p_value >= 0.001 else 'p<0.001'
                ax.text(x_pos, y_pos + y_range * 0.08, p_text,
                       ha='center', va='bottom', fontsize=7, 
                       color='gray', style='italic')
        
        # 设置标题
        aa_name = 'Serine' if aa == 'S' else 'Threonine' if aa == 'T' else 'Tyrosine'
        n_high = len(groups[f'{aa}_High'])
        n_low = len(groups[f'{aa}_Low'])
        ax.set_title(f'{aa} ({aa_name})\nHigh: n={n_high}, Low: n={n_low}', 
                    fontsize=12, fontweight='bold')
        
        # 设置轴标签
        ax.set_xlabel('')
        ax.set_ylabel('Feature Value' if idx == 0 else '', fontsize=11)
        
        # 旋转x轴标签
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
        
        # 调整y轴范围,为标注留空间
        ax.set_ylim(y_min - y_range * 0.05, y_max + y_range * 0.15)
        
        # 调整图例
        if idx == 2:
            ax.legend(title='Mobility', loc='upper right', fontsize=9)
        else:
            ax.get_legend().remove()
        
        # 添加网格
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
    
    # 添加总标题
    fig.suptitle('Biophysical Features Comparison: High vs Low Mobility (by Amino Acid Type)', 
                fontsize=14, fontweight='bold', y=1.00)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"✅ 图片已保存: {save_path}")
    
    plt.show()
    
    return fig


# ============================================================================
# 3. 执行绘图
# ============================================================================

print("\n📈 绘制带统计标注的Box Plot...")

fig = plot_sty_boxplots_with_stats(
    groups=groups,
    biophysical_features=biophysical_features,
    ks_results_df=ks_results_df,
    figsize=(20, 7),
    save_path='STY_Biophysical_Features_BoxPlot_with_Stats.png',
    dpi=300
)

# ============================================================================
# 4. 保存KS统计结果到CSV
# ============================================================================

ks_results_df.to_csv('KS_Statistics_Results.csv', index=False)
print(f"\n💾 KS统计结果已保存: KS_Statistics_Results.csv")

print("\n" + "=" * 80)
print("✅ Cell 4 完成: 统计标注已添加")
print("=" * 80)


# %%
"""
Cell 5: Feature Map可视化 (最终版)
功能: 展示top蛋白质的磷酸化位点在序列上的特征分布
更新: 直接排除CD3EAP(原始数据缺少Position信息)
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

print("=" * 80)
print("📊 Cell 5: Feature Map可视化")
print("=" * 80)

# ============================================================================
# 1. 计算每个蛋白质(基因)的high/low ratio
# ============================================================================

print("\n🔬 计算蛋白质mobility ratio...")

def calculate_protein_mobility_ratio(df_valid):
    """计算每个蛋白质的high/low mobility ratio"""
    
    protein_stats = []
    
    for gene in df_valid['Gene'].unique():
        gene_df = df_valid[df_valid['Gene'] == gene]
        
        n_total = len(gene_df)
        n_high = (gene_df['mobility_group'] == 'High').sum()
        n_low = (gene_df['mobility_group'] == 'Low').sum()
        
        # 计算ratio (避免除以0)
        if n_low > 0:
            ratio = n_high / n_low
        else:
            ratio = n_high if n_high > 0 else 0
        
        protein_stats.append({
            'Gene': gene,
            'n_total': n_total,
            'n_high': n_high,
            'n_low': n_low,
            'high_low_ratio': ratio,
            'has_high_mobility': n_high > 0
        })
    
    return pd.DataFrame(protein_stats)


protein_stats_df = calculate_protein_mobility_ratio(df_valid)

# 只保留有高移动性位点的蛋白质
protein_stats_df = protein_stats_df[protein_stats_df['has_high_mobility']].copy()

# 直接排除CD3EAP (原始数据缺少Position信息)
protein_stats_df = protein_stats_df[protein_stats_df['Gene'] != 'CD3EAP'].copy()

# 按ratio排序
protein_stats_df = protein_stats_df.sort_values('high_low_ratio', ascending=False)

print(f"\n📋 Top 20 蛋白质 (按high/low ratio排序, 已排除CD3EAP):")
print(protein_stats_df.head(20).to_string(index=False))


# ============================================================================
# 2. 提取Top蛋白质的位点数据
# ============================================================================

print("\n🎯 提取Top蛋白质的位点数据...")

def extract_top_proteins_sites(df_valid, top_genes, biophysical_features, 
                                 mobility_group='High'):
    """提取top蛋白质的高移动性位点数据"""
    
    # 筛选数据
    top_df = df_valid[
        (df_valid['Gene'].isin(top_genes)) & 
        (df_valid['mobility_group'] == mobility_group)
    ].copy()
    
    # 确保有Position列
    if 'Position' not in top_df.columns and 'sequence_position' in top_df.columns:
        top_df['Position'] = top_df['sequence_position']
    
    # 选择需要的列
    cols_to_keep = ['Gene', 'Site', 'Position', 'AA'] + biophysical_features
    top_df = top_df[cols_to_keep].copy()
    
    # 按基因和位置排序
    top_df = top_df.sort_values(['Gene', 'Position'])
    
    return top_df


# 选择top 10个蛋白质
n_top = 10
top_genes = protein_stats_df.head(n_top)['Gene'].tolist()

print(f"\n🎯 选择的Top {n_top}个蛋白质:")
print(f"\n挑选原则: 按high/low mobility ratio从高到低排序,选择高移动性位点相对富集程度最高的蛋白质。\n")
for i, gene in enumerate(top_genes, 1):
    ratio = protein_stats_df[protein_stats_df['Gene'] == gene]['high_low_ratio'].values[0]
    n_high = protein_stats_df[protein_stats_df['Gene'] == gene]['n_high'].values[0]
    n_total = protein_stats_df[protein_stats_df['Gene'] == gene]['n_total'].values[0]
    print(f"  {i:2d}. {gene:15s} - Ratio: {ratio:.2f}, High: {n_high}/{n_total}")

# 提取高移动性组的位点数据
top_sites_df = extract_top_proteins_sites(
    df_valid, top_genes, biophysical_features, mobility_group='High'
)

print(f"\n📊 提取的位点数据: {len(top_sites_df)} 个高移动性位点")


# ============================================================================
# 3. 绘制Feature Map - 方案A: 每个蛋白质一个subplot
# ============================================================================

def plot_feature_map_subplots(top_sites_df, biophysical_features, 
                               figsize=(20, 12), save_path=None, dpi=300):
    """绘制Feature Map: 每个蛋白质一个subplot"""
    
    genes = top_sites_df['Gene'].unique()
    n_genes = len(genes)
    
    # 计算subplot布局
    n_cols = 2
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    # 为每个特征分配颜色
    colors = plt.cm.tab10(np.linspace(0, 1, len(biophysical_features)))
    feature_colors = dict(zip(biophysical_features, colors))
    
    for idx, gene in enumerate(genes):
        ax = axes[idx]
        
        # 提取该蛋白质的数据
        gene_df = top_sites_df[top_sites_df['Gene'] == gene].copy()
        
        # 绘制每个特征
        for feature in biophysical_features:
            positions = gene_df['Position'].values
            values = gene_df[feature].values
            
            # 绘制线图+散点
            ax.plot(positions, values, 'o-', 
                   color=feature_colors[feature], 
                   label=feature, 
                   markersize=6, 
                   linewidth=1.5,
                   alpha=0.7)
        
        # 添加零线
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        
        # 设置标题和标签
        n_sites = len(gene_df)
        ax.set_title(f'{gene} ({n_sites} sites)', fontsize=11, fontweight='bold')
        ax.set_xlabel('Sequence Position', fontsize=10)
        ax.set_ylabel('Feature Value', fontsize=10)
        
        # 添加网格
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        # 只在第一个subplot显示图例
        if idx == 0:
            ax.legend(loc='upper left', fontsize=8, ncol=2)
    
    # 隐藏多余的subplot
    for idx in range(n_genes, len(axes)):
        axes[idx].axis('off')
    
    # 总标题
    fig.suptitle(f'Feature Map: Top {n_genes} Proteins (High Mobility Sites)', 
                fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"✅ Feature Map已保存: {save_path}")
    
    plt.show()
    
    return fig


# ============================================================================
# 4. 绘制Feature Map - 方案B: 所有蛋白质叠加
# ============================================================================

def plot_feature_map_overlay(top_sites_df, biophysical_features, 
                              figsize=(20, 10), save_path=None, dpi=300):
    """绘制Feature Map: 所有蛋白质叠加在一起"""
    
    n_features = len(biophysical_features)
    n_cols = 3
    n_rows = (n_features + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()
    
    genes = top_sites_df['Gene'].unique()
    gene_colors = plt.cm.tab20(np.linspace(0, 1, len(genes)))
    gene_color_map = dict(zip(genes, gene_colors))
    
    for idx, feature in enumerate(biophysical_features):
        ax = axes[idx]
        
        # 绘制每个蛋白质
        for gene in genes:
            gene_df = top_sites_df[top_sites_df['Gene'] == gene].copy()
            
            positions = gene_df['Position'].values
            values = gene_df[feature].values
            
            ax.plot(positions, values, 'o-', 
                   color=gene_color_map[gene], 
                   label=gene if idx == 0 else '', 
                   markersize=4, 
                   linewidth=1,
                   alpha=0.6)
        
        # 添加零线
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        
        # 设置标题和标签
        ax.set_title(feature, fontsize=11, fontweight='bold')
        ax.set_xlabel('Sequence Position', fontsize=10)
        ax.set_ylabel('Feature Value', fontsize=10)
        
        # 添加网格
        ax.grid(axis='both', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
    
    # 隐藏多余的subplot
    for idx in range(n_features, len(axes)):
        axes[idx].axis('off')
    
    # 添加图例(在第一个subplot)
    if len(genes) <= 10:
        axes[0].legend(loc='upper left', fontsize=7, ncol=2)
    
    # 总标题
    fig.suptitle(f'Feature Map: All Proteins Overlay (High Mobility Sites)', 
                fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"✅ Feature Map (Overlay)已保存: {save_path}")
    
    plt.show()
    
    return fig


# ============================================================================
# 5. 执行绘图
# ============================================================================

print("\n📈 绘制Feature Map (方案A: 每个蛋白质一个subplot)...")
fig1 = plot_feature_map_subplots(
    top_sites_df=top_sites_df,
    biophysical_features=biophysical_features,
    figsize=(20, 15),
    save_path='Feature_Map_Subplots.png',
    dpi=300
)

print("\n📈 绘制Feature Map (方案B: 所有蛋白质叠加)...")
fig2 = plot_feature_map_overlay(
    top_sites_df=top_sites_df,
    biophysical_features=biophysical_features,
    figsize=(18, 12),
    save_path='Feature_Map_Overlay.png',
    dpi=300
)


# ============================================================================
# 6. 保存数据
# ============================================================================

protein_stats_df.to_csv('Protein_Mobility_Statistics.csv', index=False)
print(f"\n💾 蛋白质统计已保存: Protein_Mobility_Statistics.csv")

top_sites_df.to_csv('Top_Proteins_Sites_Features.csv', index=False)
print(f"💾 Top蛋白质位点特征已保存: Top_Proteins_Sites_Features.csv")

print("\n" + "=" * 80)
print("✅ Cell 5 完成: Feature Map已生成")
print("   方案A: 每个蛋白质一个subplot - Feature_Map_Subplots.png")
print("   方案B: 所有蛋白质叠加 - Feature_Map_Overlay.png")
print("=" * 80)


# %%
# ===== Cell 26a: 方案A - 按Movement Score梯度细分 (修正版) =====

import pandas as pd
import numpy as np
from scipy import stats

print("="*80)
print("方案A: 按Movement Score梯度细分")
print("="*80)

# ✅ 修正: 使用final_df而不是success_df (final_df包含完整的7个特征)
# 确保使用的数据源与Cell 1一致
if 'final_df' not in globals():
    print("\n❌ 错误: final_df未定义!")
    print("请先运行Cell 1生成final_df")
    raise NameError("final_df未定义,请先运行Cell 1")

# 使用Cell 1中定义的movement_score
if 'movement_score' not in final_df.columns:
    print("\n❌ 错误: final_df中没有movement_score列!")
    print("请先运行Cell 1计算movement_score")
    raise KeyError("movement_score列不存在,请先运行Cell 1")

# [1/5] 定义细分组
print("\n[1/5] 定义细分组...")

# High mobility组细分 (使用百分比阈值)
very_high = final_df[final_df['movement_score'] >= 20].copy()
high = final_df[(final_df['movement_score'] >= 15) & (final_df['movement_score'] < 20)].copy()
moderate_high = final_df[(final_df['movement_score'] >= 10) & (final_df['movement_score'] < 15)].copy()

# Low mobility组细分
very_low = final_df[final_df['movement_score'] < 2].copy()
low = final_df[(final_df['movement_score'] >= 2) & (final_df['movement_score'] < 5)].copy()

print(f"\nHigh mobility组细分:")
print(f"  Very High (≥20%):        n={len(very_high):5,}")
print(f"  High (15-20%):           n={len(high):5,}")
print(f"  Moderate High (10-15%):  n={len(moderate_high):5,}")
print(f"  Total High (≥10%):       n={len(very_high) + len(high) + len(moderate_high):5,}")

print(f"\nLow mobility组细分:")
print(f"  Very Low (<2%):          n={len(very_low):5,}")
print(f"  Low (2-5%):              n={len(low):5,}")
print(f"  Total Low (<5%):         n={len(very_low) + len(low):5,}")

# [2/5] 检查生物物理特征列
print("\n[2/5] 检查生物物理特征列...")

# 定义期望的7个特征
expected_features = [
    'backbone_dynamics',
    'sidechain_dynamics',
    'disorder_propensity',
    'helix_propensity',
    'sheet_propensity',
    'coil_propensity',
    'earlyFolding'
]

# 检查哪些特征存在
available_features = [f for f in expected_features if f in final_df.columns]
missing_features = [f for f in expected_features if f not in final_df.columns]

print(f"\n可用特征 ({len(available_features)}个):")
for f in available_features:
    print(f"  ✓ {f}")

if missing_features:
    print(f"\n缺失特征 ({len(missing_features)}个):")
    for f in missing_features:
        print(f"  ✗ {f}")
    print("\n⚠️  将只分析可用特征")

bio_features = available_features

# [3/5] 计算各组的特征统计量
print("\n[3/5] 计算各组的特征统计量...")

groups = {
    'Very High': very_high,
    'High': high,
    'Moderate High': moderate_high,
    'Very Low': very_low,
    'Low': low
}

stats_summary = {}
for group_name, group_df in groups.items():
    stats_summary[group_name] = {
        'n': len(group_df),
        'mean': group_df[bio_features].mean(),
        'median': group_df[bio_features].median(),
        'std': group_df[bio_features].std()
    }

print("✓ 统计量计算完成")

# [4/5] 梯度效应检验 (High mobility组内部比较)
print("\n[4/5] 梯度效应检验 (High mobility组内部比较)...")

print("\n" + "="*80)
print("检验问题: 在High mobility组内部,Movement Score越高,特征值是否呈梯度变化?")
print("="*80)

print("\n梯度效应检验结果:")
print("-"*80)

gradient_results = []
for feature in bio_features:
    vh_data = very_high[feature].dropna()
    h_data = high[feature].dropna()
    mh_data = moderate_high[feature].dropna()
    
    if len(vh_data) == 0 or len(h_data) == 0 or len(mh_data) == 0:
        continue
    
    # Kruskal-Wallis检验 (3组比较)
    kw_stat, kw_p = stats.kruskal(vh_data, h_data, mh_data)
    
    # 计算均值
    vh_mean = vh_data.mean()
    h_mean = h_data.mean()
    mh_mean = mh_data.mean()
    
    # 判断趋势
    if vh_mean > h_mean > mh_mean:
        trend = "Decreasing"
    elif vh_mean < h_mean < mh_mean:
        trend = "Increasing"
    else:
        trend = "Non-monotonic"
    
    # Very High vs Moderate High直接比较
    vh_mh_stat, vh_mh_p = stats.mannwhitneyu(vh_data, mh_data, alternative='two-sided')
    
    print(f"\n{feature}:")
    print(f"  Very High mean:     {vh_mean:.4f}")
    print(f"  High mean:          {h_mean:.4f}")
    print(f"  Moderate High mean: {mh_mean:.4f}")
    print(f"  Trend:              {trend}")
    print(f"  Kruskal-Wallis p:   {kw_p:.4e}  {'***' if kw_p < 0.001 else '**' if kw_p < 0.01 else '*' if kw_p < 0.05 else 'ns'}")
    print(f"  VH vs MH p:         {vh_mh_p:.4e}")
    
    gradient_results.append({
        'Feature': feature,
        'VH_mean': vh_mean,
        'H_mean': h_mean,
        'MH_mean': mh_mean,
        'Trend': trend,
        'KW_p': kw_p,
        'VH_MH_p': vh_mh_p
    })

# [5/5] High vs Low对比 (各档分别比较)
print("\n[5/5] High vs Low对比 (各档分别比较)...")

print("\n" + "="*80)
print("检验问题: 各个High组与Low组相比,特征差异是否显著?")
print("="*80)

comparison_results = []
high_groups = [('Very High', very_high), ('High', high), ('Moderate High', moderate_high)]
low_groups = [('Very Low', very_low), ('Low', low)]

for high_name, high_df in high_groups:
    for low_name, low_df in low_groups:
        for feature in bio_features:
            high_data = high_df[feature].dropna()
            low_data = low_df[feature].dropna()
            
            if len(high_data) == 0 or len(low_data) == 0:
                continue
            
            # Mann-Whitney U检验
            stat, p = stats.mannwhitneyu(high_data, low_data, alternative='two-sided')
            
            # 效应量 (Cohen's d)
            mean_diff = high_data.mean() - low_data.mean()
            pooled_std = np.sqrt((high_data.std()**2 + low_data.std()**2) / 2)
            cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
            
            comparison_results.append({
                'High_group': high_name,
                'Low_group': low_name,
                'Feature': feature,
                'Mean_diff': mean_diff,
                'Cohens_d': cohens_d,
                'p_value': p,
                'Significant': p < 0.05
            })

# 只显示显著的结果
significant_results = [r for r in comparison_results if r['Significant']]
print(f"\n显著差异的比较 (p < 0.05): {len(significant_results)} 个")
print("-"*80)

for high_name, _ in high_groups:
    for low_name, _ in low_groups:
        comparison_name = f"{high_name} vs {low_name}"
        relevant_results = [r for r in significant_results 
                          if r['High_group'] == high_name and r['Low_group'] == low_name]
        
        if relevant_results:
            print(f"\n{comparison_name}:")
            for r in relevant_results:
                sign = '+' if r['Mean_diff'] > 0 else ''
                sig_marker = '***' if r['p_value'] < 0.001 else '**' if r['p_value'] < 0.01 else '*'
                print(f"  {r['Feature']:25s}: Δ={sign}{r['Mean_diff']:.4f}, "
                      f"d={sign}{r['Cohens_d']:.3f}, p={r['p_value']:.4e} {sig_marker}")

# 总结
print("\n" + "="*80)
print("方案A总结")
print("="*80)

print(f"\n1. 样本量分布:")
print(f"   - Very High组: {len(very_high)} {'(可能偏小)' if len(very_high) < 100 else ''}")
print(f"   - High组: {len(high)}")
print(f"   - Moderate High组: {len(moderate_high)}")

monotonic_trends = sum(1 for r in gradient_results if r['Trend'] in ['Increasing', 'Decreasing'])
significant_gradients = sum(1 for r in gradient_results if r['KW_p'] < 0.05)

print(f"\n2. 梯度效应:")
print(f"   - 呈单调趋势的特征: {monotonic_trends}/{len(bio_features)}")
print(f"   - 梯度显著的特征 (p<0.05): {significant_gradients}/{len(bio_features)}")

for high_name, _ in high_groups:
    for low_name, _ in low_groups:
        count = sum(1 for r in significant_results 
                   if r['High_group'] == high_name and r['Low_group'] == low_name)
        if count > 0:
            print(f"\n3. High vs Low差异:")
            print(f"   - {high_name} vs {low_name}显著特征: {count}/{len(bio_features)}")

print("\n" + "="*80)
print("✓ 方案A数据分析完成!")
print("="*80)

print("\n请将以上结果复制给AI进行评价。")


# %%
# ===== Cell 26b: 方案B - 按响应时间细分的数据分析 =====
# 只做统计分析,不生成图表

import pandas as pd
import numpy as np
from scipy import stats

print("=" * 80)
print("方案B: 按响应时间细分")
print("=" * 80)

# ========== 1. 定义细分组 ==========
print("\n[1/5] 定义细分组...")

# High mobility组按响应时间细分为2档
early_responders = success_df[
    (success_df['movement_score'] >= 10) & 
    (success_df['movement_best_time'].isin(['2min', '8min']))
].copy()

late_responders = success_df[
    (success_df['movement_score'] >= 10) & 
    (success_df['movement_best_time'].isin(['20min', '90min']))
].copy()

# Low mobility组保持不变
low_mobility = success_df[success_df['movement_score'] < 5].copy()

print(f"\nHigh mobility组按响应时间细分:")
print(f"  Early responders (2/8min):   n={len(early_responders):5,}")
print(f"  Late responders (20/90min):  n={len(late_responders):5,}")
print(f"  Total High (≥10%):           n={len(early_responders) + len(late_responders):5,}")

print(f"\nLow mobility组:")
print(f"  Low mobility (<5%):          n={len(low_mobility):5,}")

# 进一步细分Early和Late的时间点
print(f"\nEarly responders细分:")
early_2min = early_responders[early_responders['movement_best_time'] == '2min']
early_8min = early_responders[early_responders['movement_best_time'] == '8min']
print(f"  2min:  n={len(early_2min):5,}")
print(f"  8min:  n={len(early_8min):5,}")

print(f"\nLate responders细分:")
late_20min = late_responders[late_responders['movement_best_time'] == '20min']
late_90min = late_responders[late_responders['movement_best_time'] == '90min']
print(f"  20min: n={len(late_20min):5,}")
print(f"  90min: n={len(late_90min):5,}")

# ========== 2. 检查生物物理特征列是否存在 ==========
print("\n[2/5] 检查生物物理特征列...")

required_features = [
    'backbone_dynamics',
    'sidechain_dynamics',
    'disorder_propensity',
    'helix_propensity',
    'sheet_propensity',
    'coil_propensity',
    'earlyFolding'
]

available_features = [f for f in required_features if f in success_df.columns]
missing_features = [f for f in required_features if f not in success_df.columns]

print(f"\n可用特征 ({len(available_features)}个):")
for f in available_features:
    print(f"  ✓ {f}")

if missing_features:
    print(f"\n缺失特征 ({len(missing_features)}个):")
    for f in missing_features:
        print(f"  ✗ {f}")
    print("\n⚠️  将只分析可用特征")

# ========== 3. 计算各组的特征均值和标准差 ==========
print("\n[3/5] 计算各组的特征统计量...")

groups = {
    'Early responders (2/8min)': early_responders,
    'Late responders (20/90min)': late_responders,
    'Low mobility (<5%)': low_mobility
}

stats_summary = []

for group_name, group_df in groups.items():
    for feature in available_features:
        # 去除NaN值
        values = group_df[feature].dropna()
        
        if len(values) > 0:
            stats_summary.append({
                'Group': group_name,
                'Feature': feature,
                'N': len(values),
                'Mean': values.mean(),
                'Std': values.std(),
                'Median': values.median(),
                'Q25': values.quantile(0.25),
                'Q75': values.quantile(0.75)
            })

stats_df = pd.DataFrame(stats_summary)

print("\n✓ 统计量计算完成")

# ========== 4. Early vs Late比较 ==========
print("\n[4/5] Early vs Late响应者比较...")
print("\n" + "=" * 80)
print("检验问题: Early和Late响应者的生物物理特征是否有差异?")
print("=" * 80)

early_late_results = []

for feature in available_features:
    early_values = early_responders[feature].dropna()
    late_values = late_responders[feature].dropna()
    
    if len(early_values) > 0 and len(late_values) > 0:
        # Mann-Whitney U检验
        u_stat, p_value = stats.mannwhitneyu(early_values, late_values, alternative='two-sided')
        
        # 效应量 (Cohen's d)
        mean_diff = early_values.mean() - late_values.mean()
        pooled_std = np.sqrt((early_values.std()**2 + late_values.std()**2) / 2)
        cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
        
        early_late_results.append({
            'Feature': feature,
            'Early_mean': early_values.mean(),
            'Late_mean': late_values.mean(),
            'Mean_diff': mean_diff,
            'Cohens_d': cohens_d,
            'p_value': p_value,
            'Significant': '***' if p_value < 0.001 else ('**' if p_value < 0.01 else ('*' if p_value < 0.05 else 'ns'))
        })

early_late_df = pd.DataFrame(early_late_results)

print("\nEarly vs Late响应者比较结果:")
print("-" * 80)
for idx, row in early_late_df.iterrows():
    print(f"\n{row['Feature']}:")
    print(f"  Early mean:  {row['Early_mean']:.4f}")
    print(f"  Late mean:   {row['Late_mean']:.4f}")
    print(f"  Difference:  {row['Mean_diff']:+.4f}")
    print(f"  Cohen's d:   {row['Cohens_d']:+.3f}")
    print(f"  p-value:     {row['p_value']:.4e}  {row['Significant']}")

# ========== 5. 各组与Low mobility比较 ==========
print("\n[5/5] 各组与Low mobility比较...")
print("\n" + "=" * 80)
print("检验问题: Early和Late响应者分别与Low mobility组相比,特征差异是否显著?")
print("=" * 80)

comparison_results = []

comparisons = {
    'Early vs Low': (early_responders, low_mobility),
    'Late vs Low': (late_responders, low_mobility)
}

for comparison_name, (group1, group2) in comparisons.items():
    for feature in available_features:
        g1_values = group1[feature].dropna()
        g2_values = group2[feature].dropna()
        
        if len(g1_values) > 0 and len(g2_values) > 0:
            # Mann-Whitney U检验
            u_stat, p_value = stats.mannwhitneyu(g1_values, g2_values, alternative='two-sided')
            
            # 效应量 (Cohen's d)
            mean_diff = g1_values.mean() - g2_values.mean()
            pooled_std = np.sqrt((g1_values.std()**2 + g2_values.std()**2) / 2)
            cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
            
            comparison_results.append({
                'Comparison': comparison_name,
                'Feature': feature,
                'Group1_mean': g1_values.mean(),
                'Group2_mean': g2_values.mean(),
                'Mean_diff': mean_diff,
                'Cohens_d': cohens_d,
                'p_value': p_value,
                'Significant': '***' if p_value < 0.001 else ('**' if p_value < 0.01 else ('*' if p_value < 0.05 else 'ns'))
            })

comparison_df = pd.DataFrame(comparison_results)

for comparison in comparisons.keys():
    subset = comparison_df[comparison_df['Comparison'] == comparison]
    print(f"\n{comparison}:")
    print("-" * 80)
    for idx, row in subset.iterrows():
        print(f"  {row['Feature']:25s}: Δ={row['Mean_diff']:+.4f}, d={row['Cohens_d']:+.3f}, p={row['p_value']:.4e} {row['Significant']}")

# ========== 6. Movement Score分布比较 ==========
print("\n" + "=" * 80)
print("Movement Score分布比较")
print("=" * 80)

print(f"\nEarly responders的Movement Score:")
print(f"  Mean:   {early_responders['movement_score'].mean():.2f}%")
print(f"  Median: {early_responders['movement_score'].median():.2f}%")
print(f"  Std:    {early_responders['movement_score'].std():.2f}%")
print(f"  Range:  {early_responders['movement_score'].min():.2f}% - {early_responders['movement_score'].max():.2f}%")

print(f"\nLate responders的Movement Score:")
print(f"  Mean:   {late_responders['movement_score'].mean():.2f}%")
print(f"  Median: {late_responders['movement_score'].median():.2f}%")
print(f"  Std:    {late_responders['movement_score'].std():.2f}%")
print(f"  Range:  {late_responders['movement_score'].min():.2f}% - {late_responders['movement_score'].max():.2f}%")

# 检验Movement Score是否有差异
u_stat, p_value = stats.mannwhitneyu(
    early_responders['movement_score'].dropna(), 
    late_responders['movement_score'].dropna(), 
    alternative='two-sided'
)
print(f"\nMovement Score差异检验 (Early vs Late):")
print(f"  Mann-Whitney U p-value: {p_value:.4e}")
print(f"  结论: {'Early和Late的Movement Score有显著差异' if p_value < 0.05 else 'Early和Late的Movement Score无显著差异'}")

# ========== 7. 氨基酸分布比较 ==========
print("\n" + "=" * 80)
print("氨基酸分布比较")
print("=" * 80)

print(f"\nEarly responders的氨基酸分布:")
early_aa = early_responders['AA'].value_counts()
for aa, count in early_aa.items():
    pct = count / len(early_responders) * 100
    print(f"  {aa}: {count:4d} ({pct:5.1f}%)")

print(f"\nLate responders的氨基酸分布:")
late_aa = late_responders['AA'].value_counts()
for aa, count in late_aa.items():
    pct = count / len(late_responders) * 100
    print(f"  {aa}: {count:4d} ({pct:5.1f}%)")

# Chi-square检验氨基酸分布是否有差异
contingency_table = pd.crosstab(
    pd.concat([early_responders, late_responders])['AA'],
    pd.concat([
        pd.Series(['Early'] * len(early_responders)),
        pd.Series(['Late'] * len(late_responders))
    ])
)
chi2, p_chi, dof, expected = stats.chi2_contingency(contingency_table)
print(f"\n氨基酸分布差异检验 (Chi-square):")
print(f"  p-value: {p_chi:.4e}")
print(f"  结论: {'Early和Late的氨基酸分布有显著差异' if p_chi < 0.05 else 'Early和Late的氨基酸分布无显著差异'}")

# ========== 8. 总结 ==========
print("\n" + "=" * 80)
print("方案B总结")
print("=" * 80)

print(f"\n1. 样本量分布:")
print(f"   - Early responders: {len(early_responders)}")
print(f"   - Late responders: {len(late_responders)}")
print(f"   - 样本量比例: {len(early_responders)/len(late_responders):.2f}:1")

print(f"\n2. Early vs Late差异:")
significant_early_late = early_late_df[early_late_df['p_value'] < 0.05]
print(f"   - 显著差异的特征: {len(significant_early_late)}/{len(available_features)}")
if len(significant_early_late) > 0:
    print(f"   - 显著特征列表:")
    for idx, row in significant_early_late.iterrows():
        direction = '↑' if row['Mean_diff'] > 0 else '↓'
        print(f"     • {row['Feature']:25s} (Early {direction} Late, d={row['Cohens_d']:+.3f}, p={row['p_value']:.4e})")

print(f"\n3. 与Low mobility比较:")
early_vs_low_sig = comparison_df[(comparison_df['Comparison'] == 'Early vs Low') & (comparison_df['p_value'] < 0.05)]
late_vs_low_sig = comparison_df[(comparison_df['Comparison'] == 'Late vs Low') & (comparison_df['p_value'] < 0.05)]
print(f"   - Early vs Low显著特征: {len(early_vs_low_sig)}/{len(available_features)}")
print(f"   - Late vs Low显著特征: {len(late_vs_low_sig)}/{len(available_features)}")

print(f"\n4. Movement Score:")
print(f"   - Early和Late的Movement Score {'有显著差异' if p_value < 0.05 else '无显著差异'} (p={p_value:.4e})")

print(f"\n5. 氨基酸分布:")
print(f"   - Early和Late的氨基酸分布 {'有显著差异' if p_chi < 0.05 else '无显著差异'} (p={p_chi:.4e})")

print("\n" + "=" * 80)
print("✓ 方案B数据分析完成!")
print("=" * 80)
print("\n请将以上结果复制给AI进行评价。")



