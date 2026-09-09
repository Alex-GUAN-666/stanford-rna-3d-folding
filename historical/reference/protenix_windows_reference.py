# MODEL_TYPE 用于指定模型类型
MODEL_TYPE = 'protenix'

# VALIDATION 用于控制是否进行验证步骤
VALIDATION = False

# 导入相关的 Python 包
import os
import sys
import re
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio.PDB import PDBParser
import time

# 获取当前运行的 Python 解释器的路径
PYTHON = sys.executable
print('PYTHON', PYTHON)

# 设置本地路径
RHONET_DIR = 'C:/StanfordRNA3D/RhoFold-main'  # RhoFold 代码库路径
USALIGN = 'C:/StanfordRNA3D/tools/USalign'    # USalign 工具路径
DATA_LOCAL_DIR = 'C:/StanfordRNA3D/data'      # 数据集路径

# 确保 USalign 可执行文件具有执行权限（Windows 不需要 chmod，但确保文件存在）
if not os.path.exists(USALIGN):
    raise FileNotFoundError(f"USalign not found at {USALIGN}. Please ensure the file exists.")

# 添加 RhoFold 路径到系统路径
sys.path.append(RHONET_DIR)

# 定义 dotdict 类，允许通过点号访问字典键
class dotdict(dict):
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

# 设置 3D 图形等比例缩放
def set_aspect_equal(ax):
    x_limits = ax.get_xlim()
    y_limits = ax.get_ylim()
    z_limits = ax.get_zlim()

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    max_range = max(x_limits[1] - x_limits[0],
                    y_limits[1] - y_limits[0],
                    z_limits[1] - z_limits[0]) / 2.0

    ax.set_xlim(x_middle - max_range, x_middle + max_range)
    ax.set_ylim(y_middle - max_range, y_middle + max_range)
    ax.set_zlim(z_middle - max_range, z_middle + max_range)

# 从标签数据集中提取特定目标的真实数据
def get_truth_df(target_id):
    truth_df = LABEL_DF[LABEL_DF['target_id'] == target_id]
    truth_df = truth_df.reset_index(drop=True)
    return truth_df

# 将模型输出解析为 DataFrame
def parse_output_to_df(output, seq, target_id):
    df = []
    chain_data = []
    for i, res in enumerate(seq):
        d = dict(ID=target_id, resname=res, resid=i+1)
        for n in range(len(output)):
            d = {**d, f'x_{n+1}': round(output[n,i,0].item(),3),
                     f'y_{n+1}': round(output[n,i,1].item(),3),
                     f'z_{n+1}': round(output[n,i,2].item(),3)}
        chain_data.append(d)

    if len(chain_data) != 0:
        chain_df = pd.DataFrame(chain_data)
        df.append(chain_df)
    return df

# 从 PDB 文件解析 C1' 原子坐标
def parse_pdb_to_df(pdb_file, target_id):
    parser = PDBParser()
    structure = parser.get_structure('', pdb_file)

    df = []
    for model in structure:
        for chain in model:
            print(chain)
            chain_data = []
            for residue in chain:
                if residue.get_resname() in ['A', 'U', 'G', 'C']:
                    if "C1'" in residue:
                        atom = residue["C1'"]
                        xyz = atom.get_coord()
                        resname = residue.get_resname()
                        resid = residue.get_id()[1]

                        chain_data.append(dict(
                            ID=target_id + '_' + str(resid),
                            resname=resname,
                            resid=resid,
                            x_1=xyz[0],
                            y_1=xyz[1],
                            z_1=xyz[2],
                        ))

            if len(chain_data) != 0:
                chain_df = pd.DataFrame(chain_data)
                df.append(chain_df)
    return df

# 生成 PDB 格式的单行记录
def write_target_line(atom_name, atom_serial, residue_name, chain_id, residue_num, x_coord, y_coord, z_coord, occupancy=1.0, b_factor=0.0, atom_type='P'):
    return f'ATOM  {atom_serial:>5d}  {atom_name:<5s} {residue_name:<3s} {residue_num:>3d}    {x_coord:>8.3f}{y_coord:>8.3f}{z_coord:>8.3f}{occupancy:>6.2f}{b_factor:>6.2f}           {atom_type}\n'

# 将坐标数据写入 PDB 文件
def write_xyz_to_pdb(df, pdb_file, xyz_id=1):
    resolved_cnt = 0
    with open(pdb_file, 'w') as target_file:
        for _, row in df.iterrows():
            x_coord = row[f'x_{xyz_id}']
            y_coord = row[f'y_{xyz_id}']
            z_coord = row[f'z_{xyz_id}']
            if x_coord > -1e17 and y_coord > -1e17 and z_coord > -1e17:
                resolved_cnt += 1
                target_line = write_target_line(
                    atom_name="C1'",
                    atom_serial=int(row['resid']),
                    residue_name=row['resname'],
                    chain_id='0',
                    residue_num=int(row['resid']),
                    x_coord=x_coord,
                    y_coord=y_coord,
                    z_coord=z_coord,
                    atom_type='C',
                )
                target_file.write(target_line)
    return resolved_cnt

# 解析 USalign 输出以提取 TM-score
def parse_usalign_for_tm_score(output):
    tm_score_match = re.findall(r'TM-score=\s+([\d.]+)', output)[1]
    if not tm_score_match:
        raise ValueError('No TM score found')
    return float(tm_score_match)

# 解析 USalign 输出以提取变换矩阵
def parse_usalign_for_transform(output):
    matrix_lines = []
    found_matrix = False
    for line in output.splitlines():
        if "The rotation matrix to rotate Structure_1 to Structure_2" in line:
            found_matrix = True
        elif found_matrix and re.match(r'^\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+$', line):
            matrix_lines.append(line)
        elif found_matrix and not line.strip():
            break
    rotation_matrix = []
    for line in matrix_lines:
        parts = line.split()
        row_values = list(map(float, parts[1:]))
        rotation_matrix.append(row_values)
    return np.array(rotation_matrix)

# 调用 USalign 进行结构比对
def call_usalign(predict_df, truth_df, verbose=1):
    truth_pdb = 'truth.pdb'
    predict_pdb = 'predict.pdb'
    write_xyz_to_pdb(predict_df, predict_pdb, xyz_id=1)
    write_xyz_to_pdb(truth_df, truth_pdb, xyz_id=1)
    command = f'{USALIGN} {predict_pdb} {truth_pdb} -atom " C1\'" -m -'
    output = os.popen(command).read()
    if verbose == 1:
        print(output)
    tm_score = parse_usalign_for_tm_score(output)
    transform = parse_usalign_for_transform(output)
    return tm_score, transform

print('HELPER OK!!!')

# Protenix 模型推理流程
if MODEL_TYPE == 'protenix':
    from runner.batch_inference import get_default_runner
    from runner.inference import update_inference_configs, InferenceRunner
    from protenix.data.infer_data_pipeline import InferenceDataset

    # 设置随机种子
    np.random.seed(0)
    torch.random.manual_seed(0)
    torch.cuda.manual_seed_all(0)

    # 定义 DictDataset 类
    class DictDataset(InferenceDataset):
        def __init__(self, seq_list: list, dump_dir: str, id_list: list = None, use_msa: bool = False):
            self.dump_dir = dump_dir
            self.use_msa = use_msa
            if isinstance(id_list, type(None)):
                self.inputs = [{"sequences": [{"rnaSequence": {"sequence": seq, "count": 1}}], "name": "query"} for seq in seq_list]
            else:
                self.inputs = [{"sequences": [{"rnaSequence": {"sequence": seq, "count": 1}}], "name": i} for i, seq in zip(id_list, seq_list)]

# 配置 Protenix 模型
if MODEL_TYPE == 'protenix':
    from configs.configs_base import configs as configs_base
    from configs.configs_data import data_configs
    from configs.configs_inference import inference_configs
    from protenix.config.config import parse_configs

    configs_base["use_deepspeed_evo_attention"] = False  # 本地环境可能不需要 deepspeed
    configs_base["model"]["N_cycle"] = 10
    configs_base["sample_diffusion"]["N_sample"] = (1 if VALIDATION else 5)
    configs_base["sample_diffusion"]["N_step"] = 400
    inference_configs['load_checkpoint_path'] = 'C:/StanfordRNA3D/checkpoints/model_v0.2.0.pt'
    configs = {**configs_base, **{"data": data_configs}, **inference_configs}
    configs = parse_configs(configs, fill_required_with_null=True)
    runner = InferenceRunner(configs)

# 验证模式
if VALIDATION:
    LABEL_DF = pd.read_csv(os.path.join(DATA_LOCAL_DIR, 'train_labels.csv'))
    LABEL_DF['target_id'] = LABEL_DF['ID'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    train_df = pd.read_csv(os.path.join(DATA_LOCAL_DIR, 'train_sequences.csv'))

if MODEL_TYPE == 'protenix' and VALIDATION:
    import warnings
    warnings.filterwarnings("ignore")

    train_df['protenix_tm_score'] = None
    dataset = DictDataset(train_df.sequence, dump_dir='output', id_list=train_df.target_id, use_msa=False)
    num_data = len(dataset)

    time0 = time.time()
    for i, seq in tqdm(enumerate(train_df.sequence), total=num_data):
        if train_df.loc[i, 'protenix_tm_score'] is not None:
            continue
        if len(seq) > 300:
            continue
        target_id = train_df.loc[i, 'target_id']
        truth_df = get_truth_df(target_id)
        if sum(~np.isnan(truth_df.x_1)) < 3:
            continue

        data, atom_array, data_error_message = dataset[i]
        if data_error_message != '':
            continue
        new_configs = update_inference_configs(configs, data["N_token"].item())
        runner.update_model_configs(new_configs)
        prediction = runner.predict(data)
        prediction = prediction['coordinate'][:, data['input_feature_dict']['atom_to_tokatom_idx'] == 12]
        result = parse_output_to_df(prediction[:1], seq, target_id)[0]
        try:
            tm_score, transform = call_usalign(result, truth_df, verbose=0)
            train_df.loc[i, 'protenix_tm_score'] = tm_score
        except:
            pass
        if (time.time() - time0) > (12 * 3600 - 360):
            break
    train_df.to_csv('tm_scores.csv', index=False)
    print(train_df.protenix_tm_score.mean())
    train_df.protenix_tm_score.hist()
    plt.show()

# 提交模式
if MODEL_TYPE == 'protenix' and not VALIDATION:
    test_df = pd.read_csv(os.path.join(DATA_LOCAL_DIR, 'test_sequences.csv'))
    import warnings
    warnings.filterwarnings("ignore")
    dataset = DictDataset(test_df.sequence, dump_dir='output', id_list=test_df.target_id, use_msa=False)
    num_data = len(dataset)

    for i, seq in tqdm(enumerate(test_df.sequence), total=num_data):
        try:
            data, atom_array, data_error_message = dataset[i]
            target_id = data["sample_name"]
            assert target_id == test_df.target_id[i]
            assert data_error_message == ''
            new_configs = update_inference_configs(configs, data["N_token"].item())
            runner.update_model_configs(new_configs)
            prediction = runner.predict(data)
            prediction = prediction['coordinate'][:, data['input_feature_dict']['atom_to_tokatom_idx'] == 12]
            result = parse_output_to_df(prediction, seq, target_id)[0]
        except:
            target_id = test_df.target_id[i]
            print('Failed to predict', target_id)
            result = pd.DataFrame(columns=['ID', 'resname', 'resid',
                                           'x_1', 'y_1', 'z_1',
                                           'x_2', 'y_2', 'z_2',
                                           'x_3', 'y_3', 'z_3',
                                           'x_4', 'y_4', 'z_4',
                                           'x_5', 'y_5', 'z_5'],
                                 data=[[target_id, x, j+1] + [0.0]*15 for j, x in enumerate(seq)])
        result['ID'] = result.apply(lambda x: x.ID + '_' + str(x.resid), axis=1)
        result.to_csv('submission.csv', index=False, mode='a', header=(i==0))
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    submission_df = pd.read_csv('submission.csv')
    print(submission_df)