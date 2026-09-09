#MODEL_TYPE用于指定模型类型
MODEL_TYPE='protenix'

#VALIDATION用于控制是否进行验证步骤。
VALIDATION=False

#如果 MODEL_TYPE 是 'protenix' 且 VALIDATION 为 True，则安装以下Python库。
#get_ipython().system() 是 Jupyter Notebook 中的一个方法，用于在 Notebook 中运行系统命令。
if MODEL_TYPE=='protenix' and VALIDATION:
    get_ipython().system('pip install --no-deps protenix')
    get_ipython().system('pip install biopython')
    get_ipython().system('pip install ml-collections')
    get_ipython().system('pip install biotite==1.0.1')
    get_ipython().system('pip install rdkit')
get_ipython().system('export PROTENIX_DATA_ROOT_DIR=/kaggle/input/protenix-checkpoints')


#这部分代码执行了以下文件系统操作：
get_ipython().system(' mkdir /af3-dev') #创建一个名为 /af3-dev 的目录
get_ipython().system(' ln -s /kaggle/input/protenix-checkpoints /af3-dev/release_data') #创建一个符号链接
get_ipython().system(' ls /af3-dev/release_data/')  #列出 /af3-dev/release_data/ 目录中的内容


#导入相关的python包
import Bio
from copy import deepcopy
import pandas as pd
from Bio.PDB import Atom, Model, Chain, Residue, Structure, PDBParser
from Bio import SeqIO
import os, sys
import re
import numpy as np
import torch
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
import time


#获取当前运行的Python解释器的路径
PYTHON = sys.executable
print('PYTHON',PYTHON)

#指向 RhoFold-main 目录，这是一个可能用于RNA结构预测的代码库。
RHONET_DIR='/kaggle/input/data-for-demo-for-rhofold-plus-with-kaggle-msa/RhoFold-main'

#指向 USalign 工具的路径，USalign 是一个用于结构比对的工具。
USALIGN = '/kaggle/working//USalign'

#将 /kaggle/input/usalign/USalign 文件复制到 /kaggle/working/ 目录下。
os.system('cp /kaggle/input/usalign/USalign /kaggle/working/')
os.system('sudo chmod u+x /kaggle/working//USalign')

#添加路径到系统路径
sys.path.append(RHONET_DIR)

#存储RNA三维结构数据的目录
DATA_KAGGLE_DIR = '/kaggle/input/stanford-rna-3d-folding'


#定义了一个名为 dotdict 的类，继承自 dict
class dotdict(dict):
	__setattr__ = dict.__setitem__
	__delattr__ = dict.__delitem__

	def __getattr__(self, name):
		try:
			return self[name]
		except KeyError:
			raise AttributeError(name)

#设置3D图形等比例缩放.这个函数通常用于3D图形绘制时，确保三个轴的缩放比例一致，避免因不同轴的范围差异导致图形变形。
def set_aspect_equal(ax):       #ax参数是一个3D图形的轴对象
	x_limits = ax.get_xlim()
	y_limits = ax.get_ylim()
	z_limits = ax.get_zlim()

	#计算每个轴的平均值
	x_middle = np.mean(x_limits)
	y_middle = np.mean(y_limits)
	z_middle = np.mean(z_limits)

	#计算每个轴的范围
	max_range = max(x_limits[1] - x_limits[0],
					y_limits[1] - y_limits[0],
					z_limits[1] - z_limits[0]) / 2.0

	# 设定新的限制，以确保每个轴在同等范围
	ax.set_xlim(x_middle - max_range, x_middle + max_range)
	ax.set_ylim(y_middle - max_range, y_middle + max_range)
	ax.set_zlim(z_middle - max_range, z_middle + max_range)




#从一个大的数据集中提取特定目标的数据，通常用于处理标签数据。
def get_truth_df(target_id):
    #从全局变量 LABEL_DF中筛选出 target_id 列等于指定 target_id 的行。
    truth_df = LABEL_DF[LABEL_DF['target_id'] == target_id]
    truth_df = truth_df.reset_index(drop=True)
    return truth_df

#函数将模型的输出（3D坐标）解析为一个结构化的 DataFrame，方便后续的分析和可视化。
#参数output：模型的输出，通常是一个3D张量，形状为 (num_models, seq_length, 3)
#参数seq：序列数据，通常是一个字符串或列表，表示每个残基的名称（例如氨基酸名称）。
def parse_output_to_df(output, seq, target_id):
    df = []
    chain_data = []
    for i, res in enumerate(seq):
        d=dict(ID = target_id,
                    resname=res,
                    resid=i+1)
        #遍历模型的输出，提取坐标并更新字典，将残基数据添加到列表chain_data
        for n in range(len(output)):
            d={**d, f'x_{n+1}': round(output[n,i,0].item(),3),
                     f'y_{n+1}': round(output[n,i,1].item(),3),
                     f'z_{n+1}': round(output[n,i,2].item(),3)}
        chain_data.append(d)

    if len(chain_data)!=0:
        chain_df = pd.DataFrame(chain_data)
        df.append(chain_df)
        ##print(chain_df)
    return df

#将一个PDB（蛋白质数据银行）文件解析为一个或多个 pandas.DataFrame
#PDB文件是一种常见的生物分子结构数据格式，通常用于存储蛋白质、核酸等生物大分子的三维结构信息。
def parse_pdb_to_df(pdb_file, target_id):
    #PDBParser 是 Bio.PDB 模块中的一个类，用于解析PDB文件。
    parser = PDBParser()
    structure = parser.get_structure('', pdb_file)

    df = []
    #structure 对象包含一个或多个模型（Model）,遍历当前链中的每个残基（Residue）。
    for model in structure:
        for chain in model:
            print(chain)
            chain_data = []
            for residue in chain:
                #检查残基名称是否为 'A'（腺苷酸）、'U'（尿苷酸）、'G'（鸟苷酸）或 'C'（胞嘧啶核苷酸），这些是RNA的四种碱基。
                if residue.get_resname() in ['A', 'U', 'G', 'C']:
                    # 检查当前残基是否包含 'C1\' 原子。C1\' 是RNA残基中的一个重要原子，通常用于表示残基的位置。
                    if 'C1\'' in residue:
                        atom = residue['C1\'']
                        xyz = atom.get_coord()
                        resname = residue.get_resname()
                        resid = residue.get_id()[1]

                        #将残基的三维坐标，放入字典列表chain_data
                        #ID	resname	resid	x_1	y_1	z_1
                        chain_data.append(dict(
                            ID = target_id+'_'+str(resid),
                            resname=resname,
                            resid=resid,
                            x_1=xyz[0],
                            y_1=xyz[1],
                            z_1=xyz[2],
                        ))
                        ##print(f"Residue {resname} {resid}, Atom: {atom.get_name()}, xyz: {xyz}")

            if len(chain_data)!=0:
                chain_df = pd.DataFrame(chain_data)
                df.append(chain_df)
                ##print(chain_df)
    return df

#函数 write_target_line 用于生成一条符合PDB（蛋白质数据银行）格式的字符串。
def write_target_line(
    atom_name, atom_serial, residue_name, chain_id, residue_num, x_coord, y_coord, z_coord, occupancy=1.0, b_factor=0.0, atom_type='P'
):
    # atom_name（str）：原子名称，例如
    # "N"（氮原子）、"CA"（α碳原子）。
    # atom_serial（int）：原子的序列号。
    # residue_name（str）：残基名称，例如
    # "ALA"（丙氨酸）、"GLY"（甘氨酸）。
    # chain_id（str）：链标识符，例如
    # "A"、"B"。
    # residue_num（int）：残基编号。
    # x_coord（float）：原子的X坐标。
    # y_coord（float）：原子的Y坐标。
    # z_coord（float）：原子的Z坐标。
    # occupancy（float，可选）：占据率，默认值为
    # 1.0。
    # b_factor（float，可选）：B因子（温度因子），默认值为
    # 0.0。
    # atom_type（str，可选）：原子类型，默认值为
    # 'P'。
    return f'ATOM  {atom_serial:>5d}  {atom_name:<5s} {residue_name:<3s} {residue_num:>3d}    {x_coord:>8.3f}{y_coord:>8.3f}{z_coord:>8.3f}{occupancy:>6.2f}{b_factor:>6.2f}           {atom_type}\n'

#write_xyz_to_pdb函数的主要功能是将包含坐标信息的数据框（DataFrame）转换为PDB格式的文件
def write_xyz_to_pdb(df, pdb_file, xyz_id = 1):
    # df: 一个pandas DataFrame对象，预期包含了分子结构的坐标信息。
    # pdb_file: 字符串类型，目标PDB文件的路径，用于指定输出文件的位置和名称。
    # xyz_id: 整型，默认值为1，用于构建DataFrame中坐标列名的一部分（如'x_1', 'y_1', 'z_1'等），以此来区分不同的坐标集。
    resolved_cnt = 0
    with open(pdb_file, 'w') as target_file:    #打开目标PDB文件
        for _, row in df.iterrows():
            x_coord = row[f'x_{xyz_id}']
            y_coord = row[f'y_{xyz_id}']
            z_coord = row[f'z_{xyz_id}']
            #检查坐标有效性,过滤掉一些无效或缺失的数据
            if x_coord > -1e17 and y_coord > -1e17 and z_coord > -1e17:
                resolved_cnt += 1
                #对于每个有效坐标，调用write_target_line函数生成一条符合PDB文件格式的记录行，并将其写入到目标文件中。
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

#使用正则表达式提取 TM-score。USalign 是一种用于结构比对的工具，可以计算两个结构之间的相似性（如 TM-score）并提供变换矩阵。
def parse_usalign_for_tm_score(output):
    # Extract TM-score based on length of reference structure (second)
    tm_score_match = re.findall(r'TM-score=\s+([\d.]+)', output)[1]
    if not tm_score_match:
        raise ValueError('No TM score found')
    return float(tm_score_match)

#遍历 output 的每一行，寻找包含旋转矩阵的部分。
def parse_usalign_for_transform(output):
    # Locate the rotation matrix section
    matrix_lines = []
    found_matrix = False

    #遍历 output 的每一行，如果找到包含 "The rotation matrix to rotate Structure_1 to Structure_2" 的行，标记为找到矩阵部分。
    for line in output.splitlines():
        if "The rotation matrix to rotate Structure_1 to Structure_2" in line:
            found_matrix = True
        elif found_matrix and re.match(r'^\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+$', line):
            matrix_lines.append(line)
        elif found_matrix and not line.strip():
            break  # 如果在矩阵后遇到空行，则停止解析

    #解析旋转矩阵的值
    rotation_matrix = []
    for line in matrix_lines:
        parts = line.split()
        row_values = list(map(float, parts[1:]))  # 跳过第一列（索引）
        rotation_matrix.append(row_values)

    return np.array(rotation_matrix)

#将预测结构和真实结构的 DataFrame 转换为 PDB 文件
def call_usalign(predict_df, truth_df, verbose=1):
    truth_pdb = '~truth.pdb'
    predict_pdb = '~predict.pdb'
    write_xyz_to_pdb(predict_df, predict_pdb, xyz_id=1)
    write_xyz_to_pdb(truth_df, truth_pdb, xyz_id=1)
    #构造 USalign 命令
    command = f'{USALIGN} {predict_pdb} {truth_pdb} -atom " C1\'" -m -'
    #执行 USalign 命令并读取输出。
    output = os.popen(command).read()
    if verbose==1:
        print(output)
    #使用 parse_usalign_for_tm_score 和 parse_usalign_for_transform 解析输出，提取 TM-score 和变换矩阵。
    tm_score = parse_usalign_for_tm_score(output)
    transform = parse_usalign_for_transform(output)
    return tm_score, transform

print('HELPER OK!!!')


#定义了一个基于 Protenix 模型的推理流程，包括数据准备、模型加载和推理配置。
if MODEL_TYPE=='protenix':
    #导入了 Protenix 相关的模块和类
    from runner.batch_inference import get_default_runner
    from runner.inference import update_inference_configs, InferenceRunner
    from protenix.data.infer_data_pipeline import InferenceDataset

    #设置随机种子
    np.random.seed(0)
    torch.random.manual_seed(0)
    torch.cuda.manual_seed_all(0)

    #定义了一个名为 DictDataset 的类，继承自 InferenceDataset，用于处理推理数据。
    #这个类可以用于将用户提供的序列数据转换为模型可以接受的格式，为后续的推理任务做准备
    class DictDataset(InferenceDataset):
        # seq_list（list）：包含序列的列表，每个序列是一个字符串。
        # dump_dir（str）：输出目录路径，用于存储推理结果。
        # id_list（list，可选）：包含序列ID的列表，与 seq_list 一一对应。如果未提供，则默认为 None。
        # use_msa（bool，可选）：是否使用多序列比对（MSA）数据，默认为 False。
        def __init__(
            self,
            seq_list: list,
            dump_dir: str,
            id_list: list = None,
            use_msa: bool = False,
        ) -> None:

            self.dump_dir = dump_dir
            self.use_msa = use_msa

            #根据是否提供了 id_list，构建输入数据 self.inputs
            if isinstance(id_list,type(None)):
                self.inputs = [{"sequences": 
                                [{"rnaSequence": 
                                  {"sequence": seq, 
                                   "count": 1}}],
                                "name": "query"} for seq in seq_list]
            else:
                self.inputs = [{"sequences": 
                                [{"rnaSequence": 
                                  {"sequence": seq, 
                                   "count": 1}}],
                                "name": i} for i, seq in zip(id_list,seq_list)]


#配置和初始化 Protenix 模型的推理环境。它加载了多个配置文件，设置了模型参数，并初始化了推理运行器。
if MODEL_TYPE=='protenix':
    #导入了多个配置模块
    from configs.configs_base import configs as configs_base
    from configs.configs_data import data_configs
    from configs.configs_inference import inference_configs
    from protenix.config.config import parse_configs
    #动态设置环境变量
    configs_base["use_deepspeed_evo_attention"] = (
    os.environ.get("USE_DEEPSPEED_EVO_ATTENTION", False) == "true")
    #设置模型参数
    configs_base["model"]["N_cycle"] = 10 #10
    configs_base["sample_diffusion"]["N_sample"] = (1 if VALIDATION else 5)
    configs_base["sample_diffusion"]["N_step"] = 400
    #设置推理配置，设置推理时加载的模型检查点路径。
    inference_configs['load_checkpoint_path']='/kaggle/input/protenix-checkpoints/model_v0.2.0.pt'
    configs = {**configs_base, **{"data": data_configs}, **inference_configs}

    configs = parse_configs(
            configs=configs,
            fill_required_with_null=True,
        )
    #使用解析后的配置 configs 初始化 InferenceRunner，准备进行推理任务。
    runner=InferenceRunner(configs)


#用于在验证模式下加载和处理数据,加载训练集的标签数据和序列数据,对标签数据进行处理，提取目标ID
if VALIDATION:
    LABEL_DF = pd.read_csv('/kaggle/input/stanford-rna-3d-folding/train_labels.csv')
    LABEL_DF['target_id'] = LABEL_DF['ID'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    train_df=pd.read_csv('/kaggle/input/stanford-rna-3d-folding/train_sequences.csv')

#用于在验证模式（VALIDATION=False）下，对测试集中的序列进行推理，并生成提交文件。
#使用 Protenix 模型对训练集中的序列进行推理，并计算每个序列的 TM-score（结构相似性分数）。
if MODEL_TYPE=='protenix' and VALIDATION:
    import warnings
    warnings.filterwarnings("ignore")  

    #初始化 TM-score 列
    train_df['protenix_tm_score']=None
    #使用 DictDataset 类创建一个数据集对象
    dataset = DictDataset(train_df.sequence, dump_dir='output', id_list=train_df.target_id, use_msa=False)
    num_data = len(dataset)

    #使用 tqdm 包装 enumerate，显示进度条，total=num_data 指定总进度。
    for i, seq in tqdm(enumerate(train_df.sequence),total=num_data):
        #如果当前序列的 protenix_tm_score 已经被计算过（不为 None），则跳过该序列。
        if train_df.loc[i,'protenix_tm_score']!=None:
            continue
        #跳过过长的序列
        if len(seq)>300:
            continue
        #获取真实数据
        target_id = train_df.loc[i,'target_id']
        truth_df = get_truth_df(target_id)
        #检查真实数据的有效性
        if sum(~np.isnan(truth_df.x_1))<3:
            continue

        #获取数据和原子数组
        data, atom_array, data_error_message=dataset[i]
        #跳过有错误的数据
        if data_error_message!='':
            continue
        #更新推理配置
        new_configs = update_inference_configs(configs, data["N_token"].item())
        runner.update_model_configs(new_configs)
        #进行推理
        prediction = runner.predict(data)
        #提取预测坐标
        prediction=prediction['coordinate'][:,data['input_feature_dict']['atom_to_tokatom_idx']==12]
        #解析预测结果
        result = parse_output_to_df(prediction[:1], seq, target_id)[0]
        #计算 TM-score
        try:
            tm_score, transform = call_usalign(result, truth_df, verbose=0)
            train_df.loc[i,'protenix_tm_score']=tm_score
        except:
            pass
        if (time.time()-time0)>(12*3600-360):
            break
    #保存 TM-score 到.csv文件
    train_df.to_csv('tm_scores.csv', index=False)
    print(train_df.protenix_tm_score.mean())
    display(train_df.protenix_tm_score.hist())


# 用于在非验证模式（VALIDATION=False）下，对测试集中的序列进行推理，并生成提交文件。
if MODEL_TYPE=='protenix' and not VALIDATION:
    #加载测试数据
    test_df=pd.read_csv('/kaggle/input/stanford-rna-3d-folding/test_sequences.csv')
    #忽略警告
    import warnings
    warnings.filterwarnings("ignore")  
    #创建数据集
    dataset = DictDataset(test_df.sequence, dump_dir='output', id_list=test_df.target_id, use_msa=False)
    num_data = len(dataset)
    #遍历序列，获取数据和原子数组
    for i, seq in tqdm(enumerate(test_df.sequence),total=num_data):
        try:
            data, atom_array, data_error_message=dataset[i]
            target_id = data["sample_name"]
            assert target_id==test_df.target_id[i]
            assert data_error_message==''

            #更新推理配置
            new_configs = update_inference_configs(configs, data["N_token"].item())
            runner.update_model_configs(new_configs)

            #进行推理，提取预测坐标
            prediction = runner.predict(data)
            prediction=prediction['coordinate'][:,data['input_feature_dict']['atom_to_tokatom_idx']==12]
            #解析预测结果
            result = parse_output_to_df(prediction, seq, target_id)[0]
        except:
            #如果推理过程中出现异常，打印失败信息，并生成一个包含零坐标的 DataFrame。
            target_id==test_df.target_id[i]
            print('Failed to predict', target_id)
            result=pd.DataFrame(columns=['ID', 'resname', 'resid', 
                                         'x_1', 'y_1', 'z_1', 
                                         'x_2', 'y_2', 'z_2',
                                         'x_3', 'y_3', 'z_3', 
                                         'x_4', 'y_4', 'z_4', 
                                         'x_5', 'y_5', 'z_5'], 
                                         data=[[target_id, x, j+1] + [0.0]*15 for j, x in enumerate(seq)])
        #将预测结果解析为 DataFrame，并保存到提交文件 submission.csv 中
        result['ID']=result.apply(lambda x: x.ID + '_' + str(x.resid), axis=1)
        result.to_csv('submission.csv', index=False, mode='a', header=(i==0))
        torch.cuda.empty_cache()

    display(pd.read_csv('submission.csv'))

