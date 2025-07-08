import math, random
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import det
import time, datetime

# 设置随机种子以确保结果可复现
np.random.seed(42)

# 1. 定义系统参数
# 功率分配因子（主模拟中动态变化）
alpha_N_fixed = 0.2  # ISAC中近用户(NU)的基础功率分配
alpha_F_fixed = 0.8  # ISAC中远用户(FU)的基础功率分配

# 噪声方差
sigma_c = 0.01  # 通信噪声方差
sigma_s = 0.01  # 感知噪声方差

# 目标速率
target_rate_N = 0.8  # NU的目标速率 (bps/Hz)
target_rate_F = 0.8  # FU的目标速率 (bps/Hz)

# 蒙特卡洛模拟参数
num_trials = 10000  # 蒙特卡洛模拟试验次数

# SNR范围
SNR_dB_outage = np.arange(0, 61, 5)  # 中断概率和ECR分析的SNR范围 (dB)
SNR_dB_sensing = np.arange(-20, 21, 5)  # 感知速率分析的SNR范围 (dB)

# FDSAC参数
kappa_fixed = 0.5  # 分配给通信的固定带宽比例
mu_fixed = 0.5  # 分配给通信的固定功率比例

# 感知参数
L = 30  # 通信帧长度
M = 8  # 感知的距离单元数
lambda_mean = np.array([5, 3, 3.5, 2.5, 1.5, 2, 1, 0.5])  # R的平均特征值
num_eigen = len(lambda_mean)  # 特征值数量

# 将目标速率转换为SINR阈值
gamma_target_N = 2 ** target_rate_N - 1  # NU的SINR阈值
gamma_target_F = 2 ** target_rate_F - 1  # FU的SINR阈值

# 2. 初始化存储结果的数组
# 中断概率
OP_NU_ISAC = np.zeros(len(SNR_dB_outage))
OP_FU_ISAC = np.zeros(len(SNR_dB_outage))
OP_NU_FDSAC = np.zeros(len(SNR_dB_outage))
OP_FU_FDSAC = np.zeros(len(SNR_dB_outage))

# 感知速率
sensing_rate_ISAC_sim = np.zeros(len(SNR_dB_sensing))
sensing_rate_FDSAC_sim = np.zeros(len(SNR_dB_sensing))

# 遍历通信速率 (ECR)
ECR_NU_ISAC = np.zeros(len(SNR_dB_outage))
ECR_FU_ISAC = np.zeros(len(SNR_dB_outage))
ECR_NU_FDSAC = np.zeros(len(SNR_dB_outage))
ECR_FU_FDSAC = np.zeros(len(SNR_dB_outage))

rand_complex = lambda: random.random() + 1j * random.random()
abs2 = lambda n: n.real * n.real + n.imag * n.imag
eyeM = np.eye(M).astype(np.float32)
import os

logger = r'./log-adv.txt'
if os.path.isfile(logger):
    print('Found old logger file:')
    with open(logger, 'r', encoding='utf-8') as f:
        print(f.read())
    with open(logger, 'w', encoding='utf-8') as f:
        pass


# 改进的功率分配函数
def adaptive_power_allocation(hN, hF):
    # 计算信道增益比
    channel_ratio = abs2(hF) / (abs2(hN) + 1e-10)

    # 基础功率分配
    base_alpha_N = 0.2
    base_alpha_F = 0.8

    # 自适应调整因子 - 基于信道差异
    adjustment_factor = 0.3 * (1 - math.tanh(5 * (channel_ratio - 0.5)))

    # 应用调整
    alpha_N = np.clip(base_alpha_N + adjustment_factor, 0.1, 0.4)
    alpha_F = np.clip(base_alpha_F - adjustment_factor, 0.6, 0.9)

    # 确保功率分配有效
    total = alpha_N + alpha_F
    if total > 1:
        alpha_F = 1 - alpha_N

    return alpha_N, alpha_F

# 3. 带动态功率分配的中断概率和ECR模拟
sigma_c2 = sigma_c ** 2
AN = 0.9
AF = 0.2
_kAN = math.sqrt(AN / 2)
_kAF = math.sqrt(AF / 2)
start_time = time.time()

for idx, snr_db in enumerate(SNR_dB_outage):
    SNR_linear = 10 ** (snr_db / 10)  # 将SNR从dB转换为线性标度
    p = SNR_linear * sigma_c2  # 根据SNR值调整功率预算

    # 初始化中断概率计数器
    outage_count_NU_ISAC = 0
    outage_count_FU_ISAC = 0
    outage_count_NU_FDSAC = 0
    outage_count_FU_FDSAC = 0

    # 初始化遍历通信速率累加器
    sum_RN_ISAC = 0
    sum_RF_ISAC = 0
    sum_RN_FDSAC = 0
    sum_RF_FDSAC = 0

    for trial in range(num_trials):
        # 为NU和FU生成瑞利衰落信道
        hN = _kAN * rand_complex()
        hF = _kAF * rand_complex()
        HN = p * abs2(hN)
        HF = p * abs2(hF)

        # 基于实时信道反馈的动态功率分配
        alpha_N, alpha_F = adaptive_power_allocation(hN, hF)

        # 计算ISAC的SINR
        gamma_sic_NU_ISAC = (HN * alpha_F) / (sigma_c2 + HN * alpha_N)
        gamma_NU_ISAC = (HN * alpha_N) / sigma_c2
        gamma_FU_ISAC = (HF * alpha_F) / (sigma_c2 + HF * alpha_N)

        # 确定ISAC的中断
        if not (gamma_sic_NU_ISAC > gamma_target_F and gamma_NU_ISAC > gamma_target_N):
            outage_count_NU_ISAC += 1
        if gamma_FU_ISAC < gamma_target_F:
            outage_count_FU_ISAC += 1

        # 计算FDSAC的SINR
        gamma_NU_FDSAC = (mu_fixed * HN) / sigma_c2
        gamma_FU_FDSAC = (mu_fixed * HF) / sigma_c2

        # 确定FDSAC的中断
        if gamma_NU_FDSAC < gamma_target_N:
            outage_count_NU_FDSAC += 1
        if gamma_FU_FDSAC < gamma_target_F:
            outage_count_FU_FDSAC += 1

        # 累加ISAC的遍历通信速率
        sum_RN_ISAC += math.log2(1 + gamma_NU_ISAC)
        sum_RF_ISAC += math.log2(1 + gamma_FU_ISAC)

        # 累加FDSAC的遍历通信速率（按kappa缩放）
        sum_RN_FDSAC += kappa_fixed * math.log2(1 + gamma_NU_FDSAC)
        sum_RF_FDSAC += kappa_fixed * math.log2(1 + gamma_FU_FDSAC)

    # 计算中断概率
    OP_NU_ISAC[idx] = outage_count_NU_ISAC / num_trials
    OP_FU_ISAC[idx] = outage_count_FU_ISAC / num_trials
    OP_NU_FDSAC[idx] = outage_count_NU_FDSAC / num_trials
    OP_FU_FDSAC[idx] = outage_count_FU_FDSAC / num_trials

    # 计算遍历通信速率（试验平均值）
    ECR_NU_ISAC[idx] = sum_RN_ISAC / num_trials
    ECR_FU_ISAC[idx] = sum_RF_ISAC / num_trials
    ECR_NU_FDSAC[idx] = sum_RN_FDSAC / num_trials
    ECR_FU_FDSAC[idx] = sum_RF_FDSAC / num_trials

# 4. 感知速率模拟
sigma_s2 = sigma_s ** 2
perturbation_std = 0.1  # 扰动的标准差

for idx, snr_db in enumerate(SNR_dB_sensing):
    SNR_linear = 10 ** (snr_db / 10)  # 将SNR从dB转换为线性标度
    p = SNR_linear * sigma_s2  # 根据SNR值调整功率预算

    # 初始化感知速率累加器
    sum_Rs_ISAC = 0
    sum_Rs_FDSAC = 0

    for trial in range(num_trials):
        # 通过扰动特征值引入相关矩阵R的随机性
        lambda_perturbed = lambda_mean + perturbation_std * np.random.randn(num_eigen)
        lambda_perturbed[lambda_perturbed < 0] = 0  # 确保特征值非负

        # 构建随机相关矩阵R
        R_random = np.diag(lambda_perturbed)

        # ISAC感知速率计算
        sum_Rs_ISAC += (1 / L) * math.log2(det(eyeM + (p * L / sigma_s2) * R_random))

        # FDSAC感知速率计算
        sum_Rs_FDSAC += ((1 - kappa_fixed) / L) * math.log2(
            det(eyeM + ((1 - mu_fixed) * p * L / ((1 - kappa_fixed) * sigma_s2)) * R_random)
        )

    # 计算试验的平均感知速率
    sensing_rate_ISAC_sim[idx] = sum_Rs_ISAC / num_trials
    sensing_rate_FDSAC_sim[idx] = sum_Rs_FDSAC / num_trials

# 5. 在p=5dB时的感知速率与通信速率和的关系图
p_fixed_dB = 5
p_fixed = 10 ** (p_fixed_dB / 10) * sigma_c2  # 转换为线性标度并调整功率

# 定义mu（功率分配因子）的范围
mu_range = np.arange(0, 1.05, 0.05)  # 从0到1，步长0.05
num_mu = len(mu_range)

# 初始化存储FDSAC的感知速率和通信速率和的数组
sensing_rate_FDSAC_mu = np.zeros(num_mu)
sum_ECR_FDSAC_mu = np.zeros(num_mu)

# 初始化存储ISAC的感知速率和通信速率和的数组
sensing_rate_ISAC_mu = np.zeros(num_mu)
sum_ECR_ISAC_mu = np.zeros(num_mu)

# 对变化的mu进行模拟
perturbation_std = 0.1
for mu_idx, mu in enumerate(mu_range):
    # 初始化累加器
    sum_RN_FDSAC_mu = 0
    sum_RF_FDSAC_mu = 0
    sum_RN_ISAC_mu = 0
    sum_RF_ISAC_mu = 0
    sum_Rs_FDSAC_mu_val = 0
    sum_Rs_ISAC_mu_val = 0

    for trial in range(num_trials):
        # 为NU和FU生成瑞利衰落信道
        hN = _kAN * rand_complex()
        hF = _kAF * rand_complex()
        HN = p_fixed * abs2(hN)
        HF = p_fixed * abs2(hF)

        # 基于实时信道反馈的动态功率分配
        alpha_N, alpha_F = adaptive_power_allocation(hN, hF)

        # ISAC遍历通信速率计算
        gamma_sic_NU_ISAC = (HN * alpha_F) / (sigma_c2 + HN * alpha_N)
        gamma_NU_ISAC = (HN * alpha_N) / sigma_c2
        gamma_FU_ISAC = (HF * alpha_F) / (sigma_c2 + HF * alpha_N)

        sum_RN_ISAC_mu += math.log2(1 + gamma_NU_ISAC)
        sum_RF_ISAC_mu += math.log2(1 + gamma_FU_ISAC)

        # FDSAC遍历通信速率计算
        gamma_NU_FDSAC_mu = (mu * HN) / sigma_c2
        gamma_FU_FDSAC_mu = (mu * HF) / sigma_c2

        sum_RN_FDSAC_mu += kappa_fixed * math.log2(1 + gamma_NU_FDSAC_mu)
        sum_RF_FDSAC_mu += kappa_fixed * math.log2(1 + gamma_FU_FDSAC_mu)

        # 在p=5dB时计算当前mu的感知速率
        lambda_perturbed = lambda_mean + perturbation_std * np.random.randn(num_eigen)
        lambda_perturbed[lambda_perturbed < 0] = 0
        R_random = np.diag(lambda_perturbed)

        # ISAC感知速率计算
        sum_Rs_ISAC_mu_val += (1 / L) * math.log2(det(eyeM + (p_fixed * L / sigma_s2) * R_random))

        # FDSAC感知速率计算
        sum_Rs_FDSAC_mu_val += ((1 - kappa_fixed) / L) * math.log2(
            det(eyeM + ((1 - mu) * p_fixed * L / ((1 - kappa_fixed) * sigma_s2)) * R_random)
        )

    # 计算平均遍历通信速率
    ECR_NU_FDSAC_mu = sum_RN_FDSAC_mu / num_trials
    ECR_FU_FDSAC_mu = sum_RF_FDSAC_mu / num_trials
    sum_ECR_FDSAC_mu[mu_idx] = ECR_NU_FDSAC_mu + ECR_FU_FDSAC_mu

    ECR_NU_ISAC_mu_val = sum_RN_ISAC_mu / num_trials
    ECR_FU_ISAC_mu_val = sum_RF_ISAC_mu / num_trials
    sum_ECR_ISAC_mu[mu_idx] = ECR_NU_ISAC_mu_val + ECR_FU_ISAC_mu_val

    # 计算平均感知速率
    sensing_rate_ISAC_mu[mu_idx] = sum_Rs_ISAC_mu_val / num_trials
    sensing_rate_FDSAC_mu[mu_idx] = sum_Rs_FDSAC_mu_val / num_trials

used_time = time.time() - start_time

with open(logger, 'w', encoding='utf-8') as f:
    print('Start at:', datetime.datetime.now(), file=f)
    # 6. 绘制中断概率结果
    plt.figure(figsize=(10, 6))
    print(f'SNR_dB_outage={SNR_dB_outage.tolist()}', file=f)
    plt.semilogy(SNR_dB_outage, OP_NU_ISAC, '-o', linewidth=1.5, label='Outage Probability NU (ISAC)')
    plt.semilogy(SNR_dB_outage, OP_FU_ISAC, '-x', linewidth=1.5, label='Outage Probability FU (ISAC)')
    plt.semilogy(SNR_dB_outage, OP_NU_FDSAC, '--o', linewidth=1.5, label='Outage Probability NU (FDSAC)')
    plt.semilogy(SNR_dB_outage, OP_FU_FDSAC, '--x', linewidth=1.5, label='Outage Probability FU (FDSAC)')
    plt.xlim(SNR_dB_outage[0], SNR_dB_outage[-1])
    print(f'\tOP_NU_ISAC={OP_NU_ISAC.tolist()}', file=f)
    print(f'\tOP_FU_ISAC={OP_FU_ISAC.tolist()}', file=f)
    print(f'\tOP_NU_FDSAC={OP_NU_FDSAC.tolist()}', file=f)
    print(f'\tOP_FU_FDSAC={OP_FU_FDSAC.tolist()}', file=f)
    plt.xlabel('SNR (dB)', fontsize=12)
    plt.ylabel('Outage Probability', fontsize=12)
    plt.legend(loc='lower left', fontsize=10)
    plt.title('Outage Probability for NU and FU in ISAC and FDSAC', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('outage_probability-adv.png')
    plt.show()

    # 7. 绘制感知速率结果
    plt.figure(figsize=(10, 6))
    print(f'SNR_dB_sensing={SNR_dB_sensing.tolist()}', file=f)
    plt.plot(SNR_dB_sensing, sensing_rate_ISAC_sim, '-s', linewidth=1.5, label='Sensing Rate (ISAC) - Simulated')
    print(f'\tsensing_rate_ISAC_sim={sensing_rate_ISAC_sim.tolist()}', file=f)
    plt.plot(SNR_dB_sensing, sensing_rate_FDSAC_sim, '--o', linewidth=1.5, label='Sensing Rate (FDSAC) - Simulated')
    print(f'\tsensing_rate_FDSAC_sim={sensing_rate_FDSAC_sim.tolist()}', file=f)
    plt.xlim(SNR_dB_sensing[0], SNR_dB_sensing[-1])
    plt.xlabel('SNR (dB)', fontsize=12)
    plt.ylabel('Sensing Rate (bps/Hz)', fontsize=12)
    plt.legend(loc='upper left', fontsize=10)
    plt.title('Sensing Rate vs SNR for ISAC and FDSAC', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('sensing_rate_vs_snr-adv.png')
    plt.show()

    # 8. 绘制遍历通信速率结果
    plt.figure(figsize=(10, 6))
    print(f'SNR_dB_outage={SNR_dB_outage.tolist()}', file=f)
    plt.plot(SNR_dB_outage, ECR_NU_ISAC, '-^', linewidth=1.5, label='ECR NU (ISAC)')
    plt.plot(SNR_dB_outage, ECR_FU_ISAC, '-v', linewidth=1.5, label='ECR FU (ISAC)')
    plt.plot(SNR_dB_outage, ECR_NU_FDSAC, '--^', linewidth=1.5, label='ECR NU (FDSAC)')
    plt.plot(SNR_dB_outage, ECR_FU_FDSAC, '--v', linewidth=1.5, label='ECR FU (FDSAC)')
    plt.xlim(SNR_dB_outage[0], SNR_dB_outage[-1])
    print(f'\tECR_NU_ISAC={ECR_NU_ISAC.tolist()}', file=f)
    print(f'\tECR_FU_ISAC={ECR_FU_ISAC.tolist()}', file=f)
    print(f'\tECR_NU_FDSAC={ECR_NU_FDSAC.tolist()}', file=f)
    print(f'\tECR_FU_FDSAC={ECR_FU_FDSAC.tolist()}', file=f)
    plt.xlabel('SNR (dB)', fontsize=12)
    plt.ylabel('Ergodic Communication Rate (bps/Hz)', fontsize=12)
    plt.legend(loc='upper left', fontsize=10)
    plt.title('Ergodic Communication Rate vs SNR for NU and FU in ISAC and FDSAC', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('ergodic_communication_rate-adv.png')
    plt.show()

    # 9. 在p=5dB时绘制感知速率与通信速率和的关系
    plt.figure(figsize=(10, 6))
    plt.plot(sum_ECR_ISAC_mu, sensing_rate_ISAC_mu, '-^', linewidth=2, label='ISAC')
    plt.plot(sum_ECR_FDSAC_mu, sensing_rate_FDSAC_mu, '--o', linewidth=2, label='FDSAC')
    print(f'sum_ECR_ISAC_mu={sum_ECR_ISAC_mu.tolist()}\tsensing_rate_ISAC_mu={sensing_rate_ISAC_mu.tolist()}', file=f)
    print(f'sum_ECR_FDSAC_mu={sum_ECR_FDSAC_mu.tolist()}\tsensing_rate_FDSAC_mu={sensing_rate_FDSAC_mu.tolist()}',
          file=f)
    plt.xlabel('Sum Communication Rate (bps/Hz)', fontsize=12)
    plt.ylabel('Sensing Rate (bps/Hz)', fontsize=12)
    plt.legend(loc='upper left', fontsize=10)
    plt.title(f'Sensing Rate vs Sum Communication Rate at p = {p_fixed_dB} dB', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('sensing_vs_communication_rate-adv.png')
    plt.show()
    print('Stop at:', datetime.datetime.now(), file=f)

print(used_time)