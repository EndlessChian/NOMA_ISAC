log = open('./log.txt', encoding='utf-8').read().replace('\t', '\n').split('\n')[1:-1]
log_adv = open('./log-adv.txt', encoding='utf-8').read().replace('\t', '\n').split('\n')[1:-1]
for i in range(len(log)):
    print(i,log[i])
import numpy as np
msg = dict([s.split('=') for s in log if '=' in s])
msg_adv = dict([s.split('=') for s in log_adv if '=' in s])
for k,v in msg.items():
    msg[k] = np.array(eval(v))
    msg_adv[k] = np.array(eval(msg_adv[k]))
print(msg)
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 6. 绘制中断概率结果
plt.figure(figsize=(10, 6))
plt.semilogy(msg['SNR_dB_outage'], msg['OP_NU_ISAC'], '-o', linewidth=1.5, label='Outage Probability NU (ISAC) Original')
plt.semilogy(msg['SNR_dB_outage'], msg_adv['OP_NU_ISAC'], '-x', linewidth=1.5, label='Outage Probability NU (ISAC) Improved')
plt.semilogy(msg['SNR_dB_outage'], msg['OP_NU_FDSAC'], '--o', linewidth=1.5, label='Outage Probability NU (FDSAC) Original')
plt.semilogy(msg['SNR_dB_outage'], msg_adv['OP_NU_FDSAC'], '--x', linewidth=1.5, label='Outage Probability NU (FDSAC) Improved')
plt.xlim(msg['SNR_dB_outage'][0], msg['SNR_dB_outage'][-1])
plt.xlabel('SNR (dB)', fontsize=12)
plt.ylabel('Outage Probability', fontsize=12)
plt.legend(loc='lower left', fontsize=10)
plt.title('Outage Probability for NU and FU in ISAC and FDSAC', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.savefig('outage_probability-comp.png')
plt.show()

# 7. 绘制感知速率结果
plt.figure(figsize=(10, 6))
plt.plot(msg['SNR_dB_sensing'], msg['sensing_rate_ISAC_sim'], '-s', linewidth=1.5, label='Sensing Rate (ISAC) - Original')
plt.plot(msg['SNR_dB_sensing'], msg['sensing_rate_FDSAC_sim'], '--o', linewidth=1.5, label='Sensing Rate (FDSAC) - Original')
plt.plot(msg['SNR_dB_sensing'], msg_adv['sensing_rate_ISAC_sim'], '-s', linewidth=1.5, label='Sensing Rate (ISAC) - Improved')
plt.plot(msg['SNR_dB_sensing'], msg_adv['sensing_rate_FDSAC_sim'], '--o', linewidth=1.5, label='Sensing Rate (FDSAC) - Improved')
plt.xlim(msg['SNR_dB_sensing'][0], msg['SNR_dB_sensing'][-1])
plt.xlabel('SNR (dB)', fontsize=12)
plt.ylabel('Sensing Rate (bps/Hz)', fontsize=12)
plt.legend(loc='upper left', fontsize=10)
plt.title('Sensing Rate vs SNR for ISAC and FDSAC', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.savefig('sensing_rate_vs_snr-comp.png')
plt.show()

# 8. 绘制遍历通信速率结果
plt.figure(figsize=(10, 6))
plt.plot(msg['SNR_dB_outage'], msg['ECR_NU_ISAC'], '-^', linewidth=1.5, label='ECR NU (ISAC) Original')
plt.plot(msg['SNR_dB_outage'], msg['ECR_FU_ISAC'], '-v', linewidth=1.5, label='ECR FU (ISAC) Original')
plt.plot(msg['SNR_dB_outage'], msg_adv['ECR_NU_ISAC'], '-^', linewidth=1.5, label='ECR NU (ISAC) Improved')
plt.plot(msg['SNR_dB_outage'], msg_adv['ECR_FU_ISAC'], '-v', linewidth=1.5, label='ECR FU (ISAC) Improved')
plt.xlim(msg['SNR_dB_outage'][0], msg['SNR_dB_outage'][-1])
plt.xlabel('SNR (dB)', fontsize=12)
plt.ylabel('Ergodic Communication Rate (bps/Hz)', fontsize=12)
plt.legend(loc='upper left', fontsize=10)
plt.title('Ergodic Communication Rate vs SNR for NU and FU in ISAC and FDSAC', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.savefig('ergodic_communication_rate-comp.png')
plt.show()

# 9. 在p=5dB时绘制感知速率与通信速率和的关系
plt.figure(figsize=(10, 6))
plt.plot(msg['sum_ECR_ISAC_mu'], msg['sensing_rate_ISAC_mu'], '-^', linewidth=2, label='ISAC (Original)')
plt.plot(msg['sum_ECR_FDSAC_mu'], msg['sensing_rate_FDSAC_mu'], '--o', linewidth=2, label='FDSAC (Original)')
plt.plot(msg_adv['sum_ECR_ISAC_mu'], msg_adv['sensing_rate_ISAC_mu'], '-^', linewidth=2, label='ISAC (Improved)')
plt.plot(msg_adv['sum_ECR_FDSAC_mu'], msg_adv['sensing_rate_FDSAC_mu'], '--o', linewidth=2, label='FDSAC (Improved)')
plt.xlabel('Sum Communication Rate (bps/Hz)', fontsize=12)
plt.ylabel('Sensing Rate (bps/Hz)', fontsize=12)
plt.legend(loc='upper left', fontsize=10)
plt.title(f'Sensing Rate vs Sum Communication Rate', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.savefig('sensing_vs_communication_rate-comp.png')
plt.show()