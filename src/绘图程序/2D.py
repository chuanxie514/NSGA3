import numpy as np
import matplotlib.pyplot as plt

# 读取.data文件中的数据
data = np.genfromtxt('NSGA3_ZDT6_2D_R3.data', delimiter='\t')
x = data[:, 0]
y = data[:, 1]

# 绘制散点图
plt.scatter(x, y, c='g', marker='o')

# 设置背景网格线
plt.grid(True, linestyle='--', alpha=0.5)

# 设置标题
plt.title('ZDT6')

# 设置坐标轴标签
plt.xlabel('F1')
plt.ylabel('F2')

# 显示图形
plt.show()
