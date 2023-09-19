import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 读取.data文件中的数据
data = np.genfromtxt('NSGA3_DTLZ1_3D_R10.data', delimiter='\t')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# 创建三维坐标系
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 设置标题
plt.title('DTLZ1_3D_R10')

# 绘制散点图
ax.scatter(x, y, z, c='b', marker='o')

# 设置坐标轴标签
ax.set_xlabel('F0')
ax.set_ylabel('F1')
ax.set_zlabel('F2')

# 显示图形
plt.show()
