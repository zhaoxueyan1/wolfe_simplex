import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    # 创建一个图像窗口
    fig = plt.figure()
    # 在图像窗口添加3d坐标轴
    # ax = Axes3D(fig)

    # 使用np.linspace定义 x:范围(-10,10);个数为100
    x1 = np.linspace(0, 4, 100)
    # 定义 y:范围(-3,3);个数为50
    x2 = np.linspace(0, 4, 100)
    # 创建x-y平面网络
    x1, x2 = np.meshgrid(x1, x2)
    # 定义函数 r=1/2*(x-y)^2
    # r = 1/2*np.square(x-y)
    r = np.power(x1 - 2, 2) + 4 * np.power(x2 - 3, 2)

    xt1 = [0,
           0.164399,
           0.386601,
           0.748617,
           1.695329,
           2.02837,
           2.018618,
           2.001552,
           2.001018,
           2.000085,
           ]
    xt2 = [0,
           0.986394,
           1.961394,
           2.893567,
           3.215648,
           2.272736,
           3.272688,
           2.272834,
           3.272833,
           2.272834,
           ]

    # 将函数显示为3d  rstride 和 cstride 代表 row(行)和column(列)的跨度 get_cmap为色图分类
    # ax.plot_surface(x1, x2, r)
    plt.contourf(x1, x2, r)
    plt.contour(x1, x2, r, 100)
    plt.plot(xt1, xt2)
    # 投影

    # 显示创建的图像
    plt.show()
