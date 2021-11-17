import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    # 创建一个图像窗口
    fig = plt.figure()
    # 在图像窗口添加3d坐标轴
    # ax = Axes3D(fig)

    # 使用np.linspace定义 x:范围(-10,10);个数为100
    x1 = np.linspace(0, 3, 100)
    # 定义 y:范围(-3,3);个数为50
    x2 = np.linspace(0, 3, 100)
    # 创建x-y平面网络
    x1, x2 = np.meshgrid(x1, x2)
    # 定义函数 r=1/2*(x-y)^2
    # r = 1/2*np.square(x-y)
    r = np.power(x1 - 9/4, 2) + np.power(x2 - 2, 2)
    xt1 = [
        0.000000,
        0.675000,
        0.934338,
        1.085225,
        1.186751,
        1.223159,
        1.539390,
        1.535426,
        1.531854,
        1.528638,
        1.525743,
        1.523137,
        1.520793,
        1.518685,
        1.516789,
        1.515084,
        1.513552,
        1.512174,
        1.510936,
        1.509823
    ]
    xt2 = [
        0.000000,
        0.600000,
        0.921231,
        1.193542,
        1.410734,
        1.496913,
        2.269723,
        2.267517,
        2.265551,
        2.263800,
        2.262242,
        2.260856,
        2.259625,
        2.258531,
        2.257559,
        2.256697,
        2.255931,
        2.255252,
        2.254650,
        2.254115
    ]

    # 将函数显示为3d  rstride 和 cstride 代表 row(行)和column(列)的跨度 get_cmap为色图分类
    # ax.plot_surface(x1, x2, r)
    plt.contourf(x1, x2, r)
    plt.contour(x1, x2, r, 100)
    plt.plot(xt1, xt2)
    # 投影

    # 显示创建的图像
    plt.show()
