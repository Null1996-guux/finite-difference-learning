# Dirichlet 初边界值问题

考虑一维非齐次热传导方程 $ Dirichelet $  初边界值问题:

$$
\begin{align}
    \dfrac{\partial u}{\partial t} - a \dfrac{\partial^2 u}{\partial x^2} = f(x, t), \quad 0 < x < 1, \ 0 < t \leq T \\
    u(x, 0) = \varphi(x), \quad 0 < x \leq 1, \\
    u(0 ,t) = \alpha(t), \quad u(1, t) = \beta(t), \quad 0 < t \leq T
\end{align}
$$

考虑改问题的有限差分法, 其中 $a$ 为正常数,$ f(x, t), \varphi(x), \alpha(t), \beta(t) $为已知函数, $\varphi(0) = \alpha(0), \quad \varphi(1) = \beta(0) $称上述第二个式子为初值条件, 最后两个式子为边值条件.

# 记号以及网格剖分

考虑利用有限差分法对上述问题进行求解, 将求解区域
$$
\begin{equation} 
    \Omega = \left\lbrace (x, t) | 0 \leq x \leq 1, \ 0 \leq t \leq T \right\rbrace
\end{equation}
$$
进行网格剖分, 将区间 $[0, 1]$ 作 $m$ 等分, 将区间 $[0, T]$ 作 $n$ 等分, 并且记 $ h = \dfrac{1}{m}, \ \tau = \dfrac{T}{n}, x_i = ih, \ 0 \leq i \leq m; \ t_k = k\tau, \ 0 \leq k \leq n $分别将 $ h $ 和 $ \tau $ 称为时间步长和空间步长, 用两簇平行直线
$$
\begin{equation}
    x = x_i, \quad 0 \leq i \leq m \\
    t = t_k, \quad 0 \leq k \leq n
\end{equation}
$$
将 $ \Omega $剖分为矩形网格, 记:
$$
\begin{align}
    \Omega_h = \left\lbrace x_i | 0 \leq i \leq m\right\rbrace \\ 
    \Omega_{\tau} = \left\lbrace t_k | 0 \leq k \leq n\right\rbrace \\
    \Omega_{h \tau} = \Omega_h \times \Omega_{\tau}
\end{align}
$$
称$ \left\lbrace (x_i, t_k) \right\rbrace $为 $\color{red}{网格结点}$ ; 称在 $ t = 0, x = 0 $ 以及 $ x = 1 $ 上的结点为$ \color{red}{边界结点}$, 称其他结点为$\color{red}{内部结点}$, 称在直线 $ t = t_k $上的所有结点$ \left\lbrace (x_i, t_k ) | 0 \leq i \leq m \right\rbrace $ 为$\color{red}{第 k 层结点}$

设$\left\lbrace v_i^k | 0 \leq i \leq m, 0 \leq k \leq n\right\rbrace $为网格$ \Omega_{h \tau} $上的一个函数给出以下记号:
$$
\begin{align}
    v_i^{k + \frac{1}{2}} = \dfrac{1}{2}(v_i^{k} + v_i^{k+1}), & \quad & \delta_{t}v_i^{k + \frac{1}{2}} = \dfrac{1}{\tau}(v_i^{k+1} - v_i^k), \\
    D_tv_i^k = \dfrac{1}{\tau}(v_i^{k+1} - v_i^k), & \quad & D_{\bar{t}}v_i^k = \dfrac{1}{\tau}(v_i^k - v_i^{k-1}), \\
    D_xv_i^k = \dfrac{1}{h}(v_{i+1}^{k} - v_i^k), & \quad & D_{\bar{x}}v_i^k = \dfrac{1}{h}(v_i^k - v_{i-1}^{k}), \\
    \delta_xv_{i + \frac{1}{2}}^k = \dfrac{1}{h}(v_{i+1}^k - v_i^k), & \quad & \delta_x^2v_i^k = \dfrac{1}{h^2}(v_{i-1}^k - 2v_i^k + v_{i+1}^k ), \\
    v_i^{\bar{k}} = \dfrac{1}{2}(v_i^{k+1} + v_i^{k-1}), & \quad & D_{\hat{t}}v_i^k = \dfrac{1}{2\tau}(v_i^{k+1} - v_i^{k-1}),
\end{align}
$$
# 向前$Euler$格式

对原方程进行离散,可得到如下差分格式:
$$
\begin{align}
    D_tu_i^k - a\delta_x^2u_i^k = f(x_i, t_k), \quad & 0 \leq i \leq m-1, \ 0 \leq k \leq n-1 \\
    u_i^0 = \varphi(x_i), \quad & 0 \leq i \leq m, \\
    u_0^k = \alpha(t_k), \ u_m^{k} = \beta(t_k), \quad & 1 \leq k \leq n
\end{align}
$$
记 $ r = \dfrac{a \tau}{h^2} $称为$\color{red}{步长比}$ 差分格式写为:
$$
\begin{equation}
    u_i^{k + 1} = (1 - 2r)u_i^k +r(u_{i-1}^k + u_{i+1}^k) + \tau f(x_i, t_k), \ \ 1 \leq i \leq m-1, 0 \leq k \leq n-1
\end{equation}
$$
上式表明: 第 $k+1$ 层上的值可直接由第 $k$ 层上的值 显示表达出, 若已知第 $ k $ 层的值, 从上式可以求出第 $ k+1 $ 层的值

将上式改写为矩阵形式如下:
$$
\begin{equation} 
    \begin{pmatrix} u_1^{k+1} \\ u_2^{k+1} \\ \vdots \\ u_{m-2}^{k+1} \\ u_{m-1}^{k+1}\end{pmatrix} = 
    \begin{pmatrix} 
    1 - 2r & r &  & \\
    r & 1 - 2r & r & & \\   
    & \ddots & \ddots & \ddots &  \\
    & &  r & 1 -2r & r \\
    & & & r & 1- 2r
    \end{pmatrix}
    \begin{pmatrix} u_1^k \\ u_2^k \\ \vdots \\ u_{m-2}^k \\ u_{m-1}^k \end{pmatrix} + 
    \begin{pmatrix} 
    \tau f(x_1, t_k) + ru_0^k \\
    \tau f(x_2, t_k) \\
    \vdots \\
    \tau f(x_{m-2}, t_k) \\
    \tau f(x_{m-1}, t_k) + ru_m^k
    \end{pmatrix}
\end{equation}
$$
# 向前 $Euler$ 格式数值实验

利用向前 $ Euler $ 格式求解如下定解问题：
$$
\begin{align}
    \dfrac{\partial u}{\partial t} - \dfrac{\partial^2 u}{\partial x^2} = 0, \quad & 0 < x < 1, \ 0 < t < 1 \\
    u(x, 0) = e^x, \quad & 0 \leq x \leq 1 \\
    u(0, t) = e^t, \quad u(1, t) = e^{1+t}, \quad & 0 < t \leq 1
\end{align}
$$
该定解问题的精确解为:

$$ u(x, t) = e^{x+t} $$

## 模型参数

```python
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

class ModelData:
    def __init__(self, t, x):
        self.t0 = t[0]
        self.t1 = t[1]
        self.x0 = x[0]
        self.x1 = x[1]

    def time_grid(self, NT):
        T = np.linspace(self.t0, self.t1, NT + 1)
        tau = (self.t1 - self.t0) / NT
        return tau, T

    def space_grid(self, NS):
        X = np.linspace(self.x0, self.x1, NS + 1)
        h = (self.x1 - self.x0) / NS
        return h, X

    def solution(self, X, T):
        x, t = np.meshgrid(X, T)
        return np.exp(x+t).T

    def rhsf(self, x, t):
        return np.zeros(len(x))

    def init_solution(self, x):
        return np.exp(x)

    def left_solution(self, t):
        return np.exp(t)

    def right_solution(self, t):
        return np.exp(t + 1)

    def a(self):
        return 1
```

## 核心程序

```python
def parabolic_fdbvd(model, NT, NS, method):
    tau, T = model.time_grid(NT)
    h, X = model.space_grid(NS)
    N, M = len(T), len(X)
    uh = np.zeros((M, N))
    
    # 步长比
    r =  (model.a() * tau) / h ** 2
  
    # 初值条件
    uh[:, 0] = model.init_solution(X)

    # 左边界条件
    uh[0, :] = model.left_solution(T)

    # 右边界条件
    uh[-1, :] = model.right_solution(T)

    def forward_euler():
        assert r <= 0.5, '差分格式不收敛, 无法使用'
        c0 = np.ones(M - 2) * (1 - 2 * r)
        c1 = np.ones(M - 3) * r
        A = sparse.diags([c1, c0, c1], [-1, 0, 1], shape = (M - 2, M - 2), format = 'csr', dtype = np.float64)

        for i in range(N-1):
            bias = tau * model.rhsf(uh[1: -1, i], T[i])
            bias[0] += r * uh[0, i]
            bias[-1] += r * uh[-1, i]
            uh[1:-1, i+1] = np.dot(A.toarray(), uh[1:-1, i].reshape(-1, 1)).flatten() + bias
        
        return uh, T, X

    def backward_euler():
        pass

    if method == 'forward_euler':
        return forward_euler()

```
## 测试程序
```python
def parabolic_fdbvd_loss(solution, uh, T, X):
    emax = np.abs(solution(X, T) - uh)
    return np.max(emax), emax

def parabolic_fdbvd_test():
    t = np.array([0, 1])
    x = np.array([0, 1])
    NT = np.array([200, 800, 3200, 12800])
    NS = np.array([10, 20, 40, 80])
    emax = np.zeros(4)
    uhrec = []
    Trec = []
    Xrec = []
    error_rec = []
    
    pde = ModelData(t, x)
    uh, T, X = parabolic_fdbvd(pde, NT[0], NS[0], 'forward_euler')
    u = pde.solution(X, T)
    
    print('u = ', u, ' \n shape = ', u.shape) 
    print('')
    print('uh = ', uh, ' \n shape = ', uh.shape)
    print('')
    print('time grid = ', T)
    print('')
    print('space grid = ', X)
    
    for i in range(4):
        uh, T, X = parabolic_fdbvd(pde, NT[i], NS[i], 'forward_euler')
        uhrec.append(uh)
        Trec.append(T)
        Xrec.append(X)
        emax[i], err = parabolic_fdbvd_loss(pde.solution, uh, T, X)
        error_rec.append(err)
    print('')
        
    # 绘制误差表格
    data = {
        'NS':[NS[i] for i in range(4)],
        'NT': [NT[i] for i in range(4)],
        'Emax':[emax[i] for i in range(4)],
        'Emax(2h, 4tau) / Emax(h, tau)':['*'] + [emax[i - 1]/emax[i] for i in range(1, 4)]
    }
    
    df = pd.DataFrame(data)
    print(df)
    
    # 动画
    fig = plt.figure(figsize = (10, 5))
    ax = plt.subplot(121)
    xplt = Xrec[0]
    uhplt = uhrec[0][:, 0]
    line, = ax.plot(xplt, uhplt)
    
    def init_func():
        line.set_ydata(uhrec[0][:, 0])
        return line,
    
    def anima_func(i):
        line.set_ydata(uhrec[0][:, i])
        return line,
    
    ani = animation.FuncAnimation(fig = fig, func = anima_func, frames = 200, \
                                  init_func = init_func, interval = 10, blit = False)
    
    # 二维图
   
    ax2 = plt.subplot(122 ,projection = '3d')
    TX = np.meshgrid(Trec[0], Xrec[0])
    err = error_rec[0]
    ax2.plot_surface(TX[0], TX[1], err, cmap = 'gist_rainbow')
    plt.show()

```