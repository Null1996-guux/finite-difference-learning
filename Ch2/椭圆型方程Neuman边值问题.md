# $ Neumann $ 边值问题

考虑 $ Laplace $ 方程的 $ Neumann $ 边值问题: 
\begin{align}
    \frac{ \partial ^ 2 u }{ \partial x ^ 2 } + \frac{ \partial ^ 2 u }{ \partial y ^ 2 } = 0, & \ \ (x, y) \in \Omega \\
    \frac{ \partial u  }{ \partial n  } = g(x, y), &\ \ (x, y) \in \partial \Omega
\end{align}

其中 $ \Omega = [0, 1] \times [0, 1], \ \partial{ \Omega } $ 表示区域边界。

上述问题仅当 $ \int_{ \partial \Omega } g(x, y) ds = 0 $ 时有解

## 差分格式的建立

取等距网格, 假设网格步长 $ h = \frac{ 1 }{ M } $

### 内部节点的差分离散

\begin{equation}
    u_{ i + 1, j } + u_{ i -1, j } + u_{ i, j + 1 } + u_{i, j - 1 } - 4 u_{ i, j } = 0 \quad i, j = 1, 2, \cdots, M - 1
\end{equation}

### 边界节点的差分离散

在 $ x = 0 $处，边界条件的差分格式:
\begin{equation}
    u_{-1, j} - u_{1, j} = 2hg_{0, j}
\end{equation}

$ u_{-1, j} $ 表示虚网格节点, 令第一个式子中 $ i = 0 $:
\begin{equation}
    u_{1, j} + u_{-1, j} + u_{0, j + 1} + u_{0, j - 1} - 4u_{0, j} = 0 
\end{equation}

从上面两式中消去 $ u_{-1, j} $立得二阶精度的差分格式: 
\begin{equation}
    4 u_{ 0, j } - 2 u_{ 1, j } - u_{0, j + 1 } - u_{0, j - 1} = 2hg_{0, j}, \quad 1 \leq j \leq M - 1
\end{equation}

同样的方式, 可以得到右边界, 上边界以及下边界的二阶精度差分格式, 最终整理如下:
\begin{equation}
    \begin{cases}
        4 u_{ 0, j } - 2 u_{ 1, j } - u_{ 0, j + 1 } - u_{ 0, j - 1} = 2hg_{ 0, j }, \quad 1 \leq j \leq M - 1 \\ 
        4 u_{ M, j } - 2 u_{ M - 1, j } - u_{ M, j + 1 } - u_{ M, j - 1 } = 2hg_{ M, j }, \quad 1 \leq j \leq M - 1 \\ 
        4 u_{ i, 0 } - 2 u_{ i, 1 } - u_{ i + 1, 0 } - u_{ i - 1, 0 } = 2hg_{ i, 0 }, \quad 1 \leq i \leq M - 1 \\ 
        4 u_{ i, M } - 2 u_{ 1, M - 1 } - u_{ i + 1, M } - u_{ i - 1, M } = 2hg_{ i, M }, \quad 1 \leq i \leq M - 1 \\
    \end{cases}
\end{equation}

### 角点的差分离散

对于四个角点， 以左下角点为例:分别利用内部方程以及边界条件， 可以得到:
\begin{equation}
    \begin{cases}
        u_{ 1, 0 } + u_{ -1, 0 } + u_{ 0, 1 } + u_{ 0, -1 } - 4u_{ 0, 0 } = 0 \\
        u_{ -1, 0 } - u_{ 1, 0 } = 2hg_{ 0, 0 } , \quad u_{ 0, -1 } - u_{ 0, 1 } = 2hg_{0, 0} 
    \end{cases}
\end{equation}

从上述三个方程组中消去 $ u_{ -1, 0 }, u_{ 0, -1 } $得到

