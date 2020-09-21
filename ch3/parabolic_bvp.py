import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from matplotlib import animation
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

# 误差程序
def parabolic_fdbvd_loss(solution, uh, T, X):
    emax = np.abs(solution(X, T) - uh)
    return np.max(emax)

# 向前euluer法测试程序
def parabolic_fdbvd_test():
    print("向前Euler法测试程序")
    print("")
    t = np.array([0, 1])
    x = np.array([0, 1])
    NT = np.array([200, 800, 3200, 12800])
    NS = np.array([10, 20, 40, 80])
    emax = np.zeros(4)
    uhrec = []
    Trec = []
    Xrec = []
    
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
        emax[i] = parabolic_fdbvd_loss(pde.solution, uh, T, X)
    print('')
        
    # 误差
    data = {
        'NS':[NS[i] for i in range(4)],
        'NT': [NT[i] for i in range(4)],
        'Emax':[emax[i] for i in range(4)],
        'Emax(2h, 4tau) / Emax(h, tau)':['*'] + [emax[i - 1]/emax[i] for i in range(1, 4)]
    }
    
    df = pd.DataFrame(data)
    print(df)

    fig, ax = plt.subplots()
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
    plt.show() 

if __name__ == '__main__':
    parabolic_fdbvd_test()
