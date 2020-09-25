
import numpy as np
import matplotlib.pyplot as plt

def thomas(a, b, c, d):
    n = len(b)
    x = np.zeros_like(b)
    y = np.zeros_like(b)
    beta = np.zeros_like(b)
    beta[0] = b[0]
    y[0] = d[0]
    for i in range(1, n):
        l = a[i-1] / beta[i-1]
        beta[i] = b[i] - l*c[i-1]
        y[i] = d[i] - l*y[i-1]
    x[n-1] = y[n-1]/beta[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - c[i]*x[i+1])/beta[i]
    return x
    
class Model_Data:
    def __init__(self, lbd, rbd):
        self.lbd = lbd
        self.rbd = rbd
    
    def mesh_init(self,NS):
        x = np.linspace(self.lbd, self.rbd, NS+1)
        h = (self.rbd - self.lbd) / NS
        return x, h
    
    def solution(self,x):
        return np.exp(x)*np.sin(x)
    
    def q(self,x):
        return np.ones_like(x)
    
    def f(self,x):
        return np.exp(x)*(np.sin(x) - 2*np.cos(x))

def fd1d_bvp(model, NS):
    x, h = model.mesh_init(NS)

    # lhs matrix:
    a = -1 * np.ones(NS - 2)
    b = 2 + h**2 * model.q(x[1: -1])
    c = -1 * np.ones(NS - 2)

    # rhs vector:
    rhs = h**2 * model.f(x[1: -1])
    rhs[0] = h**2 * model.f(x[1]) + model.solution(x[0]) 
    rhs[-1] = h**2 * model.f(x[-2]) + model.solution(x[-1])

    # get numeric solution on meshgrid 
    uh = np.zeros(NS + 1)
    uh[1: NS] = thomas(a, b, c, rhs)
    uh[0] = model.solution(x[0])
    uh[-1] = model.solution(x[-1])

    return x, uh

#error function
def fd1d_bvp_error(solution, uh, x):
    NN = len(x)
    h = (x[-1] - x[0])/(NN - 1)
    u = solution(x)

    ee = u - uh
    ee2 = np.abs(u - uh)
    
    e0 = h * np.sum(ee ** 2)
    e1 = np.sum( (ee[1:] - ee[0:-1])**2 )/ h
    e1 = e1 + e0
    
    e0 = np.sqrt(e0)
    e1 = np.sqrt(e1)
    e_max = np.max(np.abs(ee))
    
    return ee2, e0,e1,e_max

# test function
def fd1d_bvp_test():
    # init data
    NS = [10, 20, 40, 80, 160]
    lbd = 0
    rbd = np.pi
    
    model = Model_Data(lbd, rbd)
    
    e_max = np.zeros(5)
    e0 = np.zeros(5)
    e1 = np.zeros(5)
    ee = [0 for i in range(5)]
    X = []
    U = []
    
    # compute error
    for i in range(5):
        x,uh = fd1d_bvp(model, NS[i])
        ee[i], e0[i], e1[i], e_max[i] = fd1d_bvp_error(model.solution, uh, x)
        X.append(x)
        U.append(uh)
    
    # exact soution on last meshgrid
    u = model.solution(X[-1])
    
    # visuable solution
    fig = plt.figure(figsize=(10, 5))
    ax = plt.subplot(121)
    lss = ['-','-.', ':', '-.', ':']
    colors = ['c', 'b', 'y', 'g', 'm'] 
    markers = ['o', '+', '<', 'x', '']
    labels = ['NS = 10', 'NS = 20', 'NS = 40', 'NS = 80', 'NS = 160']
    labels2 = ['error_NS = 10', 'error_NS = 20', 'error_NS = 40', 'error_NS = 80', 'error_NS = 160']
   
    # numeric solutions on every mesh
    ax.plot(X[4], u, color = 'r', marker = 'x', label = 'exact solution')
    for j in range(2):
        ax.plot(X[j], U[j], ls = lss[j], color = colors[j], marker = markers[j], label = labels[j])
    plt.legend() 
    ax2 = plt.subplot(1, 2, 2)
    # error curves
    for k in range(5):
        ax2.plot(X[k], ee[k], color = colors[k], marker = markers[k], label = labels2[k])
    plt.legend()
    plt.show()
    
    print('emax = ', e_max)
    print('e_0 = ', e0)
    print('e_1 = ', e1)
    for i in range(len(e_max) -1):
        print(e_max[i] / e_max[i+1])
        
if __name__ == '__main__':
    fd1d_bvp_test()
