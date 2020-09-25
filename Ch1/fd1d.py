import numpy as np
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
    

    
    
