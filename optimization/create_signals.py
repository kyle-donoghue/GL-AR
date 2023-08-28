import networkx as nx
import numpy as np

def graph_synth(n, p = .2):
    er_graph = nx.fast_gnp_random_graph(n, p)
    while not nx.is_connected(er_graph):
        er_graph = nx.fast_gnp_random_graph(n, p)
    return nx.adjacency_matrix(er_graph).toarray()


def varm_filter(D,x):
    p = D.shape[0]-1
    n = x.shape[0]
    m = x.shape[1]
    y = np.zeros((n,m))
    for i in range(p):
        y[i,:] = x[i,:]
        for k in range(1,i+1):
            y[i,:] = y[i,:] + np.dot(D[k,:,:],y[i-k+1,:].T).T
    for i in range(p,n):
        y[i,:] = x[i,:]
        for k in range(1,p+1):
            y[i,:] = y[i,:] + np.dot(D[k,:,:],y[i-k+1,:].T).T
    return y

def custom_AR_fit(x, P):
    N = len(x)
    r = np.correlate(x,x,mode='full')
    r = r[N-1:]/N
    A_r = np.zeros(P+1)
    A_r[0] = 1
    eps = np.zeros(P+1)
    eps[0] = r[0]
    alpha = np.zeros(P+1)
    alpha[0] = r[1]
    k = np.zeros(P+1)
    k[0] = -alpha[0]/eps[0]
    for i in range(1,P+1):
        A_r = A_r + np.concatenate(([0],k[i-1]*np.flipud(A_r[0:i]),np.zeros(P-i)))
        eps[i] = (1-k[i-1]**2)*eps[i-1]
        alpha[i] = r[i+1]
        for ii in range(1,i+1):
            alpha[i] = alpha[i] + r[(i+1)-ii]*A_r[ii]
        k[i] = -alpha[i]/eps[i]
    return A_r#, eps

def create_signals_AR(n, P, interval, inov, N, Fs, freqs, noise):
    D = np.zeros((P+1,n,n))
    D[0,:,:] = np.eye(n)
    D[P,:,:] = graph_synth(n)
    # print(D)
    # make D directed to remove unstableness
    D[P,:,:] = np.triu(D[P,:,:])
    D[P,:,:] = [[0,1,0,0,0],[0,0,1,0,0],[0],[],[]]
    x = np.zeros((N,n))
    #x is a sine wave
    for i in range(inov):
        for m in range(n):
            x[i,m] = np.sin(i*2*np.pi*freqs[m]/Fs)
    x = x + np.random.normal(0,noise,(N,n))
    y = varm_filter(D,x)
    coeffs = np.zeros((int(N/interval),P+1,n))
    for k in range(int(N/interval)):
        tmp = y[k*interval:(k+1)*interval,:]
        # errors = np.zeros((P+1,n))
        for i in range(n):
            coeffs[k,:,i] = custom_AR_fit(tmp[:,i], P)
    return coeffs, y, D[-1,:,:]#, errors
