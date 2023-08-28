import create_signals
from cvxopt import matrix, solvers
import matplotlib.pyplot as plt
import numpy as np
import sklearn.metrics as metrics
from scipy import sparse

def solve_Y(L, X, alpha):
    Y = np.dot(np.linalg.inv((np.eye(L.shape[0]) - alpha*L)),X)
    return Y

def solve_L(Y, N, alpha, beta, A, G, b, h, M):
    #use cvxopt to solve for vech(L)
    Y_hat = np.ravel(np.dot(Y,Y.T))
    P = 2*beta*np.dot(M.T,M)
    q = alpha*np.dot(Y_hat.T,M)
    solvers.options['show_progress'] = False
    sol = solvers.qp(matrix(P), matrix(q), matrix(G), matrix(h), matrix(A), matrix(b))
    vech_L = sol['x']
    return np.reshape(np.dot(M,vech_L), [N, N]).T

def get_u_vec(i, j, n):
    u_vec = np.zeros(n*(n+1)//2)
    pos = (j-1) * n + i - j*(j-1)//2
    u_vec[pos-1] = 1
    return u_vec
def get_T_mat(i, j, n):
    Tij_mat = np.zeros((n, n))
    Tij_mat[i-1, j-1] = Tij_mat[j-1, i-1] = 1
    return np.ravel(Tij_mat)
def create_dup_matrix(num_vertices):
    M_mat = np.zeros((num_vertices**2, num_vertices*(num_vertices + 1)//2))
    # tmp_mat = np.arange(num_vertices**2).reshape(num_vertices, num_vertices)
    for j in range(1, num_vertices+1):
        for i in range(j, num_vertices+1):
            u_vec = get_u_vec(i, j, num_vertices)
            Tij = get_T_mat(i, j, num_vertices)
            # pdb.set_trace()
            M_mat += np.outer(u_vec, Tij).T

    return M_mat

def create_eq_ineq_mat(N):
    # X = ones(N);
    # [r,c] = size(X);
    # i     = 1:numel(X);
    # j     = repmat(1:c,r,1);
    # B     = sparse(i',j(:),X(:))';
    # mat_cons1 = B*mat_obj;
    # create python code from above matlab code
    X = np.ones((N, N))
    r, c = X.shape
    i = np.arange(X.size)
    j = np.tile(np.arange(c), (r, 1))
    B = np.zeros((X.size, N))
    B[i, j.flatten('F')] = X.flatten()
    mat_cons1 = np.dot(B.T, create_dup_matrix(N))
    # for i = 1:N
    #     tmp{i} = ones(1,N+1-i);
    #     tmp{i}(1) = 0;
    # end
    # mat_cons2 = spdiags(horzcat(tmp{:})',0,N*(N+1)/2,N*(N+1)/2);
    # create python code from above matlab code
    tmp = []
    for i in range(N):
        tmp.append(np.ones(N-i))
        tmp[i][0] = 0
    mat_cons2 = np.diag(np.hstack(tmp))
    # vec_cons3 = sparse(ones(1,N*(N+1)/2)-horzcat(tmp{:}));
    # create python code from above matlab code
    vec_cons3 = np.ones(N*(N+1)//2) - np.hstack(tmp)
    # A1 = [mat_cons1;vec_cons3];
    # A2 = [mat_cons2];
    # create python code from above matlab code
    A1 = np.vstack((mat_cons1, vec_cons3))
    A2 = mat_cons2
    return A1, A2

def create_eq_ineq_vec(N):
    # b1 = [sparse(N,1);N];
    # b2 = sparse(N*(N+1)/2,1);
    # create python code from above matlab code
    b1 = np.vstack((np.zeros((N, 1)), N))
    b2 = np.zeros((N*(N+1)//2, 1))
    return b1, b2

freqs = [27, 40, 55, 59, 76]
Fs = 200
node_count = 5
[approx, y, adj] = create_signals.create_signals_AR(node_count, 15, 50, 100, 1000, Fs, freqs, .05) # 5 nodes, 15 time lags, 100 samples of inovation, 500 samples total
#innovation samples is the length of stimuli, node_count is the number of nodes, P is the number of time lags, N is the total number of samples, interval is the number of samples in each set
iter = 10
A, G = create_eq_ineq_mat(5)
b, h = create_eq_ineq_vec(5)
M = create_dup_matrix(5)
ratio = np.logspace(-3, 3, 50)
# print(ratio)
mse = np.zeros(len(ratio))
# for k, r in enumerate(ratio):
#     print(k, end='\r')
#     mean_adj = np.zeros([N,N])
#     for i in range(approx.shape[0]):
#         Y = approx[i,:,:].T
#         alpha = 1
#         beta = 1/r
#         for i in range(iter):
#             L = solve_L(Y, N, alpha, beta, A, G, b, h, M)
#             Y = solve_Y(L, approx[i,:,:].T, alpha)
#         opt_adj = np.diag(np.diag(L))-L
#         opt_adj = .5 * (opt_adj > .2) *np.ones((N,N))
#         mean_adj[:,:] += opt_adj/approx.shape[0]
#     mse[k] = np.mean((mean_adj - adj)**2)/.5
# print(mse)
# print(ratio[np.argmin(mse)])


np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
ratio = 10
mean_adj = np.zeros([node_count,node_count])
for i in range(approx.shape[0]):
    Y = approx[i,:,:].T
    alpha = 1
    beta = 1/ratio
    for i in range(iter):
        L = solve_L(Y, node_count, alpha, beta, A, G, b, h, M)
        Y = solve_Y(L, approx[i,:,:].T, alpha)
    opt_adj = np.diag(np.diag(L))-L
    opt_adj = .5 * (opt_adj > .2) *np.ones((node_count,node_count))
    mean_adj[:,:] += opt_adj/approx.shape[0]
print(mean_adj)
print()
print(adj)

plt.plot(y)
plt.show()
# print(metrics.f1_score(np.ravel(adj), np.ravel(mean_adj), average='binary'))
# print(metrics.f1_score(np.ravel(adj), np.ravel(mean_adj), average='micro'))
# print(metrics.f1_score(np.ravel(adj), np.ravel(mean_adj), average='macro'))
# print(metrics.f1_score(np.ravel(adj), np.ravel(mean_adj), average='weighted'))