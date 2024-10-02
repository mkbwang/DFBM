import torch
import pickle


def log1exp(mymatrix: torch.Tensor):
    """
    prevent numerical issues when calculating log(1+exp(x))
    :param mymatrix:
    :return:
    """
    output = torch.zeros_like(mymatrix)
    large_filter = mymatrix > 10
    # small_filter = mymatrix < -10
    middle_filter = torch.logical_and(mymatrix <= 10, mymatrix >= -10)
    output[large_filter] = mymatrix[large_filter]
    output[middle_filter] = torch.log(1 + torch.exp(mymatrix[middle_filter]))

    return output


def loglik(X: torch.Tensor, Y: torch.Tensor):
    """
    :param X: binary matrix
    :param Y: estimated log link matrix
    :return: log likelihood
    """
    Q = 2*X-1
    intermediate = Q*Y
    loglik_mat = intermediate - log1exp(intermediate)
    return torch.sum(loglik_mat)


def neg_loglik_grad(X: torch.Tensor, Y: torch.Tensor):
    """
    :param X: binary matrix
    :param Y: estimated log link matrix
    :return: gradient of negative log likelihood over Y
    """
    Q = 2*X-1
    return -Q/(1+torch.exp(Q*Y))


def lowrank(Y: torch.Tensor, k: int):
    """
    use SVD to get a lower rank version of input matrix
    :param Y:
    :param k:
    :return:
    """
    nsample = Y.shape[0]
    nfeature = Y.shape[1]
    vec_ones = torch.ones(size=(nsample, ))
    mu = Y.mean(dim=0)
    if k == 1:
        U = torch.zeros(size=(nsample, 1))
        D = 0
        Vt = torch.zeros(size=(nfeature, ))
        output = torch.ger(vec_ones, mu)
    else:  # k > 1
        new_residual_Y = Y - torch.ger(vec_ones, mu) # center each column
        all_U, all_D, all_Vt = torch.linalg.svd(new_residual_Y, full_matrices=False)
        U = all_U[:, :(k - 1)]
        D = all_D[:(k - 1)]
        Vt = all_Vt[:(k - 1), :]
        lowrank_residual_Y = U @ torch.diag(D) @ Vt
        output = torch.ger(vec_ones, mu) + lowrank_residual_Y
    return mu, U, D, Vt, output


def logisticSVD(X: torch.Tensor, k: int, conv=1e-4, max_iter=1000):
    """
    :param X: binary matrix
    :param k: preset degrees of freedom
    :param conv: convergence criteria
    :param max_iter: maximum number of iterations
    :return: estimated SVD output (U and D and V)
    """
    dimensions = X.shape
    nsample = dimensions[0]
    nfeature = dimensions[1]
    num_entries = nsample * nfeature
    # initial values
    Q = 2 * X - 1
    full_Y = 4*Q
    mu, U, D, Vt, Y = lowrank(full_Y, k)
    ell = loglik(X, Y)
    ell_trace = torch.zeros(size=(max_iter+1,), dtype=torch.float, device="cpu")
    total_iter = 1
    ell_trace[0] = ell
    # print(f'Initial log likelihood: {ell}')
    for j in torch.arange(1, max_iter+1):
        # check if convergence is met
        full_Y = Y - 4 * neg_loglik_grad(X, Y)
        mu, U, D, Vt, Y = lowrank(full_Y, k)
        ell_trace[j] = loglik(X, Y)
        # print(f'Log likelihood after iteration {j}: {ell_trace[j]}')
        if j==max_iter or torch.abs(ell_trace[j] - ell_trace[j-1])/num_entries < conv:
            total_iter = j+1
            break
    ell_trace = ell_trace[:total_iter]
    # expected_P = torch.sigmoid(Y)
    output = {"mu":mu, "U":U, "D":D,  "Vt":Vt, "loglik":ell_trace}

    return output

