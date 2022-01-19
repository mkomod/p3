import torch
import torch.nn as nn
import torch.distributions as dist


class NeuralPCA(nn.Module):
    """
        Basic Neural PCA
    """
    def __init__(self, input_dim, encoding_dim):
        super(NeuralPCA, self).__init__()
        self.l1 = nn.Linear(input_dim, encoding_dim)
        self.l2 = nn.Linear(encoding_dim, input_dim)

    def forward(self, x):
        x = self.l1(x)
        x = self.l2(x)
        return x

    def subspace(self, x):
        x = self.l1(x)
        return x



def loss(batch, net):
    J = 0.0

    for x in batch:
        J += torch.sum(torch.pow(x - net(x), 2.0))
    
    return J


if __name__=="__main__":
    torch.manual_seed(1)

    X_dist = dist.MultivariateNormal(torch.zeros(5), torch.eye(5))
    X = X_dist.sample((50, ))
    
    net = NeuralPCA(5, 4)
    opt = torch.optim.SGD(params=net.parameters(), lr=1e-3)

    for i in range(1, 100):
        J = loss(X, net)

        opt.zero_grad()
        J.backward()
        opt.step()

        if i % 50 == 0:
            print(J)
    
    U, S, V = torch.svd(torch.cov(X.T))
    X_proj = torch.stack([torch.sum(i * U[:, 0:4], dim=1) for i in X @ U[:, 0:4]])

    print(X[1, ])
    print(net(X[1, ]))
    print(X_proj[1, ])
    print(torch.sum(torch.pow(X - X_proj, 2.0)))


