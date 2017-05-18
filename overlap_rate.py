import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from itertools import combinations
from scipy.stats import multivariate_normal

class BiGauss(object):
    """docstring for BiGauss"""
    def __init__(self, mu1, mu2, Sigma1, Sigma2, pi1, pi2, steps = 100):
        super(BiGauss, self).__init__()
        self.mu1      = mu1
        self.mu2      = mu2
        self.Sigma1   = Sigma1
        self.Sigma2   = Sigma2
        self.pi1      = pi1
        self.pi2      = pi2
        self.biGauss1 = multivariate_normal(mean = self.mu1, cov = self.Sigma1, allow_singular = True)
        self.biGauss2 = multivariate_normal(mean = self.mu2, cov = self.Sigma2, allow_singular = True)
        self.steps    = steps
        self.inv_Sig1 = -inv(self.Sigma1)
        self.inv_Sig2 = -inv(self.Sigma2)

        # variables to calculate RC
        self.A_1 = self.inv_Sig1[0][0]
        self.B_1 = self.inv_Sig1[0][1]
        self.C_1 = self.inv_Sig1[1][0]
        self.D_1 = self.inv_Sig1[1][1]
        self.A_2 = self.inv_Sig2[0][0]
        self.B_2 = self.inv_Sig2[0][1]
        self.C_2 = self.inv_Sig2[1][0]
        self.D_2 = self.inv_Sig2[1][1]

    def pdf(self, x):
        return self.pi1 * self.biGauss1.pdf(x) + self.pi2 * self.biGauss2.pdf(x)

    def OLR(self):
        e      = math.sqrt((self.mu1[0] - self.mu2[0])**2 + (self.mu1[1] - self.mu2[1])**2) / float(self.steps)
        x_step = e*(self.mu1[0]-self.mu2[0]) # each step for x
        y_step = e*(self.mu1[1]-self.mu2[1]) # each step for y
        p_x    = self.mu1[0] - x_step

        while self.RC(p_x) == None:
            p_x = p_x - x_step
        
        p_y   = self.RC(p_x)
        p     = [p_x, p_y]
        p_pre = self.mu1
        p_min = min(self.pdf(p), self.pdf(p_pre))
        p_max = max(self.pdf(p), self.pdf(p_pre))
        index = 0
        while index < self.steps:
            if self.RC(p[0] - x_step) != None:
                p_next = [p[0] - x_step, self.RC(p[0] - x_step)] # next point on ridge curve
                if self.pdf(p) > self.pdf(p_pre) and self.pdf(p) > self.pdf(p_next):
                    p_max = self.pdf(p)
                if self.pdf(p) < self.pdf(p_pre) and self.pdf(p) < self.pdf(p_next):
                    p_min = self.pdf(p)
            p_pre = p
            p     = p_next
            index += 1

        pdf_mu1 = self.pdf(self.mu1)
        pdf_mu2 = self.pdf(self.mu2)
        return p_min / min(pdf_mu1, pdf_mu2) if p_min < min(pdf_mu1, pdf_mu2) else 1.0

    # get y given x, satisfying (x,y) is on ridge curve (RC)
    def RC(self, x):
        E = self.A_1 * (x - self.mu1[0])
        F = self.C_1 * (x - self.mu1[0])
        G = self.A_2 * (x - self.mu2[0])
        H = self.C_2 * (x - self.mu2[0])

        I = E * self.D_2 - F * self.B_2
        J = H * self.B_1 - G * self.D_1
        K = self.B_1 * self.D_2 - self.B_2 * self.D_1
        M = F * G - E * H

        P = K
        Q = I + J - K * (self.mu2[1] + self.mu1[1])
        S = -(M + I * self.mu2[1] + J * self.mu1[1])

        if Q**2 - 4*P*S < 0:
            return None

        y = max((-Q + math.sqrt(Q**2 - 4*P*S)) / (2*P), (-Q - math.sqrt(Q**2 - 4*P*S)) / (2*P))

        return y



if __name__ == '__main__':

    # Single sample
    # you can change whatever data you like
    pi1 = 0.5
    pi2 = 1.0 - pi1
    mu1 = [0,0]
    mu2 = [3,0]
    Sigma1 = [[1,0],[0,1]]
    Sigma2 = [[2.17,1.82], [1.82,2.17]]
    Bi = BiGauss(mu1, mu2, Sigma1, Sigma2, pi1, pi2)
    print('the OLR of this mixture Gaussian: {}'.format(Bi.OLR()))


    # OLR's changing over distance between two means
    step = 0.1
    pi1 = 0.5
    OLRs = []
    x = []
    mu1 = [0,0]
    mu2 = [0,0]
    Sigma1 = [[1,0],[0,1]]
    Sigma2 = [[2.17,1.82], [1.82,2.17]]
    while mu2[0] <= 8.0:
        pi2 = 1.0-pi1
        Bi = BiGauss(mu1, mu2, Sigma1, Sigma2, pi1, pi2)
        OLRs.append(Bi.OLR())
        x.append(mu2[0])
        mu2[0] += step
    plt.plot(x, OLRs)
    plt.xlabel('Distance between two means')
    plt.show()


    # OLR's changing over pi1
    step = 0.05
    pi1 = 0.0
    OLRs = []
    x = []
    mu1 = [0,0]
    mu2 = [3,0]
    Sigma1 = [[1,0],[0,1]]
    Sigma2 = [[2.17,1.82], [1.82,2.17]]
    while pi1 <= 1.0:
        pi2 = 1.0-pi1
        Bi = BiGauss(mu1, mu2, Sigma1, Sigma2, pi1, pi2)
        OLRs.append(Bi.OLR())
        x.append(pi1)
        pi1 += step
    print(OLRs)
    plt.plot(x, OLRs)
    plt.xlabel('Coefficient(%) of {}'.format(r"$\pi_1$"))
    plt.show()
