import numpy as np
from call_rseq import rseq

EPSVAL = 1.e-20

class Basis():
    def __init__(self,nr):
        self.g = np.zeros((2 * nr),dtype = np.float64)
        self.gfac = np.zeros((nr),dtype = np.float64)
        self.e = 0.
        self.val = 0.
        self.slo = 0.
        self.l = 0
        self.node = 0
        self.nre = 0
        self.open = False

    def make_basis(self, a, b, emin, emax, lvalsh, node, nr, rofi, slo, vofi, val):
        self.l = lvalsh
        self.node = node
        self.slo = slo
        self.val = val
        iz = 0
        if slo == 0. :
            self.open = True

        self.e, self.g, self.gfac, self.nre, self.slo, self.val = rseq(a, b, self.e, emin, emax, self.g, self.gfac, self.l, self.node, nr, self.nre, rofi, self.slo, vofi, self.val, iz)
        self.g[:nr] = self.g[:nr] * np.sqrt(self.gfac)
        if self.open:
            print ("l = {:>2}  \"open\" slo = {:>9.6f} val = {:>9.6f} energy = {:>9.6f}".format(self.l,self.slo,self.val,self.e))
        else:
            print ("l = {:>2} \"close\" slo = {:>9.6f} val = {:>9.6f} energy = {:>9.6f}".format(self.l,self.slo,self.val,self.e))
        """
        with open ("wavefunc.dat", mode = "w") as fw_w :
            fw_w.write("# rfoi, wavefunction\n")
            for i in range (nr):
                fw_w.write("{:>13.8f}".format(rofi[i]))
                fw_w.write("{:>13.8f}".format(self.g[i]))
                fw_w.write("\n")
        """

