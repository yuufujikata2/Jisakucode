import numpy as np
from call_rseq import rseq

EPSVAL = 1.e-20

def make_basis(LMAX,node_open,node_close,nr,a,b,rofi,vofi):

    all_basis = []

    for lvalsh in range (LMAX):
        l_basis = []
        # for close channel
        val = 0.
        slo = -1.
        for nc   in range(node_close):
            basis = Basis(nr)
            emin = -10.
            emax = 1000.
            basis.solve_rseq(a,b,emin,emax,lvalsh,nc,nr,rofi,slo,vofi,val)
            l_basis.append(basis)
        # for open channel
        val = 1.
        slo = 0.
        for no in range(node_open):
            basis = Basis(nr)
            emin = -10.
            emax = 1000.
            basis.solve_rseq(a,b,emin,emax,lvalsh,no,nr,rofi,slo,vofi,val)
            l_basis.append(basis)



        all_basis.append(l_basis)

    with open ("wavefunc_2.dat", mode = "w") as fw_w :
        fw_w.write("#r,   l, node, open or close = ")
        for l_basis in all_basis:
            for nyu_basis in l_basis:
                fw_w.write(str(nyu_basis.l))
                fw_w.write(str(nyu_basis.node))
                if nyu_basis.open:
                    fw_w.write("open")
                else:
                    fw_w.write("close")
                fw_w.write("    ")
        fw_w.write("\n")

        for i in range(nr):
            fw_w.write("{:>13.8f}".format(rofi[i]))
            for l_basis in all_basis:
                for nyu_basis in l_basis:
                    fw_w.write("{:>13.8f}".format(nyu_basis.g[i]))
            fw_w.write("\n")

    return all_basis

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

    def solve_rseq(self, a, b, emin, emax, lvalsh, node, nr, rofi, slo, vofi, val):
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

