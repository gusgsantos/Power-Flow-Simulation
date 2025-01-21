from numpy import flatnonzero as find, array, exp, conj, shape, concatenate, linalg, zeros, r_, inf, angle
from scipy.sparse import hstack, vstack
from scipy.sparse.linalg import spsolve
from dSbus_dV import dSbus_dV


def newtonpf(nbus, bus, branch, sparse_ybus):

    ref = [i for i in range(nbus) if bus['type'][i] == 2]
    pv = [i for i in range(nbus) if bus['type'][i] == 1]
    pq = [i for i in range(nbus) if bus['type'][i] == 0]
    pvpq = r_[pv, pq]

    # set up indexing for updating V
    # The generator status is not taken into consideration.
    gbus = sorted(ref + pv)
    nref = len(ref)
    npv = len(pv)
    npq = len(pq)
    ngen = len(gbus)

    j1 = 0
    j2 = npv  # j1:j2 - V angle of pv buses
    j3 = j2
    j4 = j2 + npq  # j3:j4 - V angle of pq buses
    j5 = j4
    j6 = j4 + npq  # j5:j6 - V mag of pq buses

    ppopt = dict()
    ppopt['pf_tol'] = 1e-8
    ppopt['pf_max_it'] = 1e-8

    pf_it = 0

    Va = array(bus['voltage_angle'])
    Vm = array(bus['voltage_magnitude'])

    V0 = Vm * exp(1j * Va)

    pf_it = 0

    Sbus = (array(bus['active_generation']) + 1j * array(bus['reactive_generation'])
            ) - (array(bus['active_demand']) + 1j * array(bus['reactive_demand']))

    while True:

        pf_it += 1

        V = Vm * exp(1j * Va)

        I = sparse_ybus*V

        mis = (V*conj(I)-Sbus)

        F = r_[mis[pv].real,
               mis[pq].real,
               mis[pq].imag]

        # evalute F(x)
        if linalg.norm(F, inf) < ppopt['pf_tol']:
            break

        [dS_dVm, dS_dVa] = dSbus_dV(sparse_ybus, V)  # Complex Derivatives

        J11 = dS_dVa[array([pvpq]).T, pvpq].real
        J12 = dS_dVm[array([pvpq]).T, pq].real
        J21 = dS_dVa[array([pq]).T, pvpq].imag
        J22 = dS_dVm[array([pq]).T, pq].imag

        J = vstack([
            hstack([J11, J12]),
            hstack([J21, J22])
        ], format="csr")

        x = concatenate((Va[pvpq], Vm[pq]), axis=0)

        dx = spsolve(J, F)

        x -= dx
        x = x.reshape(x.size,)

        Va[pv] = x[j1:j2]
        Va[pq] = x[j3:j4]
        Vm[pq] = x[j5:j6]

    return V, Va, Vm
