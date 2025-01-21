import os
from math import radians as rad
from cmath import exp as exp
from numpy import ndarray, zeros, pi
from scipy.sparse import csr_matrix

systemName = 'ieee14bus_fault.txt'


class systemData:
    def __init__(
        self, systemName,
        system: str = "",
        Sbase: int = 100,
        nbus: int = 0,
        bus: dict = dict(),
        nbranch: int = 0,
        branch: dict = dict(),
        ybus: ndarray = [],
        sparse_ybus=[]
    ):

        # Inicialization
        self.system = os.path.splitext(systemName)[0]
        self.Sbase = Sbase
        self.nbus = nbus
        self.bus = bus
        bus['num'] = list()
        bus['name'] = list()
        bus['type'] = list()
        bus['voltage_magnitude'] = list()
        bus['voltage_angle'] = list()
        bus['active_generation'] = list()
        bus['reactive_generation'] = list()
        bus['reactive_generation_min'] = list()
        bus['reactive_generation_max'] = list()
        bus['active_demand'] = list()
        bus['reactive_demand'] = list()
        bus['capacitor_reactor'] = list()
        self.nbranch = nbranch
        self.branch = branch
        branch['from'] = list()
        branch['to'] = list()
        branch['r_value'] = list()
        branch['x_value'] = list()
        branch['b_value'] = list()
        branch['tap'] = list()
        branch['tap_angle'] = list()
        self.ybus = ybus
        self.sparse_ybus = sparse_ybus

    def add_values(self):

        Sbase = 100

        folder = os.path.join(os.path.dirname(__file__), 'Systems')
        fname = os.path.join(folder, self.system + '.pwf')

        f = open(f'{fname}', 'r', encoding='latin-1')

        lines = f.readlines()
        f.close()

        i = 0
        while lines[i].strip() != 'FIM':
            i += 1
            if 'DBAR' in lines[i].strip():
                i += 1
                while '99999' != lines[i][0:5]:
                    if lines[i][0] != '(':
                        self.nbus += 1
                        self.bus['num'].append(int(lines[i][0:5].strip()))
                        self.bus['name'].append(lines[i][10:22].strip())
                        if lines[i][7:8] != " ":
                            self.bus['type'].append(int(lines[i][7:8].strip()))
                        else:
                            self.bus['type'].append(0)
                        self.bus['voltage_magnitude'].append(
                            int(lines[i][24:28].strip())/1000)
                        try:
                            self.bus['voltage_angle'].append(
                                float(lines[i][28:32].strip()) * pi/180)
                        except:
                            self.bus['voltage_angle'].append(0.0)
                        try:
                            self.bus['active_generation'].append(
                                float(lines[i][32:37].strip())/Sbase)
                        except:
                            self.bus['active_generation'].append(0.0)
                        try:
                            self.bus['reactive_generation'].append(
                                float(lines[i][37:42].strip())/Sbase)
                        except:
                            self.bus['reactive_generation'].append(0.0)
                        if self.bus['type'][-1] != 0:
                            try:
                                self.bus['reactive_generation_min'].append(
                                    float(lines[i][42:47].strip())/Sbase)
                            except:
                                self.bus['reactive_generation_min'].append(
                                    -99.99)
                            try:
                                self.bus['reactive_generation_max'].append(
                                    float(lines[i][47:52].strip())/Sbase)
                            except:
                                self.bus['reactive_generation_max'].append(
                                    99.99)
                        else:
                            self.bus['reactive_generation_min'].append(0.0)
                            self.bus['reactive_generation_max'].append(0.0)
                        try:
                            self.bus['active_demand'].append(
                                float(lines[i][58:63].strip())/Sbase)
                        except:
                            self.bus['active_demand'].append(0.0)
                        try:
                            self.bus['reactive_demand'].append(
                                float(lines[i][63:68].strip())/Sbase)
                        except:
                            self.bus['reactive_demand'].append(0.0)
                        try:
                            self.bus['capacitor_reactor'].append(
                                float(lines[i][68:73].strip())/Sbase)
                        except:
                            self.bus['capacitor_reactor'].append(0.0)
                    i += 1
            if 'DLIN' in lines[i].strip():
                i += 1
                while '99999' != lines[i][0:5]:
                    if lines[i][0] != '(':
                        self.nbranch += 1
                        self.branch['from'].append(int(lines[i][0:5].strip()))
                        self.branch['to'].append(int(lines[i][10:15].strip()))
                        try:
                            self.branch['r_value'].append(
                                round(float(lines[i][20:26].strip())/100, 5))
                        except:
                            self.branch['r_value'].append(0.0)
                        self.branch['x_value'].append(
                            round(float(lines[i][26:32].strip())/100, 5))
                        try:
                            self.branch['b_value'].append(
                                round(float(lines[i][32:38].strip())/Sbase, 5))
                        except:
                            self.branch['b_value'].append(0.0)
                        try:
                            self.branch['tap'].append(
                                float(lines[i][38:43].strip()))
                        except:
                            self.branch['tap'].append(0.0)
                        try:
                            self.branch['tap_angle'].append(
                                float(lines[i][53:58].strip()))
                        except:
                            self.branch['tap_angle'].append(0.0)
                    i += 1

        self.ybus = zeros(shape=[self.nbus, self.nbus], dtype='complex')
        for i in range(self.nbranch):
            j = 0
            while True:
                if self.branch['from'][i] == self.bus['num'][j]:
                    from_bus = j
                    break
                j += 1
            j = 0
            while True:
                if self.branch['to'][i] == self.bus['num'][j]:
                    to_bus = j
                    break
                j += 1
            rx_branch = complex(
                real=self.branch['r_value'][i], imag=self.branch['x_value'][i])
            b_branch = complex(real=0, imag=self.branch['b_value'][i]/2)
            if self.branch['tap'][i] == 0:
                self.ybus[from_bus, to_bus] -= 1/rx_branch
                self.ybus[to_bus, from_bus] -= 1/rx_branch
                self.ybus[from_bus, from_bus] += 1/rx_branch + b_branch
                self.ybus[to_bus, to_bus] += 1/rx_branch + b_branch
            else:
                tap = self.branch['tap'][i] * \
                    exp(1j * rad(self.branch['tap_angle'][i]))
                self.ybus[from_bus, to_bus] -= (1/rx_branch)/tap.conjugate()
                self.ybus[to_bus, from_bus] -= (1/rx_branch)/tap
                self.ybus[from_bus, from_bus] += (1/rx_branch +
                                                  b_branch) / (tap * tap.conjugate())
                self.ybus[to_bus, to_bus] += (1/rx_branch) + b_branch

        for i in range(self.nbus):
            if self.bus['capacitor_reactor'] != 0:
                self.ybus[i, i] += complex(real=0,
                                           imag=self.bus['capacitor_reactor'][i])

        self.sparse_ybus = csr_matrix(self.ybus)
