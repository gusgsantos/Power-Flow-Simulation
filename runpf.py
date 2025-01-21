# Developed by Gustavo Santos
# e-mail: g.gustavo.santos@gmail.com

from system_Data import systemData
from newtonpf import newtonpf


systemName = '39bus_system.pwf'
systemName = 'ieee14.pwf'

# Reading system data
d = systemData(systemName)
d.add_values()

# Power flow
V, Vm, Va = newtonpf(d.nbus, d.bus, d.branch, d.sparse_ybus)

teste = 1
