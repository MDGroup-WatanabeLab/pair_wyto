# plot hamiltonian for various timestep

import matplotlib.pyplot as plt

LEGEND = ["dt = 1.0fs", "dt = 0.5fs", "dt = 0.1fs"]

# read log.lammps
with open("log.lammps", "r") as f:
    lines = f.readlines()

# get tables for plot of hamiltonian-time graph
tables = []
data_start_index = -1
target_index = 2
number_of_MDsim = 0
while target_index < len(lines):
    if "Per MPI rank memory allocation" in lines[target_index - 2]:
        columns = lines[target_index - 1].split()
        time_index = columns.index("Time")
        econserve_index = columns.index("Econserve")
        table = []
        number_of_MDsim = 1
        while "Loop time of" not in lines[target_index]:
            raw_record = list(map(float, lines[target_index].split()))
            time = raw_record[time_index]
            econserve = raw_record[econserve_index]
            table.append((time, econserve))
            target_index += 1
        tables.append(table)
    else:
        target_index += 1

# plot
fig = plt.figure(figsize=(16.0, 9.0))
ax_nvt = fig.add_subplot(1, 1, 1)
x_nvt_1, y_nvt_1 = zip(*tables[0])
x_nvt_05, y_nvt_05 = zip(*tables[1])
x_nvt_01, y_nvt_01 = zip(*tables[2])
ax_nvt.plot(x_nvt_1, y_nvt_1)
ax_nvt.plot(x_nvt_05, y_nvt_05)
ax_nvt.plot(x_nvt_01, y_nvt_01)
ax_nvt.legend(LEGEND)
ax_nvt.set_title("NVT", fontdict={'fontsize':30, 'fontweight':'bold'})
ax_nvt.set_xlabel("Time[ps]", fontdict={'fontsize':20, 'fontweight':'bold'})
ax_nvt.set_ylabel("Econserve[eV]", fontdict={'fontsize':20, 'fontweight':'bold'})

plt.tight_layout()
plt.savefig("hamiltonian.png")