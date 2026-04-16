import subprocess

import matplotlib.pyplot as plt
import numpy as np

# temporary plot for

try:
    vstr = (
        "OpenShieldHIT "
        + subprocess.check_output(
            ["git", "describe", "--tags", "--dirty", "--always"],
            stderr=subprocess.DEVNULL,
        )
        .decode()
        .strip()
    )
except subprocess.CalledProcessError:
    vstr = "OpenShieldHIT (unknown version)"
# Load data from the file, skipping the first 4 header lines
data = np.loadtxt("tests/cases/00_minimal/NB_msh.dat", skiprows=4)
data2 = np.loadtxt("../shieldhit/test/00_minimal/NB_msh.dat", skiprows=4)

# Extracting columns
z = data[:, 2]
energy = data[:, 3]
z2 = data2[:, 0]
energy2 = data2[:, 1]

# normalize energy to the max energy found in the data
energy = energy / max(energy)
energy2 = energy2 / max(energy2)

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(z, energy, marker="o", linestyle="-", label="OpenShieldHIT")
plt.plot(z2, energy2, marker="x", linestyle="--", label="SHIELD-HIT12A")
plt.title("150 MeV protons on Water\n" + vstr, fontsize=12)
plt.xlabel("Depth (cm)")
plt.ylabel("Energy")
# position the legend in the upper left corner
plt.legend(loc="upper left")
plt.grid(True)
plt.tight_layout()
plt.show()
