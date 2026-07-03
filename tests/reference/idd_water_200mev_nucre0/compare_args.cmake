# Tolerances for the auto-compared idd.dat curve (col 4 = total Dose).  Provisional
# while issue #212 (p+A elastic) is open: the all-particle dose is exactly what we
# are scrutinising, so treat any committed reference as advisory until the physics
# lands.
set(COMPARE_ARGS --tol-integral 0.05 --tol-bin 0.05 --tol-peak-mm 0.5)
