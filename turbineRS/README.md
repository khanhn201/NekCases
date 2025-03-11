# TurbineRS

# Inner vs Outer
- Different usrdat2 to scale and project the boundaries and reassign boundary type
- ~Different usrbc~
- ~Different usrcheck (inner have rotation and torqz, can move to udf?)~


# Boundaries ID
- 'mv ': 1
- 'W  ': 2
- 'v  ': 3
- 'O  ': 4
- 'int': 5

# Run
`nohup mpiexec -np 2 nekrs --setup turbine.sess > logfile 2>&1 &`
