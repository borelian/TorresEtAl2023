# %%
from steelcc import *
from steelccGurobi import *
import sys



if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: run.py file_name useMill useGrade addCritPair useStarting")
        sys.exit(-1)

    sol = SteelccSolution()
    print ("###### LOADING PROBLEM ###########")
    sol.loadSchoppyOutput(sys.argv[1])
    print(sol.getCosts())
    sol.validateSolution()
    print ("###### FORMULATING PROBLEM ###########")
    formulateModel(sol,
     useMillCuts = int(sys.argv[2]),
     useGradeCuts = int(sys.argv[3]),
     addCriticalPairs = int(sys.argv[4]))
    if int(sys.argv[5]):
      print ("###### INITIAL SOLUTION ###########")
      useStartingSolution(sol)
    print ("###### OPTIMIZING PROBLEM ###########")
    optimize(sol,21600)
    print(sol.getCosts())
    sol.validateSolution()



