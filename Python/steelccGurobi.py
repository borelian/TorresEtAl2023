import gurobipy as gp
from gurobipy import GRB
from steelcc import *


def Callback(model, where):
    tolerance = 1e-6
    MaxOrderWidth = model._maxOrderW
    
    # Lazy cuts for Mill edge cuts    
    if (model._useLazyCuts) & (where == GRB.Callback.MIPSOL):
        varX = model.cbGetSolution(model._varX)
        varV = model.cbGetSolution(model._varV)
        varB = model.cbGetSolution(model._varB)
        for c in model._coils:
            if model._coilsEdges[c] == 'M':
                for k in model._casters:
                    for s in model._slots:
                        if varV[k,s] + tolerance < model._coilsWidth[c] * varX[k,s,c] + varB[c]*MIN_TRIMLOSS_MILL_ALTERNATIVE:
                            model.cbLazy(model._varV[k,s] >= model._coilsWidth[c] * model._varX[k,s,c] + model._varB[c]*MIN_TRIMLOSS_MILL_ALTERNATIVE)
                        if varV[k,s] - tolerance > MAX_TRIMLOSS_MILL + model._coilsWidth[c] * varX[k,s,c] + varB[c]*(MAX_TRIMLOSS_MILL_ALTERNATIVE - MAX_TRIMLOSS_MILL) + (1.0 - varX[k,s,c])*MaxOrderWidth:
                            model.cbLazy(model._varV[k,s] <= MAX_TRIMLOSS_MILL + model._coilsWidth[c] * model._varX[k,s,c] + model._varB[c]*(MAX_TRIMLOSS_MILL_ALTERNATIVE - MAX_TRIMLOSS_MILL) + (1.0 - model._varX[k,s,c])*MaxOrderWidth)



def formulateModel(sol, useMillCuts = 1, useGradeCuts = 4,  addCriticalPairs = 1):

    # useMillCuts
    #   0 : ignore special mill cut
    #   1 : enforce special mill cut as constraints
    #   2 : enforce special mill cut as lazy contraints
    # 
    # useGradeCuts
    #   0 : grade change at slot variables (required for heats)
    #   1 : (0) + with subset grades cuts to improve LB
    #   2 : grade change between grades variable + min grade change per g
    #   3 : (2) + improved grade change per caster
    #   4 : (3) + with subset of grades
    #   Note: previously (4) was (3) and (4) does not exists

    # addCriticalPairs
    #   0 : no cuts
    #   1 : Cuts of critical pairs forbiding widths too low/high at each slot
    #   2 : Cuts of critical pairs forbidding pairs of coils (|C|^2 constraints)

    #Formulate the model in Python, 

    #Sets the number of casters to be used
    sol.numCasters = sol.numCasters

    #Sets the number of slots to be used
    sol.numSlots = int(len(sol.coils)/sol.numCasters)

    #Computes critic pairs for cuts
    sol.computeCriticPairs()

    #Creates the model on Gurobi
    m = gp.Model("SteelccModel")
    m.Params.Threads = 2
    
    #Creates parameters for constraints
    casters = list(range(1,sol.numCasters+1))
    slots = range(1,sol.numSlots+1)
    coils = list(sol.coils)
    grades = list(sol.getGrades())
    
    #Computes maximum values for degradation and OrderWidth
    MaxOrderWidth = max(sol.coils[c].orderWidth for c in coils)

    #Defining global variables for lazy-constraints function (callbacks)
    m._casters = casters
    m._coils = coils
    m._slots = slots
    m._coilsEdges = {x : sol.coils[x].edge for x in coils}
    m._coilsWidth = {x : sol.coils[x].orderWidth for x in coils}
    m._maxOrderW = MaxOrderWidth
    m._criticPairs = sol.criticPairs
    m._sortedDecWidth = sol.sortedDecWidth
    m._sortedIncWidthUB = sol.sortedIncWidthUB
    m._numSlots = sol.numSlots
    

    # grade changes
    gradeChanges = {}
    for g1 in grades:
        for g2 in grades:
            if g1 != g2 :
                gradeChanges[g1,g2] = getGradeCost(g1,g2)
    
    m._GChgKeys = gradeChanges.keys()
    
    #Defining variables
    V = m.addVars(casters, slots, name="CastWidth")
    X = m.addVars(casters, slots, coils, vtype=GRB.BINARY, name="Assignation")


    #Adding constraints

# 1) Caster width constraints:
    P = m.addVars(slots, ub=MAX_EXTRA_INTERCASTER_WIDTH, name = "WidthDiff")

    m.addConstrs( (V[k,s] >= V[k,s+1] for k in casters for s in range(1,sol.numSlots)), name = "DecreasingWidth") # Decreasing width
    m.addConstrs( (V[k,s+1] >= V[k,s] - MAX_CASTER_DECREASE for k in casters for s in range(1,sol.numSlots)), name = "CasterDecrease") # Avoiding large decreases in width
    if len(casters)==2:
        m.addConstrs( (V[1,s] - V[2,s] <= ALLOWABLE_INTERCASTER_WIDTH + P[s] for s in slots), name = "IntercasterDiff1") # Inter-caster excess diff
        m.addConstrs( (V[2,s] - V[1,s] <= ALLOWABLE_INTERCASTER_WIDTH + P[s] for s in slots), name = "IntercasterDiff2") # Inter-caster excess diff

    # Maximum diff between casters set as an upper bound for the variable P

    #Cuts to improve LB
    m.addConstrs( (V[k,s] >= gp.quicksum(sol.coils[c].orderWidth * X[k,s,c] for c in coils) for k in casters for s in slots), name = "CasterWidthLB") #Lower bound on caster width
    m.addConstrs( (V[k,s] <= gp.quicksum(sol.coils[c].getWidthUB() * X[k,s,c] for c in coils) for k in casters for s in slots), name = "CasterWidthUB") #Upper bound on caster width

# 4) Roller constraints

    Ybar = m.addVars(slots, vtype=GRB.BINARY, name="RollerChange")
    CumDegrad = m.addVars(slots, ub=MAX_ROLLER_DEGRADATION, name="CummulativeDegradation")

    m.addConstr(Ybar[1]==1)  # Initial roller change
    # min Gauge for roller change
    m.addConstrs(Ybar[s] <= gp.quicksum( X[k,s,c] for c in coils if sol.coils[c].gauge>=MIN_GAUGE_ROLLER_CHANGE) for k in casters for s in slots)

    m.addConstr(CumDegrad[1] == gp.quicksum(sol.coils[c].degradation*X[k,1,c] for k in casters for c in coils))

    # minimum degradation (useful if roller change)
    m.addConstrs(CumDegrad[s] >= gp.quicksum( sol.coils[c].degradation * X[k,s,c] for k in casters for c in coils) for s in slots)  

    # Cumulative degradation (unless roller change)
    m.addConstrs(CumDegrad[s] >= CumDegrad[s-1] + gp.quicksum( sol.coils[c].degradation * X[k,s,c] for k in casters for c in coils) - MAX_ROLLER_DEGRADATION*Ybar[s] for s in slots if s>1)

    # min Roller change (to improve LB)
    nRollers = math.ceil(sum(sol.coils[c].degradation/MAX_ROLLER_DEGRADATION for c in coils))
    m.addConstr( gp.quicksum(Ybar[s] for s in slots) >= nRollers) 


# 6) Valid assignments
    m.addConstrs( (gp.quicksum(X[k,s,c] for k in casters for s in slots) == 1 for c in coils), name="ValidAssign1")
    m.addConstrs( (gp.quicksum(X[k,s,c] for c in coils) == 1 for k in casters for s in slots), name="ValidAssign2")

# 7) gauge penalization (c-1s are used in tunder in order to get values from tunder list indexed in 0 first)
    GaugeJump = m.addVars(casters, slots, name="GaugeJump")
    m.addConstrs( GaugeJump[k,s] >= 25000* (gp.quicksum(sol.coils[c].getGaugeMin()*X[k,s,c] for c in coils) - gp.quicksum(sol.coils[cc].gauge*X[k,s+1,cc] for cc in coils))  for k in casters for s in range(1,sol.numSlots))

# 7) grade penalization
    GradeChange = m.addVars(casters, slots,  vtype=GRB.BINARY,name="GradeChange")

    m.addConstrs( GradeChange[k,s] >= gp.quicksum(X[k,s,c] for c in coils if sol.coils[c].grade == g) + gp.quicksum(X[k,s+1,cc] for cc in coils if sol.coils[cc].grade != g) - 1 for k in casters for s in range(1,sol.numSlots) for g in grades)

    if useGradeCuts < 2:
        GradeChangeCost = m.addVars(casters, slots, lb=0, name="GradeChangeCost")
        m.addConstrs( GradeChangeCost[k,s] >= gradeChanges[g1,g2] * (gp.quicksum(X[k,s,c] for c in coils if sol.coils[c].grade == g1) + gp.quicksum(X[k,s+1,cc] for cc in coils if sol.coils[cc].grade == g2) - 1) for k in casters for s in range(1,sol.numSlots) for (g1,g2) in gradeChanges.keys())

    if useGradeCuts >= 2:
        GChg = m.addVars(casters, slots, gradeChanges.keys(), vtype=GRB.BINARY, name="GradeChangeSep")
        # define variable
        m.addConstrs( GChg[k,s,g1,g2] >= gp.quicksum(X[k,s,c] for c in coils if sol.coils[c].grade == g1) + gp.quicksum(X[k,s+1,cc] for cc in coils if sol.coils[cc].grade == g2) - 1 for k in casters for s in range(1,sol.numSlots) for (g1,g2) in gradeChanges.keys())
        # redundancy (to eliminate previous variable)
        m.addConstrs( GradeChange[k,s] == gp.quicksum(GChg[k,s,g1,g2] for (g1,g2) in gradeChanges.keys()) for k in casters for s in range(1,sol.numSlots))
 



# 9) Mill edge cut contraints. Can be formulated as Lazy
    m._useLazyCuts = False
    m.Params.LazyConstraints = 0
    if useMillCuts > 0:
        B = m.addVars(coils, vtype=GRB.BINARY, name="AlternativeMillWidth")
        if useMillCuts == 1:
            m.addConstrs( V[k,s] >= sol.coils[c].orderWidth * X[k,s,c] + B[c]*MIN_TRIMLOSS_MILL_ALTERNATIVE for k in casters for s in slots for c in coils if sol.coils[c].edge == 'M' )
            m.addConstrs( V[k,s] <= MAX_TRIMLOSS_MILL + sol.coils[c].orderWidth * X[k,s,c] + B[c]*(MAX_TRIMLOSS_MILL_ALTERNATIVE - MAX_TRIMLOSS_MILL) + (1 - X[k,s,c])*MaxOrderWidth for k in casters for s in slots for c in coils if sol.coils[c].edge == 'M' )
        else:
            m.Params.LazyConstraints = 1
            m._useLazyCuts = True


            
# 11) Heats constraints
    useHeatChg = True
    sortedW = np.sort([sol.coils[c].weight * (1/2000) for c in sol.coils.keys()])
    MAX_NUM_HEAT_CHANGES = np.ceil(
            np.sum(sortedW[::-1][:sol.numSlots]) / MIN_HEAT_WEIGHT) -1
    MIN_NUM_HEAT_CHANGES = np.ceil(
        np.sum(sortedW[:sol.numSlots]) / MAX_HEAT_WEIGHT) -1

    if useHeatChg:
        HeatChg = m.addVars(casters, slots, vtype=GRB.BINARY, name="HeatChange")
        CumHeat = m.addVars(casters, slots, lb=0, ub=MAX_HEAT_WEIGHT, name="cumtonHeat")

        # first slot
        m.addConstrs(CumHeat[k,1] == gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,1,c] for c in coils) for k in casters)
        # cumulative (unless change of heat)
        m.addConstrs(CumHeat[k,s] >= CumHeat[k,s-1] + gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) - MAX_HEAT_WEIGHT*HeatChg[k,s] for k in casters for s in slots if s>1)
        m.addConstrs(CumHeat[k,s] <= CumHeat[k,s-1] + gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils)  for k in casters for s in slots if s>1)
        # restart (if change of heat)
        m.addConstrs(CumHeat[k,s] <=  gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) + MAX_HEAT_WEIGHT*(1-HeatChg[k,s]) for k in casters for s in slots)
        # allow intermixing
        m.addConstrs(CumHeat[k,s] + (MAX_HEAT_WEIGHT-CumHeat[k,s-1] )>=  gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) for k in casters for s in slots if s>1)
        # force minimum tonn for heat change
        m.addConstrs( MIN_HEAT_WEIGHT * HeatChg[k,s] <= CumHeat[k,s-1] for k in casters for s in slots if s > 1)

        # if grande change, force heat change and forbid intermixing 
        if useGradeCuts < 2:
            m.addConstrs(  GradeChange[k,s-1] <= HeatChg[k,s] for k in casters for s in slots if s > 1)
            m.addConstrs(CumHeat[k,s] >=  gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) - MAX_HEAT_WEIGHT*(1-GradeChange[k,s-1]) for k in casters for s in slots if s>1)
        else:
            m.addConstrs(  gp.quicksum(GChg[k,s-1,g1,g2] for (g1,g2) in gradeChanges.keys()) <= HeatChg[k,s] for k in casters for s in slots if s > 1)
            m.addConstrs(CumHeat[k,s] >=  gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) - MAX_HEAT_WEIGHT*(1-gp.quicksum(GChg[k,s-1,g1,g2] for (g1,g2) in gradeChanges.keys())) for k in casters for s in slots if s>1)
    else:
        MAX_NUM_HEAT = np.ceil(np.sum(np.sort([sol.coils[c].weight*(1/2000) for c in sol.coils.keys()])[sol.numSlots:])/MAX_HEAT_WEIGHT)
        MAX_CUM_HEAT_WEIGHT = MAX_NUM_HEAT*MAX_HEAT_WEIGHT
        NumHeats = m.addVars(casters, slots, vtype=GRB.INTEGER, lb=1, ub=MAX_NUM_HEAT, name="numHeatsInSeq")
        CumHeat = m.addVars(casters, slots, lb=0, ub=MAX_CUM_HEAT_WEIGHT, name="cumtonHeat")

        # first slot
        m.addConstrs(CumHeat[k,1] == gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,1,c] for c in coils) for k in casters)
        # cumulative (unless change of heat)
        m.addConstrs(CumHeat[k,s] <= CumHeat[k,s-1] + gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils)  for k in casters for s in slots if s>1)
        m.addConstrs(CumHeat[k,s] >= CumHeat[k,s-1] + gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) - MAX_CUM_HEAT_WEIGHT*GradeChange[k,s-1] for k in casters for s in slots if s>1)
        # restart cumHeat if grade change
        m.addConstrs(CumHeat[k,s] >= gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) for k in casters for s in slots if s>1)
        
        m.addConstrs(CumHeat[k,s] <=  gp.quicksum(sol.coils[c].weight*(1/2000)* X[k,s,c] for c in coils) + MAX_CUM_HEAT_WEIGHT*(1-GradeChange[k,s-1]) for k in casters for s in slots if s>1)

        # force min and max tonn at grade change
        # Idea: gap of +1 (enough for any combination) except at grande changes
        m.addConstrs( MAX_HEAT_WEIGHT * NumHeats[k,s] >= CumHeat[k,s] for k in casters for s in slots)

        m.addConstrs( CumHeat[k,s] >= MIN_HEAT_WEIGHT * (NumHeats[k,s] - 1 + GradeChange[k,s]) for k in casters for s in slots)

# Bounds
    m.addConstr(gp.quicksum(GradeChange[k,s] for k in casters for s in slots) >= len(grades) -2)
    # sort by weight
    m.addConstrs(gp.quicksum(HeatChg[k,s] for s in slots) >= MIN_NUM_HEAT_CHANGES for k in casters)
    m.addConstrs(gp.quicksum(HeatChg[k,s] for s in slots) <= MAX_NUM_HEAT_CHANGES for k in casters)


#10) Improving grade change cost

    # Extra vars detailing grade change
    if useGradeCuts >= 2:
        # improving LB on number of grade changes
        for g in grades:
                
            m.addConstr( gp.quicksum(GChg[k,s,g,g2] for k in casters for s in slots for g2 in sol.getGrades() if (g,g2) in gradeChanges.keys()) + gp.quicksum( X[k,sol.numSlots,c] for k in casters for c in coils if sol.coils[c].grade == g ) >= 1 )
                
            m.addConstr( gp.quicksum(GChg[k,s,g2,g] for k in casters for s in slots for g2 in sol.getGrades() if (g2,g) in gradeChanges.keys()) + gp.quicksum( X[k,1,c] for k in casters for c in coils if sol.coils[c].grade == g ) >= 1 )

    if useGradeCuts == 1 or useGradeCuts >= 3:
        
        # adding extra variables
        Z = m.addVars(grades, casters, vtype=GRB.BINARY, name="GradeAppearsInK")
        bigM_GC = {g: min(sol.numSlots,sum(sol.coils[c].grade == g for c in sol.coils)) for g in grades} 
        m.addConstrs( Z[g,k]*bigM_GC[g] >= gp.quicksum(X[k,s,c] for s in slots for c in coils if sol.coils[c].grade == g) for k in casters for g in grades)
        m.addConstrs( gp.quicksum(Z[g,k] for k in casters) >= 1 for g in grades)

        # if useGradeCuts>=3 then I can improve lower bound with Z variables
        if useGradeCuts >= 3:
            for g in grades:
                m.addConstrs( gp.quicksum(GChg[k,s,g,g2] for s in slots for g2 in sol.getGrades() if (g,g2) in gradeChanges.keys()) + gp.quicksum( X[k,sol.numSlots,c] for c in coils if sol.coils[c].grade == g ) >= Z[g,k] for k in casters)
                
                m.addConstrs( gp.quicksum(GChg[k,s,g2,g] for s in slots for g2 in sol.getGrades() if (g2,g) in gradeChanges.keys()) + gp.quicksum( X[k,1,c] for c in coils if sol.coils[c].grade == g ) >= Z[g,k] for k in casters)

        # set of grade combinations
        alpha = computeAlpha(sol)
        import itertools
        for L in range(2, len(grades)+1):
            for subset in itertools.combinations(grades, L):
                for phi in itertools.permutations(subset, len(subset)):
                    #print("φ subset:", phi, end=' ')
                    notin = []
                    for i in grades:
                        if i not in phi:
                            notin.append(i)
                #print("Not in φ:",notin)
                
                #Cuts (#46-#47 in pdf)
                    # CASE 1: only GradeChange
                    if useGradeCuts == 1:
                        m.addConstrs( gp.quicksum(GradeChangeCost[k,s] for s in slots) >= alpha[phi] * (1 - gp.quicksum((1 - Z[g1,k]) for g1 in phi) - gp.quicksum( Z[g2,k] for g2 in notin)) for k in casters)
                    if useGradeCuts == 4:
                        for g in phi:
                        # CASE >1 : with grade change detail vars
                            m.addConstrs( gp.quicksum(GChg[k,s,g2,g] for s in slots for g2 in phi if g2 != g) +
                                        gp.quicksum(X[k,1,c] for c in coils if sol.coils[c].grade == g) >=
                                        (1 - gp.quicksum((1 - Z[g1,k]) for g1 in phi) - gp.quicksum(Z[g2,k] for g2 in notin))
                                        for k in casters)
                            
                            m.addConstrs( gp.quicksum(GChg[k,s,g,g2] for s in slots for g2 in phi if g2 != g) +
                                        gp.quicksum(X[k,sol.numSlots,c] for c in coils if sol.coils[c].grade == g) >=
                                        (1 - gp.quicksum((1 - Z[g1,k]) for g1 in phi) - gp.quicksum(Z[g2,k] for g2 in notin))
                                        for k in casters)
                        

# Critical Pairs Cuts
    if addCriticalPairs:
        nCritic = len(sol.criticPairs)
        Zcritic = m.addVars(range(nCritic), casters, slots, name="criticChange")
        for i in range(nCritic):
            (idxSmall,idxLarge) = sol.criticPairs[i]
            m.addConstrs( gp.quicksum(X[k,s,c] for c in sol.sortedDecWidth[0:idxSmall+1]) <= Zcritic[i,k,s] for k in casters for s in slots )
            m.addConstrs( gp.quicksum(X[k,s,c] for c in sol.sortedIncWidthUB[0:idxLarge+1]) <= 1-Zcritic[i,k,s] for k in casters for s in slots )
            m.addConstrs( Zcritic[i,k,s+1] <= Zcritic[i,k,s] for k in casters for s in slots if s < len(slots)-1 )

            if addCriticalPairs == 2:
                m.addConstrs( gp.quicksum(X[k,ss,sol.sortedDecWidth[idxSmall]] for ss in slots if ss >= s) <= Zcritic[i,k,s] for k in casters for s in slots)
                m.addConstrs( gp.quicksum(X[k,ss,sol.sortedIncWidthUB[idxLarge]] for ss in slots if ss <= s) <= 1-Zcritic[i,k,s] for k in casters for s in slots)
            elif addCriticalPairs == 3:
                m.addConstrs( gp.quicksum(X[k,ss,c] for ss in slots if ss >= s) <= Zcritic[i,k,s] for k in casters for s in slots for c in sol.sortedDecWidth[0:idxSmall+1])
                m.addConstrs( gp.quicksum(X[k,ss,c] for ss in slots if ss <= s) <= 1-Zcritic[i,k,s] for k in casters for s in slots for c in sol.sortedIncWidthUB[0:idxLarge+1])

# Objective Function
    if  useGradeCuts >= 2:
        GradeChangeObj = sum( gradeChanges[g1,g2]*GChg[k,s,g1,g2] for k in casters for s in slots for (g1,g2) in gradeChanges.keys())
    else:
        GradeChangeObj = sum( GradeChangeCost[k,s] for k in casters for s in slots)
    GaugeJumpObj = sum( GaugeJump[k,s] for k in casters for s in slots)
    IntCasterDiffObj = sum(sol.penaltyWidthDifference*P[s] for s in slots)
    b_bar = sum(sol.coils[c].length for c in sol.coils.keys())/len(sol.coils) #av length
    TrimLossObj = sum(sol.penaltyTrimLoss*( V[k,s]-sum(sol.coils[c].orderWidth * X[k,s,c] for c in coils) )*3.5433*b_bar*0.2817929*(1/2000) for k in casters for s in slots)
    RollChangesObj = sum( sol.penaltyRollChange*Ybar[s] for s in slots)
    
    m.setObjective(GradeChangeObj+GaugeJumpObj+IntCasterDiffObj+TrimLossObj+RollChangesObj)

# finish formulation
    m.update()

    #Storing the model
    sol.m = m

    #Storing variables
    m._varX = X
    m._varV = V
    m._varYbar = Ybar
    if useGradeCuts < 2:
        m._varGC = GradeChange
    else:
        m._varGC = GChg
    m._varHeatChg = HeatChg


def useStartingSolution(sol):
    # fix initial value from Schoppy Output
    fixInitial = sol.m.addConstrs(sol.m._varX[k,s,sol.schedule[(k,s)]] == 1 for (k,s) in list(sol.schedule))
    optimize(sol, int(1e10))
    if sol.m.status == 3:
        print("ERROR: No solution. Infeasible.  Computting IIS")
        sol.m.computeIIS()
        sol.m.write('infeasible.ilp')
    else:
        print("Initial solution with value "+str(sol.m.ObjBound))
        for v in sol.m.getVars():
            v.start = v.x
    sol.m.remove(fixInitial)
    sol.m.update()


def optimize(sol, timeLimit=86400):
    sol.m.Params.timeLimit = timeLimit
    sol.m.Params.MIPFocus = 2 # from gurobi params tune
    sol.m.Params.Symmetry = 2 # from gurobi params tune
    if sol.m._useLazyCuts:
        sol.m.optimize(Callback)
    else:
        sol.m.optimize()
    if sol.m.status == 3:
        print("ERROR: No solution. Infeasible")
    else:
        solV = sol.m.getAttr('x',sol.m._varV)
        solX = sol.m.getAttr('x',sol.m._varX)
        solY = sol.m.getAttr('x',sol.m._varYbar)
        #solH = sol.m.getAttr('x',sol.m._varHbar)

        #Update roller use
        sol.rollerChanges = []
        for s in range(2,sol.numSlots+1):
            if solY[s] > 0:
                sol.rollerChanges.append(s)

        #Replace width
        sol.castWidth = {}
        for (k,s) in solV.keys():
            sol.castWidth[(k,s)]= solV[k,s]

        sol.schedule={}
        #Update schedule
        for k in sol.m._casters:
            for s in sol.m._slots:
                for c in sol.m._coils:
                    if solX[k,s,c]>0.999:
                        sol.schedule[(k,s)] = c

def reoptimizeSolution(sol):
    formulateModel(sol)
    fixInitial = sol.m.addConstrs(sol.m._varX[k,s,sol.schedule[(k,s)]] == 1 for (k,s) in list(sol.schedule))
    optimize(sol, int(1e10))


