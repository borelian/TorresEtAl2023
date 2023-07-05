import csv
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import matplotlib.lines as mlines
import numpy as np
import math

from parameters import *


# Set of valid coil edge types
COIL_EDGES = {'M', 'C', 'HRB'}

# Set of valid coil grade types
COIL_GRADES = {'Grade_0', 'Grade_1', 'Grade_2', 'Grade_3', 'Grade_4', 'Grade_5',
       'Grade_6', 'Grade_7', 'Grade_8', 'Grade_9', 'Grade_10', 'Grade_11',
       'Grade_12', 'Grade_13', 'Grade_14', 'Grade_15', 'Grade_16', 'Grade_17',
       'Grade_18', 'Grade_19', 'Grade_20', 'Grade_21', 'Grade_22', 'Grade_23',
       'Grade_24', 'Grade_25', 'Grade_26', 'Grade_27', 'Grade_28', 'Grade_29',
       'Grade_30', 'Grade_31', 'Grade_32', 'Grade_33', 'Grade_34', 'Grade_35',
       'Grade_36', 'Grade_37', 'Grade_38', 'Grade_39', 'Grade_40', 'Grade_41',
       'Grade_42', 'Grade_43', 'Grade_44', 'Grade_45', 'Grade_46', 'Grade_47',
       'Grade_48', 'Grade_49', 'Grade_50', 'Grade_51', 'Grade_52', 'Grade_53',
       'Grade_54', 'Grade_55', 'Grade_56', 'Grade_57', 'Grade_58', 'Grade_59',
       'Grade_60'}

# Grade change costs from ($ per ton)
# tuple <grade_from, grade_to> to value
filename = 'grades_hash.csv'
GRADE_COST = {}
with open(filename) as fileData:
    reader = csv.reader(fileData, delimiter='\t')
    header = reader.__next__()
    for line in reader:
        for i in range(1,len(header)):
            if line[i] != ".":
                GRADE_COST[(header[i],line[0])] = float(line[i])

#%%
def getGradeCost(grade_from, grade_to):
    cost = GRADE_COST.get((grade_from,grade_to))
    if grade_from == grade_to:
        return 0
    elif cost == None:
        return COST_FOR_FORBIDDEN_GRADES
    else:
        return cost * MIN_HEAT_SIZE + MIN_GRADE_CHANGE_PENALTY

def GetInitialSolution(filename):
    """ read Schoppy output to construct an array with its solution"""
    # header: 
    # [0]Coil [1]Grade [2]Gauge [3]OrderWidth [4]Weight [5]CoilLength 
    # [6]EdgeCode [7]CastWidth [8]Caster [9]RollerCampaign [10]Slot
    Schoppy_Assignment = []
    with open(filename) as fileData:
        reader = csv.reader(fileData, delimiter='\t')
        for line in reader:
            if reader.line_num > 22:
                #print(line[8], line[10], line[0])
                Schoppy_Assignment.append([line[8], line[10], line[0]])

    return Schoppy_Assignment

def getGaugeAllowed(coilGauge):
    if coilGauge >= HEAVY_GAUGE:
        return 0.5 * coilGauge
    elif coilGauge >= MEDIUM_GAUGE:
        return 0.75 * coilGauge
    else:
        return 0.9 * coilGauge
        

def computeAlpha(sol):
    
        import itertools

        tipos = list(sol.getGrades())

        alpha = {}

        for L in range(2, len(tipos)+1):
            for subset in itertools.combinations(tipos, L):
                #print("Combination:")
                #print(subset)
                min_cost = 9999
                min_comb = 0
                for subset2 in itertools.permutations(subset, len(subset)):
                    #print("  Posibilities:")
                    #print("  ",subset2)
                    cost = getGradeCost(subset2[0],subset2[1])
                    #print("  Cost:",cost)
                    if cost<min_cost:
                        min_cost = cost
                        min_comb = subset2
                #print("    The best possible combination is",min_comb,'with cost:',min_cost)

                for subset2 in itertools.permutations(subset, len(subset)):

                    alpha[subset2] = min_cost
        print(alpha)
        return alpha



#%%
class Coil:
    """a class for coils that need to be casted"""

    def __init__(self, coilId, coilEdge, coilGrade, coilOrderWidth, coilGauge, coilWeight, coilLength, coilDegradation = None):
        self.id = coilId
        self.orderWidth = coilOrderWidth
        self.gauge = coilGauge
        self.weight = coilWeight
        self.length = coilLength
        if coilEdge not in COIL_EDGES:
            raise ValueError("invalid Edge type")
        self.edge = coilEdge
        if coilGrade not in COIL_GRADES:
            raise ValueError("invalid Grade type")
        self.grade = coilGrade
        if coilDegradation == None:
            self.degradation =  max(1,-716.68*self.gauge**3 + 424.91*self.gauge**2 - 84.595*self.gauge + 6.7878)
        else:
            self.degradation = coilDegradation
     
    def getGaugeMin(self):
        coilGauge = self.gauge
        min_gauge_next = 0
        
        if coilGauge >= HEAVY_GAUGE:
            min_gauge_next = 0.5 * coilGauge
        elif coilGauge >= MEDIUM_GAUGE:
            min_gauge_next = 0.75 * coilGauge
        else:
            min_gauge_next = 0.9 * coilGauge
        
        return(min_gauge_next)

    def getWidthUB(self):
        coilOrderWidth = self.orderWidth
        coilEdge = self.edge
        wUB = 0
        
        if coilEdge == 'M':
            wUB = coilOrderWidth + MAX_TRIMLOSS_MILL_ALTERNATIVE
        if coilEdge == 'HRB':
            wUB = coilOrderWidth + MAX_TRIMLOSS_HRB
        if coilEdge == 'C':
            wUB = coilOrderWidth + MAX_TRIMLOSS_CUT
        
        return(wUB)

 

            

#%%
class SteelccProblem:
    """class with the Steelcc problem to solve"""

    def __init__(self):
        self.coils={}
        # Penalty for width differences between casters [$/in]
        self.penaltyWidthDifference = PENALTY_WIDTH_DIFFERENCE
        # Penalty for width differences trim loss [$/in]
        self.penaltyTrimLoss = PENALTY_TRIM_LOSS
        # Penalty for roll change [$]
        self.penaltyRollChange = PENALTY_ROLL_CHANGE
        # Penalty for excessive gauge decrease
        self.penaltyGaugeDecrease = PENALTY_GAUGE_DECREASE

    def addCoil(self, coil):
        """add coil to the instance, indexed by coilId"""
        self.coils[coil.id]= coil
   
    def getGrades(self):
        """get set of grades that appears in the instance"""
        grades =  set()
        for c in self.coils.keys():
            grades.add(self.coils[c].grade)
        return grades

    def getNumCoils(self):
        """get number of coils"""
        return len(self.coils)
    
    def getNumGrades(self):
        """get number of different coil grades"""
        grades = self.getGrades()
        return len(grades)
    
    def getMinOrderWidth(self):
        """return minimum among all order widths"""

    def getSubproblem(self, numCoils):
        subset = np.random.choice(len(self.coils), numCoils, replace=False)
        listCoils = list(self.coils.keys())
        newProb = SteelccProblem()
        for i in subset:
            newProb.coils[listCoils[i]] = self.coils[listCoils[i]]
        return newProb

    def getHalfproblem(self):
        newProb = SteelccProblem()
        newlength = len(self.coils.keys())/2
        if newlength%2 == 1:
            newlength -= 1
        for i in range(1,int(newlength)+1):
            newProb.coils[i] = self.coils[i]
        return newProb
    
    def computeCriticPairs(self):
        self.sortedDecWidth = sorted(self.coils, key = lambda x : self.coils[x].orderWidth, reverse=True)
        self.sortedIncWidthUB = sorted(self.coils, key = lambda x : self.coils[x].getWidthUB())
        self.criticPairs = []
        i = 0
        while i < len(self.coils):
            j=-1
            while self.coils[self.sortedDecWidth[i]].orderWidth > self.coils[self.sortedIncWidthUB[j+1]].getWidthUB():
                j=j+1
            while self.coils[self.sortedDecWidth[i+1]].orderWidth > self.coils[self.sortedIncWidthUB[j]].getWidthUB():
                i=i+1
            if j > -1:
                self.criticPairs.append((i,j))
                i = i+1
            else:
                break


#%%
class SteelccSolution(SteelccProblem):
    """class that provide solution to an instance"""

    def __init__(self):
        SteelccProblem.__init__(self)
        # default
        self.numCasters = 0
        self.numSlots = 0
        # assignment of coils to caster-slots
        # a dictionary of tuples <caster, slot> to coilId assigned
        self.schedule = {}
        # a dictionary of tuples <caster, slot> to cast width
        self.castWidth = {}
        # roller changes (list of slots when roller is changed)
        self.rollerChanges = []
        # to store the model
        self.m = 0

    def importProblem(self, problem):
        self.coils= problem.coils
        self.penaltyWidthDifference = problem.penaltyWidthDifference
        self.penaltyTrimLoss = problem.penaltyTrimLoss
        self.penaltyRollChange = problem.penaltyRollChange
        self.penaltyGaugeDecrease = problem.penaltyGaugeDecrease
    
    def computesForCuts(self):
        
        coils = list(self.coils)
        #Create a dictionary with coil ids and order widths, in descending order
        orderWidth_aux = {self.coils[i].id : self.coils[i].orderWidth for i in coils}
        DoW = {k: v for k, v in sorted(orderWidth_aux.items(), key=lambda item: item[1], reverse=True)}
        
        #Create a dictionary with coil ids and the upper bound on order widths, in ascending order
        UBorderWidth_aux = {self.coils[i].id : self.coils[i].getWidthUB() for i in coils}
        AoW_ub = {k: v for k, v in sorted(UBorderWidth_aux.items(), key=lambda item: item[1], reverse=False)}

        #Create lists to index ids for coils that have to be first because its width is greater than others UB order width
        firsts_id = []
        lasts_id = []
        for i in range(0,self.numSlots*2):
            if (list(DoW.values())[i]) > (list(AoW_ub.values())[i]):
                firsts_id.append((list(DoW.keys())[i]))
                lasts_id.append((list(AoW_ub.keys())[i]))
            else:
                break
         
        return firsts_id, lasts_id


    def getCosts(self):
        """evaluate the cost of the solution"""
        # dictionary of each type of cost
        cost = {}
        # evalute trim cost
        trimCost = 0
        gaugeCost = 0
        gradeCost = 0
        casterDiffWidthCost = 0
        for (k,s) in list(self.schedule):
            # get coil assigned to this position
            coil = self.coils[self.schedule[(k,s)]]
            # compute trim difference
            trimCost += (self.castWidth[(k,s)] - coil.orderWidth) * 3.5433 * coil.length * 0.2817929 * (1/2000)
            # compute gauge and grade jump
            if s < self.numSlots:
                nextCoil = self.coils[self.schedule[(k,s+1)]]
                # compute gauge cost
                if coil.gauge >= HEAVY_GAUGE:
                    gaugeCost += max(0,coil.gauge * 0.5 - nextCoil.gauge)
                elif coil.gauge >= MEDIUM_GAUGE:
                    gaugeCost += max(0,coil.gauge * 0.75 - nextCoil.gauge)
                else:
                    gaugeCost += max(0,coil.gauge * 0.90 - nextCoil.gauge)
                #conpute grade cost
                gradeCost += getGradeCost(coil.grade, nextCoil.grade)
        # evaluate caster width difference
        if self.numCasters == 2:
            for s in range(1,self.numSlots+1):
                diffWidth = abs(self.castWidth[1,s]-self.castWidth[2,s])
                casterDiffWidthCost += max(0,diffWidth - ALLOWABLE_INTERCASTER_WIDTH)
        else:
            casterDiffWidthCost = 0

        ## construct dictionary with al costs       
        cost['grade_change_cost'] = gradeCost
        cost['gauge_change_cost'] = self.penaltyGaugeDecrease * gaugeCost
        cost['width_difference_cost'] = self.penaltyWidthDifference * casterDiffWidthCost
        cost['trim_loss_cost'] = self.penaltyTrimLoss * trimCost
        cost['roller_usage_cost'] = COST_ROLLER_CHANGE*(len(self.rollerChanges)+1)
        cost['objective_function_value'] = cost['grade_change_cost'] + cost['gauge_change_cost'] + cost['width_difference_cost'] + cost['trim_loss_cost'] + cost['roller_usage_cost'] 
        return cost


    def assign(self, coil, caster, slot, width):
        """asign a coil to a caster-slot with a given casted width"""
        self.schedule[(caster,slot)] = coil
        self.castWidth[(caster,slot)] = width

    def changeCastWidth(self, caster, slot, width):
        self.castWidth[(caster, slot)] = width

                
    def loadSchoppyOutput(self,filename, maxCoils = 1000, numCaster=2):
        rollers = {}
        with open(filename) as fileData:
            reader = csv.reader(fileData, delimiter=',')
            slot=0
            flag=0
            for line in reader:
                if reader.line_num > 1:
                    if flag==0:
                        if int(line[8])==2:
                            slot=0
                            flag=1
                    slot = slot+1
                    # read coil information
                    coil = Coil(reader.line_num-1,line[5].split()[0],line[0],float(line[2]),float(line[1]),float(line[3]),float(line[4]))
                    self.addCoil(coil)
                    #coilId, coilEdge, coilGrade, coilOrderWidth, coilGauge, coilLength, coilDegradation = None
                    #print("Coil",reader.line_num-1,"EDGE",line[5],"GRADE:",line[0],"OrderW",float(line[2]),"gauge",float(line[1]),"Weight",float(line[3]),"Length",float(line[4]))
                    # read assignment to caster-slot and its width
                    self.assign(reader.line_num-1,int(line[8]),slot,float(line[7]))
                    # read roller in each slot
                    rollers[slot] = int(line[9])
                    if self.getNumCoils() == maxCoils:
                        break
        self.numCasters =numCaster 
        self.numSlots = int(self.getNumCoils()/numCaster)
        for s in range(2,self.numSlots+1):
            if rollers[s] > rollers[s-1]:
                self.rollerChanges.append(s)

    def getGaugesFromCaster(self, caster):
        """ returns a list with gauges of coils in caster, sorted by slot"""
        out = []
        for s in range(self.numSlots):
            coil = self.schedule[(caster,s+1)]
            out.append(self.coils[coil].gauge)
        return out

    def getOrderWidthFromCaster(self, caster):
        """ returns a list with orderWidths of coils in caster, sorted by slot"""
        out = []
        for s in range(self.numSlots):
            coil = self.schedule[(caster,s+1)]
            out.append(self.coils[coil].orderWidth)
        return out
    
    def getCastWidthFromCaster(self,caster):
        return [self.castWidth[(caster,s)] for s in range(1,self.numSlots+1)]

    def getGradeFromCaster(self, caster):
        """ returns a list with orderWidths of coils in caster, sorted by slot"""
        out = []
        for s in range(self.numSlots):
            coil = self.schedule[(caster,s+1)]
            out.append(self.coils[coil].grade)
        return out   

    def validateCoilWidth(self):
        returnValue = True

        for (k,s) in list(self.schedule):
            # get coil assigned to this position
            coil = self.coils[self.schedule[(k,s)]]
            # get cated width for the coil
            castedWidth = self.castWidth[(k,s)]

            # check minimum width
            if castedWidth < coil.orderWidth - CHECK_TOLERANCE:
                print("ERROR: coil ", coil.id, " of edge ", coil.edge,
                 " and order width ", str(coil.orderWidth),
                 " has been casted with width ", str(castedWidth),
                 " which is too narrow.")
                returnValue = False
            
            # check maximum width for HRB
            if coil.edge == "HRB":
                if castedWidth > coil.orderWidth + MAX_TRIMLOSS_HRB + CHECK_TOLERANCE:
                    print("ERROR: coil ", coil.id, " of edge ", coil.edge,
                     " has been casted with width ", str(castedWidth),
                     " which is outside the limits ",
                     "[", str(coil.orderWidth), ",", str(coil.orderWidth + MAX_TRIMLOSS_HRB), "].")
                    returnValue = False
            # check maximum width for HRB
            elif coil.edge == "C":
                if castedWidth > coil.orderWidth + MAX_TRIMLOSS_CUT + CHECK_TOLERANCE:
                    print("ERROR: coil ", coil.id, " of edge ", coil.edge,
                     " has been casted with width ", str(castedWidth),
                     " which is outside the limits ",
                     "[", str(coil.orderWidth), ",", str(coil.orderWidth + MAX_TRIMLOSS_CUT), "].")
                    returnValue = False
            # check maximum width for MILL
            else:
                if castedWidth > coil.orderWidth + MAX_TRIMLOSS_MILL + CHECK_TOLERANCE:
                    # potentially, it can be casted with alternative width limits
                    if castedWidth < coil.orderWidth + MIN_TRIMLOSS_MILL_ALTERNATIVE:
                        print("ERROR: coil ", coil.id, " of edge ", coil.edge,
                         " has been casted with width ", str(castedWidth),
                         " which is outside the limits ",
                         "[", str(coil.orderWidth), ",", str(coil.orderWidth + MAX_TRIMLOSS_MILL), "] or ",
                         "[", str(coil.orderWidth+MIN_TRIMLOSS_MILL_ALTERNATIVE), ",", str(coil.orderWidth + MAX_TRIMLOSS_MILL_ALTERNATIVE), "].")
                        returnValue = False
                    elif castedWidth > coil.orderWidth + MAX_TRIMLOSS_MILL_ALTERNATIVE + CHECK_TOLERANCE:
                        print("ERROR: coil ", coil.id, " of edge ", coil.edge,
                         " has been casted with width ", str(castedWidth),
                         " which is outside the limits ",
                         "[", str(coil.orderWidth), ",", str(coil.orderWidth + MAX_TRIMLOSS_MILL), "] or ",
                         "[", str(coil.orderWidth+MIN_TRIMLOSS_MILL_ALTERNATIVE), ",", str(coil.orderWidth + MAX_TRIMLOSS_MILL_ALTERNATIVE), "].")
                        returnValue = False
        return returnValue
        
    def validateCasterWidth(self):
        returnValue = True

        for s in range(1,self.numSlots):
            # check consecutive difference
            for k in range(1,self.numCasters+1):
                if self.castWidth[(k,s)] < self.castWidth[(k,s+1)] - CHECK_TOLERANCE:
                    print("ERROR: caster ", k, "slot ", s, " has width ", self.castWidth[(k,s)],
                        " but next slot has width ",  self.castWidth[(k,s+1)],
                        " which is larger")
                    returnValue = False
                if self.castWidth[(k,s)] > self.castWidth[(k,s+1)] + MAX_CASTER_DECREASE + CHECK_TOLERANCE :
                    print("ERROR: caster ", k, "slot ", s, " has width ", self.castWidth[(k,s)],
                        " but next slot has width ",  self.castWidth[(k,s+1)],
                        " which is too narrow")
                    returnValue = False

        # check intercaster difference
        if (self.numCasters == 2):
            for s in range(1,self.numSlots+1):
                if abs(self.castWidth[(1,s)] - self.castWidth[(2,s)]) > ALLOWABLE_INTERCASTER_WIDTH + MAX_EXTRA_INTERCASTER_WIDTH + CHECK_TOLERANCE:
                    print("ERROR: caster width difference at slot ", s, " is too wide",
                    " caster 1 is ", str(self.castWidth[(1,s)]), 
                    " and caster 2 is ", str(self.castWidth[(2,s)]))
                    returnValue = False
        
        # return
        return returnValue

    def validateRollerChange(self):
        returnValue = True

        # check valid gauge for new roller
        for s in [1]+self.rollerChanges:
            for k in range(1,self.numCasters+1):
                coil = self.coils[self.schedule[(k,s)]]
                if coil.gauge < MIN_GAUGE_ROLLER_CHANGE - CHECK_TOLERANCE:
                    print("ERROR: new roller in slot", s, " but coil ", coil.id, " in caster", k, " as gauge ", coil.gauge,
                     "which is too narrow to start a new roller campaign")
                    returnValue = False
        
        # check degradation
        start = 1
        for r in self.rollerChanges+[self.numSlots+1]:
            degradation=0
            for s in range(start, r):
                degradation += self.coils[self.schedule[(1,s)]].degradation
                if (self.numCasters == 2):
                    degradation += self.coils[self.schedule[(2,s)]].degradation

            if degradation > MAX_ROLLER_DEGRADATION + CHECK_TOLERANCE:
                print("ERROR: cumulated degradation for roller between slots", start, "and", r-1, "is", degradation)
                returnValue = False
            start = r
        return returnValue



    def validateHeat(self):

        maxHeatSizeCheck = math.ceil(MAX_HEAT_WEIGHT/(MAX_HEAT_WEIGHT-MIN_HEAT_WEIGHT))

        for k in range(1,self.numCasters+1):
            heatOk = False
            heatTon = self.coils[self.schedule[(k,1)]].weight/2000
            for s in range(2,self.numSlots):
                if self.coils[self.schedule[(k,s)]].grade != self.coils[self.schedule[(k,s-1)]].grade :
                    ## Grade change
                    for idx in range(1,maxHeatSizeCheck):
                        if (heatTon > MIN_HEAT_WEIGHT*idx) and (heatTon < MAX_HEAT_WEIGHT*idx):
                            heatOk = True
                            # print("Current heat of size ", heatTon, " fit in ", idx)
                    if heatTon >  MIN_HEAT_WEIGHT*maxHeatSizeCheck:
                            heatOk = True
                            # print("Current heat of size ", heatTon, " large enough")
                    if heatOk == False:
                        print("ERROR: current heat has incompatible tonnage ", heatTon)
                        return False
                    
                    heatOk = False
                    heatTon = 0
                heatTon += self.coils[self.schedule[(k,s)]].weight/2000
        return True











    def validateSolution(self):
        out = True
        out = self.validateCoilWidth() & out
        out = self.validateCasterWidth() & out
        out = self.validateRollerChange() & out
        out = self.validateHeat() & out
        return out


    def validateGradeChange(self):
        returnValue = True

        for k in {1,2}:
            for s in range(1,self.numSlots):
                # get consecutive coils
                coil = self.coils[self.schedule[(k,s)]]
                nextcoil = self.coils[self.schedule[(k,s+1)]]
                if getGradeCost(coil.grade, nextcoil.grade) == COST_FOR_FORBIDDEN_GRADES:
                    print("WARNING: coil ", nextcoil.id, " of grade ", nextcoil.grade,
                     " is casted after coil ", coil.id, "of grade ", coil.grade, 
                     " at caster", k, "slot", s+1,
                     " which is an invalid grade transition")
                    returnValue = False 
        return returnValue

    def validateExcesiveGaugeDecrease(self):
        returnValue = True

        for k in {1,2}:
            for s in range(1,self.numSlots):
                # get consecutive coils
                coil = self.coils[self.schedule[(k,s)]]
                nextcoil = self.coils[self.schedule[(k,s+1)]]
                if coil.gauge >= HEAVY_GAUGE:
                    if nextcoil.gauge < coil.gauge * 0.15 - CHECK_TOLERANCE:
                        print("WARNING: coil ", nextcoil.id, " of gauge ", nextcoil.gauge,
                         " is casted after coil ", coil.id, "of gauge ", coil.gauge,
                         " at caster", k, "slot", s+1,
                         " which is too thin for a HEAVY gauge ",
                         " (","%.2f" % (100*(1-(nextcoil.gauge/coil.gauge))),"% decrease)")
                        returnValue = False 
                elif coil.gauge >= MEDIUM_GAUGE:
                    if nextcoil.gauge < coil.gauge * 0.70 - CHECK_TOLERANCE:
                        print("WARNING: coil ", nextcoil.id, " of gauge ", nextcoil.gauge,
                         " is casted after coil ", coil.id, "of gauge ", coil.gauge, 
                         " at caster", k, "slot", s+1,
                         " which is too thin for a MEDIUM gauge ",
                         " (","%.2f" % (100*(1-(nextcoil.gauge/coil.gauge))),"% decrease)")
                        returnValue = False 
                else:
                    if nextcoil.gauge < coil.gauge * 0.90 - CHECK_TOLERANCE:
                        print("WARNING: coil ", nextcoil.id, " of gauge ", nextcoil.gauge,
                         " is casted after coil ", coil.id, "of gauge ", coil.gauge, 
                         " at caster", k, "slot", s+1,
                         " which is too thin for a LIGHT gauge ",
                         " (","%.2f" % (100*(1-(nextcoil.gauge/coil.gauge))),"% decrease)")
                        returnValue = False 
        return returnValue


    def validateSolutionStrict(self):
        return self.validateSolution() & self.validateGradeChange() #& self.validateExcesiveGaugeDecrease()




    def plotSolution(self, omitGrades=False, saveFile=False, saveFileName='figure.pdf'):
        maxPriY = 1.15*self.coils[max(self.coils, key=lambda x: self.coils[x].orderWidth)].orderWidth
        minPriY = 0.6*self.coils[min(self.coils, key=lambda x: self.coils[x].orderWidth)].orderWidth
        maxAltY = 2.5*self.coils[max(self.coils, key=lambda x: self.coils[x].gauge)].gauge
        minAltY = 0

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        # plot order width
        ax1.plot(range(1,self.numSlots+1),self.getOrderWidthFromCaster(1), 'k+', markersize=4)
        if (self.numCasters == 2):
            ax1.plot(range(1,self.numSlots+1),self.getOrderWidthFromCaster(2), 'ko',  markersize=4, fillstyle='none')
        # plot cast width
        ax1.plot(range(1,self.numSlots+1),self.getCastWidthFromCaster(1), 'k-', linewidth=2)
        if (self.numCasters == 2):
            ax1.plot(range(1,self.numSlots+1),self.getCastWidthFromCaster(2), 'k--', linewidth=2)
        # plot gauges
        #         if coilGauge >= HEAVY_GAUGE:
        #    min_gauge_next = 0.5 * coilGauge
        #elif coilGauge >= MEDIUM_GAUGE:
        #    min_gauge_next = 0.75 * coilGauge
        #else:
        #    min_gauge_next = 0.9 * coilGauge
        lstyle = ['-','--'] 
        for k in range(1,self.numCasters+1):
            gauge = self.getGaugesFromCaster(k)
            for s in range(2,self.numSlots+1):
                if gauge[s-1] < getGaugeAllowed(gauge[s-2]):
                    ax2.plot([s-1,s],[gauge[s-2],gauge[s-1]],'r', linestyle=lstyle[k-1])
                else:
                    ax2.plot([s-1,s],[gauge[s-2],gauge[s-1]],'b', linestyle=lstyle[k-1])

        #ax2.plot(range(1,self.numSlots+1),self.getGaugesFromCaster(1), 'b-')
        #if (self.numCasters == 2):
        #    ax2.plot(range(1,self.numSlots+1),self.getGaugesFromCaster(2), 'b--')

        # plot roller changes
        ax1.vlines(self.rollerChanges, minPriY, 0.91*maxPriY, colors='purple')
        
        
        ax1.set_ylim(minPriY, maxPriY)
        ax1.set_ylabel('Coil Order Width (inches)')
        ax2.set_ylim(minAltY, maxAltY)
        ax2.set_ylabel('Gauge (inches)', color='b')
        ax1.set_xlabel('Slot')

        self.getGradeFromCaster(1)
        listGrades = list(self.getGrades())
        nameGrades = {}
        cntGrade = 1
        for g in listGrades:
            if omitGrades:
                nameGrades[g]='Grade '+str(cntGrade)
                cntGrade += 1
            else:
                nameGrades[g]=g
        markers = []
        for g in listGrades:
            pos = [i+1 for i,x in enumerate(self.getGradeFromCaster(1)) if x == g]
            ax1.plot(pos,[0.93*maxPriY]*len(pos), 's', markersize=4, color=list(mcd.TABLEAU_COLORS)[listGrades.index(g)])
            if (self.numCasters == 2):
                pos = [i+1 for i,x in enumerate(self.getGradeFromCaster(2)) if x == g]
                ax1.plot(pos,[0.95*maxPriY]*len(pos), 's', markersize=4, color=list(mcd.TABLEAU_COLORS)[listGrades.index(g)])
            # set handle for legend
            mark = mlines.Line2D([], [], color=list(mcd.TABLEAU_COLORS)[listGrades.index(g)], marker='s', linestyle='None', markersize=4, label=nameGrades[g])
            markers.append(mark)
        ax1.legend(handles=markers, loc='upper right', ncol=len(listGrades), borderpad=0, frameon=False)

        if  saveFile:
            fig.savefig(saveFileName)
        plt.show()


    def writeSolution(self, filename):
        with open(filename, 'w') as fileData:
            writer = csv.writer(fileData, delimiter='\t')
            writer.writerow(['Caster', 'Slot', 'CoilID', 'Width'])
            for (k,s) in list(self.schedule):
                writer.writerow([k,s,self.schedule[(k,s)], self.castWidth[(k,s)]])
    
    def GetNovSolution(self, filename):
        """ read output to construct an array with its solution"""
        # header:
        with open(filename) as fileData:
            reader = csv.reader(fileData, delimiter='\t')
            for line in reader:
                if reader.line_num > 1:
                    self.assign(int(line[2]), int(line[0]), int(line[1]), float(line[3]))            
        
#%%
# Example
# sol = SteelccSolution()
# sol.loadSchoppyOutput('sch00.txt')
# sol.validateSolutionStrict()
# sol.plotSolution()
# plt.rcParams["figure.figsize"] = [15, 5]
# sol.getCosts()


# %%


# %%
