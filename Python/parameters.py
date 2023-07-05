# Maximum allowable trim loss for CUT edge type [in]
MAX_TRIMLOSS_CUT = 5 + 0.5
# Maximum allowable trim loss for HRB edge type [in]
MAX_TRIMLOSS_HRB = 0.5
# Maximum allowable trim loss for MILL edge type [in]
MAX_TRIMLOSS_MILL = 0.5
# ALTERNATIVE Maximum allowable trim loss for MILL edge type [in]
MAX_TRIMLOSS_MILL_ALTERNATIVE = 6 + 0.5
# ALTERNATIVE MINIMUM  allowable trim loss for MILL edge type [in]
MIN_TRIMLOSS_MILL_ALTERNATIVE = 1.5

# Maximum cast width decrease between adjacent coils [in]
MAX_CASTER_DECREASE = 1.5 + 1.73
# Maximum difference between cast width of casters [in]
MAX_EXTRA_INTERCASTER_WIDTH = 4
# Allowable difference in widths between casters without penalty [in]
ALLOWABLE_INTERCASTER_WIDTH = 2

# Threshold gauge to be considered a heavy, medium or light coil
HEAVY_GAUGE = 0.4
MEDIUM_GAUGE = 0.123

# Minimum gauge for a new roller
MIN_GAUGE_ROLLER_CHANGE = 0.155
# Maximum total degradation
MAX_ROLLER_DEGRADATION = 100

# Grade cost parameters
COST_FOR_FORBIDDEN_GRADES = 70
MIN_HEAT_SIZE = 150
MIN_GRADE_CHANGE_PENALTY = 0.066666 

# Roller cost
COST_ROLLER_CHANGE = 13.33333

# CHECK TOLERANCE
CHECK_TOLERANCE = 1e-4


# Penalty for width differences between casters [$/in]
PENALTY_WIDTH_DIFFERENCE = 6.66666
# Penalty for width differences trim loss [$/in]
PENALTY_TRIM_LOSS = 1
# Penalty for roll change [$]
PENALTY_ROLL_CHANGE = 13.33333
# Penalty for excessive gauge decrease
PENALTY_GAUGE_DECREASE = 166.66666

# Heat wight limit
MAX_HEAT_WEIGHT = 175
MIN_HEAT_WEIGHT = 140

