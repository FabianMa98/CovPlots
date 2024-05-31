"""
Constants 
"""

# Debating whether this should just be an option for argparsing
WINDOW = 2000
MUTATION_SCALE = 30
ALPHA = 0.8
LINEWIDTH = 1.5

# CONSTANT FOR LINE LOC LOCATION
# This implies that the maximum coverage in our window should be higher than 6.5
# Throw error if this is not the case (for e.g with message: No significant coverage found?)
FACTOR = 6.5
