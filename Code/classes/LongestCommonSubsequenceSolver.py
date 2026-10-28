from Profile import *
from BinnedHistory import *
from Bin import *

class LongestCommonSubsequenceSolver:
    """
    Class with functions for aligning patients to each other or to a profile.
    """
    def __init__(self, MATCH_REWARD=-2, MISMATCH_PENALTY=1, INDEL_PENALTY=0.1):
        """
        Construct a 'LongestCommonSubsequenceSolver' object. 

        Parameters:
        MATCH_REWARD (int or float): Reward for matching events. Should be negative.
        MISMATCH_PENALTY (int or float): Penalty for an event not occuring in both bins. Should be positive.
        INDEL_PENALTY (int or float): -INDEL_PENALTY is used as the distance of two empty bins. Should be positive.
        """
        self.MATCH_REWARD = MATCH_REWARD
        self.MISMATCH_PENALTY = MISMATCH_PENALTY
        self.INDEL_PENALTY = INDEL_PENALTY
    
    def getEmptyMatrix(self, patient1_last_bin_index):
        """
        Constructs an empty matrix (a dictionary of dictionaries).

        Parameter:
        patient1_last_bin_index (int): length of patient1's treatment in bin_lengths
        """
        Matrix = {}
        # generate other rows with a zero in the 0th position
        for i in range(patient1_last_bin_index+1):
            Matrix[i] = {}
        return Matrix

    def initializeDirectionMatrix(self, history1_length, history2_length):
        """
        Constructs an empty direction matrix and fills in top row and leftmost column.
        """
        Matrix = self.getEmptyMatrix(history1_length)   
        
        numColumns = history2_length + 1
        for i in range(numColumns):
            Matrix[0][i] = 2

        numRows = history1_length + 1
        for i in range(numRows):
            Matrix[i][0] = 1

        return Matrix
    
    def initializeDistanceMatrixWithProfile(self, history1:BinnedHistory, profile:Profile):
        """
        Constructs an empty distance matrix for comparing a binned history and profile,
        and fills in first row and column. 
        """
        Matrix = self.getEmptyMatrix(history1.get_length()) 
        numColumns = profile.get_length() + 1
        
        Matrix[0][0] = 0
        for i in range(numColumns):
            if i == 0:
                continue
            prev_value = 0
            if i > 0:
                prev_value = Matrix[0][i-1]
            penalty = profile.get_total_events(i-1)*self.MISMATCH_PENALTY 
            if penalty == 0:
                penalty = self.INDEL_PENALTY # no events missed
            Matrix[0][i] = prev_value + penalty
            
        numRows = history1.get_length() + 1           
        for i in range(numRows):
            if i==0:
                continue
            prev_value = 0
            if i > 0 :
                prev_value = Matrix[i-1][0]

            penalty = len(history1.get_bin(i-1).event_names)*self.MISMATCH_PENALTY
            if penalty == 0:
                penalty = self.INDEL_PENALTY # no events missed
            Matrix[i][0] = prev_value + penalty

        return Matrix
    
    def initializeDistanceMatrix(self, history1:BinnedHistory, history2:BinnedHistory):
        """
        Constructs an empty distance matrix for comparing two binned histories,
        and fills in first row and column. 
        """
        Matrix = self.getEmptyMatrix(history1.get_length())
        numColumns = history2.get_length() +1
        Matrix[0][0] = 0

        for i in range(numColumns):
            if i == 0:
                continue
            prev_value = 0
            if i > 0 :
                prev_value = Matrix[0][i-1]
            penalty = len(history2.get_bin(i-1).event_names)*self.MISMATCH_PENALTY
            if penalty == 0:
                penalty = self.INDEL_PENALTY # no events missed
            Matrix[0][i] = prev_value + penalty

        numRows = history1.get_length() + 1      
        for i in range(numRows):
            if i==0:
                continue
            prev_value = 0
            if i > 0 :
                prev_value = Matrix[i-1][0]

            penalty = len(history1.get_bin(i-1).event_names)*self.MISMATCH_PENALTY
            if penalty == 0:
                penalty = self.INDEL_PENALTY # no events missed
            Matrix[i][0] = prev_value + penalty

        return Matrix
    
    def printMatrix(self, matrix):
        """
        Prints a matrix row by row in the form index:value 
        """
        for index, row in matrix.items():
            rowString = ""
            for index, value in row.items():
                rowString += f"{index}:{round(value, 2)}\t"
            print(rowString)

    def alignWithProfile(self, history1:BinnedHistory, profile:Profile):
        """
        Aligns a binned history with a profile. 

        Returns:
        alignment: aligned binned history
        profile_indices: list of indices and indels ("-") representing how the position of 
        aligned profile looks relative to the original profile.
        distance: distance between the binned history and profile
        """
        distanceMatrix = self.initializeDistanceMatrixWithProfile(history1, profile) 
        directionMatrix = self.initializeDirectionMatrix(history1.get_length(), profile.get_length())
        profile_length = profile.get_length()

        for i in range(1, history1.get_length()+1):  
            columnRange = (1, profile.get_length()+1)
            for j in range(columnRange[0], columnRange[1]):  
                
                bin1 = history1.get_bin(i-1)
                match = profile.get_distance(j-1, bin1, self.MATCH_REWARD, self.MISMATCH_PENALTY, self.INDEL_PENALTY)
  
                down_penalty = len(history1.get_bin(i-1).event_names)*self.MISMATCH_PENALTY
                if down_penalty == 0:
                    down_penalty = self.INDEL_PENALTY

                right_penalty = profile.get_indel_penalty(j-1, self.MISMATCH_PENALTY, self.INDEL_PENALTY)

                distanceMatrix[i][j] = min([(distanceMatrix[i-1][j] + down_penalty),
                                            (distanceMatrix[i][j-1] + right_penalty), 
                                            (distanceMatrix[i-1][j-1] + match)])
                    
                if distanceMatrix[i][j] == (distanceMatrix[i][j-1] + right_penalty):
                    directionMatrix[i][j] = 2
                elif distanceMatrix[i][j] == (distanceMatrix[i-1][j] + down_penalty):
                    directionMatrix[i][j] = 1
                elif distanceMatrix[i][j] == (distanceMatrix[i-1][j-1] + match):
                    directionMatrix[i][j] = 0

        distance = self.make_distance_relative(unadjusted_distance=distanceMatrix[history1.get_length()][profile_length], history1=history1, history2=profile)
        alignment, profile_indicies = self.backtrack_with_profile(history1, profile, directionMatrix)
        
        return alignment, profile_indicies, distance
    
    def make_distance_relative(self, unadjusted_distance, history1:BinnedHistory, history2):
        """
        Makes distance relative using num events and gaps in the two histories.
        """
        if type(history2) == BinnedHistory:
            total_events = history1.get_num_events() + history2.get_num_events()
            total_gaps = history1.get_num_gaps() + history2.get_num_gaps()
            max_distance = total_events*self.MISMATCH_PENALTY + total_gaps*(self.INDEL_PENALTY)
            min_distance = total_events*self.MATCH_REWARD + total_gaps*(-self.INDEL_PENALTY)
            if max_distance == min_distance:
                print("error: max_distance == min_distance")
            distance = abs(min_distance-unadjusted_distance)/(max_distance-min_distance)
            if round(distance, 3) > 1.0:
                print(f"error adjusting distance: {round(distance, 3)} > 1.0")
            # print(f"max_distance: {max_distance}, min_distance: {min_distance}, unadjusted_distance: {unadjusted_distance}, distance: {distance}")
            
        elif type(history2) == Profile:
            h1_events = history1.get_num_events()
            h1_gaps = history1.get_num_gaps()
            min_distance = history2.get_min_dist(self.MATCH_REWARD, self.MISMATCH_PENALTY, self.INDEL_PENALTY) + h1_events*self.MATCH_REWARD + h1_gaps*-self.INDEL_PENALTY
            max_distance = history2.get_max_dist(self.MATCH_REWARD, self.MISMATCH_PENALTY, self.INDEL_PENALTY) + h1_events*self.MISMATCH_PENALTY + h1_gaps*self.INDEL_PENALTY
            distance = abs(min_distance-unadjusted_distance)/(max_distance-min_distance)
            # if distance > 1:
            #     print(f"error adjusting distance: {round(distance, 3)} > 1")

        return distance
    
    def getDistanceandAlignment(self, history1:BinnedHistory, history2:BinnedHistory, returnAlignment=False):
        """
        Aligns two binned histories. 

        Parameters:
        history1 (BinnedHistory): first binned history
        history2 (BinnedHistory): second binned history
        returnAlignment (bool): if True, returns the aligned binned history

        if returnAlignment, returns alignment1, alignment2, distance
        otherswise, return distance
        """
        if history1 == history2:
            return 0

        distanceMatrix = self.initializeDistanceMatrix(history1, history2) 
        directionMatrix = self.initializeDirectionMatrix(history1.get_length(), history2.get_length())
        
        for i in range(1, history1.get_length()+1):   

            # determine columns to fill in for given row
            columnRange = (1, history2.get_length()+1)
    
            for j in range(columnRange[0], columnRange[1]):  
                
                bin1 = history1.get_bin(i-1)
                bin2 = history2.get_bin(j-1)
                
                # assume mismatch
                match = None

                # neither bins are blank
                if (bin1 != False) and (bin2 != False):
                    # match distance will be determined by comparing the content of the 2 bins
                    match = bin1.get_distance(bin2, self.MATCH_REWARD, self.MISMATCH_PENALTY, self.INDEL_PENALTY)
                
                down_penalty = len(history1.get_bin(i-1).event_names)*self.MISMATCH_PENALTY
                if down_penalty == 0:
                    down_penalty = self.INDEL_PENALTY
                right_penalty = len(history2.get_bin(j-1).event_names)*self.MISMATCH_PENALTY
                if right_penalty == 0:
                    right_penalty = self.INDEL_PENALTY
                
                distanceMatrix[i][j] = min([(distanceMatrix[i-1][j] + down_penalty),
                                        (distanceMatrix[i][j-1] + right_penalty), 
                                        (distanceMatrix[i-1][j-1] + match)])
                if returnAlignment:
                    if distanceMatrix[i][j] == (distanceMatrix[i][j-1] + right_penalty):
                        directionMatrix[i][j] = 2
                    elif distanceMatrix[i][j] == (distanceMatrix[i-1][j] + down_penalty):
                        directionMatrix[i][j] = 1
                    elif distanceMatrix[i][j] == (distanceMatrix[i-1][j-1] + match):
                        directionMatrix[i][j] = 0

        distance = self.make_distance_relative(unadjusted_distance=distanceMatrix[history1.get_length()][history2.get_length()], history1=history1, history2=history2)
    
        if returnAlignment:
            alignment1, alignment2 = self.backtrack(history1, history2, directionMatrix)
            return alignment1, alignment2, distance
        else:
            return distance

    def backtrack_with_profile(self, history1:BinnedHistory, profile:Profile, directionMatrix):
        """
        Backtracks through direction matrix to get aligned binned history and profile indicies

        Returns:
        alignment: history1 aligned
        profile_indicies: list of indices and indels ("-") representing how the position of 
        aligned profile looks relative to the original profile.
        """
        profile_indicies = []
        history1_aligned_bins = []
        i = history1.get_length()
        j = profile.get_length()

        while i != 0 or j != 0: # continue the loop until we reach [0,0]
            direction = directionMatrix[i][j]
            if direction == 1:
                # go up, add character from s2 to alignment2 and INDEL to alignment1
                profile_indicies.insert(0, "-")
                history1_aligned_bins.insert(0, history1.get_bin(i-1))
                i -= 1
            elif direction == 2:
                indel_bin = Bin()
                indel_bin.event_names.add("-")
                # go left, add character from s1 to alignment1 and INDEL to alignment2
                history1_aligned_bins.insert(0, indel_bin)
                profile_indicies.insert(0, profile.get_bin(j-1))
                j -= 1
            elif direction == 0:
                # match (or SUB), go diagonally up and add character to alignment1 and 2
                profile_indicies.insert(0, profile.get_bin(j-1))
                history1_aligned_bins.insert(0, history1.bins[i-1])
                i -= 1
                j-= 1
        
        alignment = BinnedHistory(bins=history1_aligned_bins)
        return alignment, profile_indicies
    
    def backtrack(self, history1:BinnedHistory, history2:BinnedHistory, directionMatrix):
        """
        Backtracks through direction matrix to get both aligned binned histories 

        Returns:
        alignment: history1 aligned
        profile_indicies: list of indices and indels ("-") representing how the position of 
        aligned profile looks relative to the original profile.
        """        
        history1_aligned_bins = []
        history2_aligned_bins = []
        i = history1.get_length()
        j = history2.get_length()
        indel_bin = Bin()
        indel_bin.add_event(Event("-", 0, 0))
        
        while i != 0 or j != 0: # continue the loop until we reach [0,0]
            direction = directionMatrix[i][j]
            if direction == 1:
                # go up, add character from s2 to alignment2 and INDEL to alignment1
                history2_aligned_bins.insert(0, indel_bin) # TODO: use append instead if needed to speed things up
                history1_aligned_bins.insert(0, history1.bins[i-1])
                
                i -= 1
            elif direction == 2:
                # go left, add character from s1 to alignment1 and INDEL to alignment2
                history1_aligned_bins.insert(0, indel_bin)
                history2_aligned_bins.insert(0, history2.bins[j-1])
                j -= 1
            elif direction == 0:
                # match (or SUB), go diagonally up and add character to alignment1 and 2
                history2_aligned_bins.insert(0, history2.bins[j-1])
                history1_aligned_bins.insert(0, history1.bins[i-1])
                i -= 1
                j-= 1
        
        alignment1 = BinnedHistory(bins=history1_aligned_bins)
        alignment2 = BinnedHistory(bins=history2_aligned_bins)

        return alignment1, alignment2
