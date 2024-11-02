import pandas as pd
from collections import OrderedDict
import pandas as pd
import numpy as np
from datetime import datetime
from math import ceil
import itertools
from itertools import permutations
from datetime import date, datetime
import math
import random
import mmh3
from datetime import datetime
from bitarray import bitarray
import warnings
warnings.filterwarnings("ignore")
FALSEPROB = 0.001 # False positive probability for the Bloom Filter

# Define topologies for the network
# Order of stations for each line and topology
# A part of the network that is a line network, no common node ot ring
line_topology = ['Lls Centrum', 'Lls Schouwbrug', 'Lls Botterbrug','Lls Gondel Zuid', 'Lls Gondel Noord', 'Lls Punter','Lls Punter Zuid','Lls Wkc Jol', 'Lls Jol West', 'Lls Galjoenbrug', 'Lls Oostvaarders']
intersect_topology = ['Lls Centrum', 'Lls Agora','Lls Centrumbrug','Lls Lelycentre','Lls Vlieterbrug', 'Lls Zwb Koploper','Lls Oostkaap','Lls Atolplaza', 'Lls Bongerd','Lls Plantagebrug','Lls Geulbrug','Lls Gorsbrug']
ring_topology = ['Lls Centrum','Lls Ziekenhuis','Lls Politiebrug','Lls Slotenlaan','Lls Oostzeestr','Lls Zwartezeestr','Lls Duinbee','Lls Westhoven', 'Lls Bingerden','Lls Nienoord','Lls Horst', 'Lls Wkc Kamp','Lls Wold', 'Lls Kamp' ]

# User input for choosing topology
topology_choice = input("Which topology do you want to run (1: line_topology, 2: intersect_topology, or 3: ring_topology) ? ")

# Assign corresponding list to 'STATIONS' and 'LINE' based on user's choice
if topology_choice == '1':
    STATIONS = line_topology
    LINE = ['\\16002']
    GROUND_TRUTH = ['\\16002']

elif topology_choice == '2':
    STATIONS = intersect_topology
    LINE = ['\\16005']
    GROUND_TRUTH = ['\\16005']
elif topology_choice == '3':
    STATIONS = ring_topology
    LINE = ['\\16006']
    GROUND_TRUTH = ['\\16006']
else:
    print("Invalid choice. Defaulting to line_topology.")
    stations = line_topology


################################
##		Bloom Filter class	  ##
################################

def get_size(n, p):
	''' 
	Return the size of bit array(m) to used using
	following formula
	'''
	#return BFLEN
	m = -(n * math.log(p))/(math.log(2)**2)
	return int(m)
	
	
def get_hash_count(m, n):
	'''
	Return the hash function(k) to be used using
	following formula
	'''
	#return NHASH
	k = (m/n) * math.log(2)
	return int(k)

class BloomFilter ( object ):

	'''
		Class for Bloom filter, using murmur3 hash function
	'''
	def __init__(self, length, numOfHashes):
		
		'''
		items_count : int
			Number of items expected to be stored in bloom filter
		fp_prob : float
			False Positive probability in decimal
		'''
		# Size of bit array to use
		self.size = length
		# number of hash functions to use
		self.hash_count = numOfHashes
		# Bit array of given size
		self.bit_array = bitarray(self.size)
		# initialize all bits as 0
		self.bit_array.setall(0)
		
	def add(self, item):
		'''
		Add an item in the filter
		'''
		digests = []
		for i in range(self.hash_count):
			# create digest for given item.
			# i work as seed to mmh3.hash() function
			# With different seed, digest created is different
			digest = mmh3.hash(item, i) % self.size
			digests.append(digest)
			# set the bit True in bit_array
			self.bit_array[digest] = True
		return self.bit_array

	def check(self, item):
		'''
		Check for existence of an item in filter
		'''
		for i in range(self.hash_count):
			digest = mmh3.hash(item, i) % self.size
			if self.bit_array[digest] == False:
				return False
		return True
###################################
##### calling Bloom filter   ######
###################################

def transfer_to_bfs(card_ids_dict):
    card_id_bfs = {}
    for key, values in card_ids_dict.items():
        bloomf = BloomFilter(BFLEN, NHASH)
        for value in values:
            bloomf.add(str(value))
        card_id_bfs[key] = bloomf.bit_array

    return card_id_bfs
#################################################
###                intersetion                ### 
#################################################
def bitwise_and(a, b):
    result_arr = a.copy()
    result_arr &= b
    return result_arr
#################################################
###                Union                      ### 
#################################################
def bitwise_or(a, b):
    result_or = a.copy()
    result_or |= b
    return result_or
#################################################
###                Cardinality                ### 
#################################################

def cardinality(union_lst):
    t = sum(int(bit) for bit in union_lst)
    return int(-1 * (BFLEN / NHASH) * math.log(1 - t/BFLEN))
#####################################
####		 Accuracy			 ####
#####################################	
def calculate_accuracy(gt, estimated):
    accuracy = []
    for i in range(len(gt)):
        if gt[i] == 0:
            if estimated[i] == 0:
                acc = 1.0  # Perfect accuracy if both are zero
            else:
                acc = 0.0  # No accuracy if gt is zero but estimated is not
        else:
            acc = max(1 - (abs(estimated[i] - gt[i]) / gt[i]), 0)
        accuracy.append(round(acc, 2))  # Format the accuracy value to two decimal places
    return accuracy
##################################################################################
# Creating two check in and check out table for distinction between check in /out
##################################################################################
def split_df_by_table_in_table_out(filtered_routes):
    # This table selects only the columns related to check-in (card_ID, check-in time, location of check-in, and the route)
    check_in_table = filtered_routes[['card_ID', 'check_in', 'loc1', 'Routes']]
    check_in_table = check_in_table.rename({"card_id": "card_ID", "loc1": "station", "check_in": "epoch", 'Routes':"line"}, axis=1)
    
    # Create the check-out table
    check_out_table = filtered_routes[['card_ID', 'check_out', 'loc2', 'Routes']]
    check_out_table = check_out_table.rename({"card_id": "card_ID", "loc2": "station", "check_out": "epoch", 'Routes':"line"}, axis=1)
    # Combine the check-in and check-out tables into a dictionary for easy access
    dfs_in_out = {'check_in': check_in_table, 'check_out': check_out_table}    # The dictionary keys 'check_in' and 'check_out' can be used to retrieve the respective tables
    return dfs_in_out
####################################################################################
# Transfering card ids to bloom filters for each station
####################################################################################
def extract_card_ids(STATIONS, filtered_df):
    card_ids_dict = {}

    for station in STATIONS:
        # Check which columns are present in the DataFrame. This is important to determine 
        # how to filter the DataFrame based on the available information.
        has_loc1 = 'loc1' in filtered_df.columns
        has_loc2 = 'loc2' in filtered_df.columns
        # If both 'loc1' and 'loc2' columns are present, it indicates a detection scenario.
        # In this case, extract card IDs where the station is mentioned in either 'loc1' or 'loc2'.
        if has_loc1 and has_loc2:
            # For Detection scenario: Both loc1 and loc2 are present
            card_ids = filtered_df[(filtered_df['loc1'] == station) | (filtered_df['loc2'] == station)]['card_ID'].tolist()
        else:
            # If only 'loc1' or 'loc2' is present, it indicates a scenario where we need to distinguish 
            card_ids = filtered_df[filtered_df['station'] == station]['card_ID'].tolist()
        card_ids_dict[station] = card_ids

    return card_ids_dict 
###########################################################################
#                        Size based on GT:
###########################################################################
def calculate_checkout_counts(df, STATIONS):
    checkout_counts = []
    
    for index, station in enumerate(STATIONS[1:], start=1):
        # Get all previous stations as a list of possible check-in locations
        possible_check_in_stations = STATIONS[:index]
        # Filter dataframe for rows where loc2 is the current station
        # and loc1 is in the list of possible check-in stations
        filtered_df = df[(df['loc2'] == station) & (df['loc1'].isin(possible_check_in_stations))]
        # Count the number of check-outs at the current station
        checkout_count = len(filtered_df)
        # Store the count in the dictionary with the station as the key
        checkout_counts.append(checkout_count)

    return checkout_counts

########################################################
########################################################
                   # Main
########################################################
########################################################

df = pd.read_csv('lelystad_data.csv', sep=';')
travelerIDSet = random.sample(range(100000000000000), len(df)) # Generate random traveler IDs
df['card_ID']=travelerIDSet # Assign random IDs to 'card_ID' column

# Filter DataFrame for rows where either 'loc1' or 'loc2' matches the stations in the selected topology
filtered_df = df[df['loc1'].isin(STATIONS) | df['loc2'].isin(STATIONS)] # looking at nodes we are interested not the whole nodes.
########### # Initialize Bloom Filter parameters###################
BFLEN = get_size(len(df), FALSEPROB)
NHASH = get_hash_count(BFLEN, len(df))
########################### Calculate ground truth (GT) for the checkout counts at each station###################################
filter_line = filtered_df[filtered_df['Routes'].isin(LINE)]
GT = calculate_checkout_counts(filter_line, STATIONS)
print("GT:",GT)
##############################################################################
# First Scenario Results: we only have Detections  (A ∪ B ∪ C ∪ D ∪ E) ∩ ( F )
##############################################################################
card_ids_dict = extract_card_ids(STATIONS,filtered_df)# Extract card IDs and transfer them to Bloom filters for each station
card_id_bfs= transfer_to_bfs(card_ids_dict)

output_IDS_bfs = []
for station in reversed(STATIONS):  # Using "F" as a starting point, intersect it with the union of all stations before it
    method_1=[]
    union = bit_bfs = bitarray([False] * BFLEN)
    for prev_station_bfs in STATIONS[:STATIONS.index(station)]:
        union = bitwise_or(union, card_id_bfs[prev_station_bfs]) # (A ∪ B ∪ C ∪ D ∪ E)
    inter_bfs= bitwise_and(card_id_bfs[station], union)     # Intersection of the current station's Bloom filter with the union of previous stations
    method_1.append(cardinality(inter_bfs))
        
    if station==STATIONS[0]: # If the first station is reached, stop the loop because it is first source 
        break        
    output_IDS_bfs.extend(method_1) 
output_IDS_bfs = output_IDS_bfs[::-1] # Reverse the list to maintain the correct order
        
print("Estimated size first scenario: " ,output_IDS_bfs)
first_result= calculate_accuracy(GT, output_IDS_bfs)
print("First scenario Accuracy: ", first_result)
##############################################################################
# Second Scenario Results: we have distinction between check in and out  
##############################################################################
tables_dict= split_df_by_table_in_table_out(filtered_df) # Split the filtered data into separate tables for check-ins and check-outs
check_in_table = tables_dict['check_in'] # Separate table for check-ins
check_out_table = tables_dict['check_out']  
card_ids_dict_in = extract_card_ids(STATIONS,check_in_table)# Extract card IDs from check-in and check-out tables and convert them to Bloom filters
card_ids_dict_out = extract_card_ids(STATIONS,check_out_table)
card_id_bfs_in= transfer_to_bfs(card_ids_dict_in)# Convert the lists of card IDs to Bloom filters for both check-in and check-out
card_id_bfs_out= transfer_to_bfs(card_ids_dict_out)

output_IDS_bfs_2 = []
for station in reversed(STATIONS):  
    scenario_2=[]
    union = bit_bfs = bitarray([False] * BFLEN)
    for prev_station_bfs in STATIONS[:STATIONS.index(station)]:
        union = bitwise_or(union, card_id_bfs_in[prev_station_bfs]) 

    inter_bfs_2= bitwise_and(card_id_bfs_out[station], union)
    scenario_2.append(cardinality(inter_bfs_2))
        
    if station==STATIONS[0]:
        break        
    output_IDS_bfs_2.extend(scenario_2) 
output_IDS_bfs_2 = output_IDS_bfs_2[::-1]
        
print("Estimated size second scenario: " ,output_IDS_bfs_2)
scenod_scenario_result= calculate_accuracy(GT, output_IDS_bfs_2)
print("second scenario Accuracy: ", scenod_scenario_result)
##############################################################################
# Third Scenario Results: we have Lines, distinction between check in and out  
############################################################################## 
# Filtering the check-in and check-out tables based on the specific line chosen by the user
check_in_table_line = check_in_table[check_in_table['line'].isin(LINE)]# Filter the check-in table for records in the selected LINE
check_out_table_line = check_out_table[check_out_table['line'].isin(LINE)]
card_ids_dict_in_line = extract_card_ids(STATIONS,check_in_table_line)
card_ids_dict_out_line = extract_card_ids(STATIONS,check_out_table_line)
card_id_bfs_in_line= transfer_to_bfs(card_ids_dict_in_line)
card_id_bfs_out_line= transfer_to_bfs(card_ids_dict_out_line)

output_IDS_bfs_3 = []
for station in reversed(STATIONS):  # Using "F" as a starting point, intersect it with the union of all stations before it
    scenario_3=[]
    union = bit_bfs = bitarray([False] * BFLEN)
    for prev_station_bfs in STATIONS[:STATIONS.index(station)]:
        union = bitwise_or(union, card_id_bfs_in_line[prev_station_bfs]) # (A ∪ B ∪ C ∪ D ∪ E)

    inter_bfs_3= bitwise_and(card_id_bfs_out_line[station], union)
    scenario_3.append(cardinality(inter_bfs_3))
        
    if station==STATIONS[0]:
        break        
    output_IDS_bfs_3.extend(scenario_3) 
output_IDS_bfs_3 = output_IDS_bfs_3[::-1]
        
print("Estimated size third scenario: " ,output_IDS_bfs_3)
third_scenario_result= calculate_accuracy(GT, output_IDS_bfs_3)
print("third scenario Accuracy: ", third_scenario_result)





