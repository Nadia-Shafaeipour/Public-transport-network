# Lelystad_paper
Balancing Precision and Privacy: Using Bloom Filters and Homomorphic Encryption to Count Travelers

# Introduction
This document explains a Python script designed for analyzing public transportation networks using Bloom Filters. The script evaluates different network topologies - line, intersect, and ring  under three different scenarios to count the movement of travelers between nodes. 
Code Structure
The script is structured as follows:
- Importing Libraries: Essential Python libraries such as pandas, numpy, and others are imported.
- Network Topologies Definition: Three types of topologies are defined - line, intersect, and multi paths.
- User Input for Topology Selection: The script prompts the user to select a topology for analysis.
- Bloom Filter Implementation: Functions and a class to implement the Bloom Filter data structure are defined.
- Main Script: The core logic for data processing, applying Bloom Filters, and scenario analysis.
- Scenario-Based Analysis: The script performs different analyses for each scenarios.
-Accuracy of counting method: To evaluate the counting method's accuracy, we compare the estimated number of travelers derived from the method against the actual counts, referred to as the ground truth, for each scenario.

# Detailed Explanation
1. Importing Libraries: The script starts by importing required libraries. These libraries provide functionalities for data manipulation (pandas), mathematical operations (math, numpy), and specific data structures (bitarray, collections).

2. Network Topologies: Three arrays define the station names for each of the three network topologies. These arrays are used to simulate different network structures and the order of stations on Lelystad map.

3. User Input: The script asks the user to choose one of the three predefined topologies for analysis.

4. Bloom Filter Implementation: The Bloom Filter, a space-efficient probabilistic data structure, is implemented using murmur3 hash function. It is used to estimate the presence of card identifiers of travelers with a controlled error rate.

5. Main Script Logic: After reading the input data, the script filters it based on the selected topology. Bloom Filters are then used to estimate the number of travelers at each station. This section contains the core logic for processing the data and applying the Bloom Filters.

6. Scenario-Based Analysis: The script runs three different analyses for each scenario. These scenarios include simple detection of presence at stations, distinguishing between check-ins and check-outs, and considering specific lines in the network.

# How to Use
running the Jupyter notebook file "lelystad_network.ipynb".
View Results: Observe the output that includes estimated sizes at each station per scenario and accuracy measurements.
