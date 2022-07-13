# KEGGMutli
# Language: Python
# Input: TXT
# Output: SCREEN
# Tested with: PluMA 1.1, Python 3.6
# Dependency: BioPython 1.79

PluMA plugin to take a CSV file of multiomics data (microbe and metabolite) and report
connections that appear in the KEGG (Kanehisa and Goto, 2000) database.

The plugin accepts as input a tab-delimited file of keyword-value pairs:
csvfile: Network of connections.  Names of microbes and metabolites are assumed to be in the header
keggfile: CSV file that maps kegg identifiers to microbes and metabolites in the CSV

Connections (other than microbe-microbe) that have a positive edge weight in the network
and have at least one pathway in the KEGG database will be output to the screen.
