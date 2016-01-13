# csToGraph - CellSet to Graph

# A tool for converting CellSet segmentation into a network representation.

# Running:

# ./python csTograph segmentation.xml

# Author:
# M.Delmans

import sys
import xml.etree.ElementTree as et
import networkx as nx
from networkxgmml import XGMMLWriter

inFile = sys.argv[1]
outFileS = inFile.split('.xml')[0] + '.xgmml'
edgeFileS = inFile.split('.xml')[0] + '.edgeLst'  

tree = et.parse(inFile)
outFile = open(outFileS, 'w')

root = tree.getroot()

xmlCells = root[1]

nCells = int(xmlCells.attrib['count'])

walls = dict()

for cell in xmlCells:
	
	cellId = int(cell.attrib['id'])
	
	for wall in cell[0]:
		
		wallId = int(wall.attrib['id'])
		
		if not wallId in walls:
			walls[wallId] = []

		walls[wallId].append(cellId)


g = nx.Graph()

g.add_nodes_from( range(nCells) )

for wall in walls:
	if len(walls[wall]) == 2:
		g.add_edge( walls[wall][0], walls[wall][1] )


XGMMLWriter(outFile, g, 'Gemma')

nx.write_edgelist(g, edgeFileS)
