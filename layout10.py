from gdsCAD import *
import shapely
import shapely.ops
import shapely.affinity
from shapely.geometry import *
import os
import numpy as np
import math
import copy
import operator
from descartes import PolygonPatch
from lxml import etree as ET
from materials import Materials
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.text
import matplotlib.lines
import matplotlib.transforms as transforms
import matplotlib.cm
from matplotlib.backends.backend_pdf import PdfPages
from random import shuffle

pp = PdfPages('output.pdf')
cellListNode = ET.Element("CellList")

def myShow(element,name):
	"""
	Display the object
	
	:returns: the display Axes.
	
	"""
	fig = plt.figure()
	ax = plt.gca()
	ax.set_aspect('equal')
	ax.set_title(name)
	ax.margins(0.1)	   
	
	artists=element.artist()
	for a in artists:
		a.set_transform(a.get_transform() + ax.transData)
		if isinstance(a, matplotlib.patches.Patch):
			ax.add_patch(a)
		elif isinstance(a, matplotlib.lines.Line2D):
			ax.add_line(a)
		else:
			ax.add_artist(a)
	
	ax.autoscale(True)
	plt.draw()
	pp.savefig()
	plt.close()
	# plt.show()
	return ax

def p2npa(point):
	return np.array([point.x,point.y])

def cartCoord(rad,ang):
	return np.array([rad*math.cos(ang),rad*math.sin(ang)])

def pieShape(radius,ang,resolution):
	points = [np.array([0,0])]
	for i in range (0,resolution+1):
		points.append(cartCoord(radius,float(i)/float(resolution)*ang+math.pi/2))
	return Polygon(points)

def starShape(radius,thickness,number):
	shape = Point(0,0)
	for i in range (0,number):
		point = cartCoord(radius,float(i)/float(number)*math.pi)
		line = LineString([(point[0],point[1]),(-point[0],-point[1])])
		shape = shape.union(line.buffer(thickness,cap_style=2))
	return shape

def uniteReflection(shape):
	return shape.union(shapely.affinity.scale(shape, xfact=-1.0,origin=(0,0)))
	
def distance(pointA,pointB):
	return math.hypot(pointB[0] - pointA[0], pointB[1] - pointA[1])

def boundingBoxPoints(shape):
	boundingBox = shape.bounding_box
	return [boundingBox[0],[boundingBox[0][0],boundingBox[1][1]],boundingBox[1],[boundingBox[1][0],boundingBox[0][1]]]

def shape_in_shape(defaultCellObj,position,shapeB):
	poly = Polygon(shapeB.points)
	polyB = shapely.affinity.translate(defaultCellObj.getBoundingBox(),position[0],position[1])
	return poly.contains(polyB)

def shapeIntersect(defaultCellObj,position,shapeB):
	polyA = shapely.affinity.translate(defaultCellObj.getBoundingBox(),position[0],position[1])
	polyB = Polygon(boundingBoxPoints(shapeB))
	return polyA.intersects(polyB)
	
def getHollowShapePoints(shape,layer=None):
	element = core.Elements()
	try:
		shape_iter = iter(shape)
		# print 'is iterable'
		for geom in shape.geoms:
			element.add(getGDSpolygon(geom,layer))
		return element
	except TypeError, te:
		# print 'is not iterable'
		return getGDSpolygon(shape,layer)

def getGDSpolygon(shape,layer=None):
	element = core.Elements()
	points = shape.exterior.coords
	for interior in shape.interiors:
		pointsB = interior.coords
		index = min(range(len(pointsB)), key=lambda i: distance(points[0],pointsB[i]))
		pointsB = [pointsB[i - (len(pointsB)-index)] for i in range(len(pointsB))]
		points = np.vstack((points,pointsB,pointsB[0],points[0]))
	return core.Boundary(points,layer)

def labelRange(label1,label2):
	labels = []	
	j1 = ord(label1[0]) - ord('A')
	i1 = ord(label1[1]) - ord('A')
	j2 = ord(label2[0]) - ord('A')
	i2 = ord(label2[1]) - ord('A')
	imin = min(i1,i2)
	imax = max(i1,i2)
	jmin = min(j1,j2)
	jmax = max(j1,j2)
	for i in range(imin,imax+1):
		for j in range(jmin,jmax+1):
			labels.append(chr(ord('A')+j)+chr(ord('A')+i))
	return labels

def polyLabelRange(labelCoupleList):
	labels = []
	for labelCouple in labelCoupleList:
		if isinstance(labelCouple,basestring):
			labels.extend([labelCouple])
		else:
			labels.extend(labelRange(labelCouple[0],labelCouple[1])) #Make label range on the line from "DB"
	return list(set(labels))
	
def labelAdd(label1,(i,j)):
	j1 = ord(label1[0]) - ord('A')
	i1 = ord(label1[1]) - ord('A')
	j2 = j1+j
	i2 = i1+i
	return chr(ord('A')+j2)+chr(ord('A')+i2)

def labelRangeFromPos(label1,(i,j)):
	return labelRange(label1, labelAdd(label1,(i,j)))

def symetricLabel(label,colNum):
	j1 = ord(label[0]) - ord('A')
	i1 = colNum-1-ord(label[1]) + ord('A')
	return chr(ord('A')+j1)+chr(ord('A')+i1)
		
def fill(labelList,cellList,random=False,symetrise=False,colNum=0):
	specialCellList = []
	def myRand():
		return 0.54321
	if random: shuffle(cellList,myRand)
	if len(cellList)>len(labelList):
		print "Warning: the label list is too small: somes cells will not be assigned"
		print str(len(cellList))+"cells for "+str(len(labelList))+" positions"
	print labelList
	for label in  labelList:
		i = labelList.index(label) % len(cellList)
		specialCellList.append((label,cellList[i]))
		if symetrise:
			specialCellList.append((symetricLabel(label,colNum),cellList[i]))
	return specialCellList
	
LAYOUT_layer = 0
DRIE1_layer = 1
DRIE2_layer = 2
SPUT1_layer = 3
SPUT2_layer = 4
OTS_MASK_layer = 6
CHANETCH_layer = 5
GLASSETCH_layer = 7

class Cavity(object):
	position = None
	packingWidth = 300
	packing = True
	outsiteLength = 0
	totalArea = None
	refCell = None
	shape = None
	
	def __init__(self, refCell):
		self.refCell = refCell

class DispCav(Cavity):
	name="SquareDispenserCavity"
	dispCavSize = 1600
	
	def make(self):
		cell = core.Cell(self.name)
		self.totalArea = 0
		self.outsiteLength = 0
		# dispCav = core.Elements()
		dispCavShape = Point((0,0)).buffer(self.dispCavSize/2,resolution=1).envelope
		self.outsiteLength = self.outsiteLength + dispCavShape.exterior.length
		self.totalArea = self.totalArea + dispCavShape.area
		if self.packing:
			dispCavShape = dispCavShape.difference(dispCavShape.buffer(-self.packingWidth))
		self.shape = dispCavShape
		dispCav = getHollowShapePoints(dispCavShape, layer=DRIE1_layer)
		dispCav = utils.translate(dispCav,self.position)
		cell.add(dispCav)
		
		return cell, self.toElement()

	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("size",str(self.dispCavSize))
		root.set("perimeter",str(self.outsiteLength))
		root.set("area",str(self.totalArea))
		return root

class OptCav(Cavity):
	name="RoundOpticalCavity"
	radius = 1000
	optCavArea = 0
	otsMaskOffset = 500

	def make(self):
		cell = core.Cell(self.name)
		optCavShape = Point((0,0)).buffer(self.radius)
		if self.packing:
			optCavShape = optCavShape.difference(optCavShape.buffer(-self.packingWidth))
		self.outsiteLength = optCavShape.exterior.length
		self.shape = optCavShape
		optCav = getHollowShapePoints(optCavShape, layer=DRIE1_layer)
		optCav = utils.translate(optCav,self.position)
		cell.add(optCav)
		cellBox = shapely.geometry.box(-self.refCell.dicingWidth/2.,-self.refCell.cellHeight-self.refCell.dicingWidth/2.,self.refCell.cellWidth+self.refCell.dicingWidth/2.,self.refCell.dicingWidth/2.)
		if self.otsMaskOffset>0:
			otsMaskShape = Point(self.position).buffer(self.radius+self.otsMaskOffset)
			cellBox = cellBox.difference(otsMaskShape)
		otsMask = getHollowShapePoints(cellBox, layer=OTS_MASK_layer)
		cell.add(otsMask)
		return 	cell, self.toElement()

	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("radius",str(self.radius))
		root.set("perimeter",str(self.outsiteLength))
		return root		

class ThermIsolCavs(OptCav):
	name="ThermIsolCavs"
	dispCavSize = 1600
	wallThickness = 300
	trenchThichness = 300
	optCavCenter = np.array([0,0])
	dispCavCenter = np.array([0,-2000])
	chanAreaWidth = 500
	bridgeWidth = 60
	# topBridgeAngles = [-1/2.,0,1/2.]
	topBridgeAngles = np.array([])
	botBridgeAngles = np.array([-4/6.,-2/6.,0,4/6.,2/6.])
	isolationCavs = True
	
	def make(self):
		cavities = core.Elements()
		
		optCavShape = Point(self.optCavCenter).buffer(self.radius,resolution=50)
		optCavMask = Polygon(optCavShape.exterior)
		if self.packing:
			optCavShape = optCavShape.difference(optCavShape.buffer(-self.packingWidth))
		self.outsiteLength = optCavShape.exterior.length
		cavities.add(getHollowShapePoints(optCavShape, layer=DRIE1_layer))
		
		dispCavShape = Point(self.dispCavCenter).buffer(self.dispCavSize/2,resolution=50)#.envelope
		dispCavMask = Polygon(dispCavShape.exterior)
		if self.packing:
			dispCavShape = dispCavShape.difference(dispCavShape.buffer(-self.packingWidth))
		cavities.add(getHollowShapePoints(dispCavShape, layer=DRIE1_layer))

		#Channels
		chanSpacing = 400
		chanWidth = 5
		chanLine = LineString([self.optCavCenter,self.dispCavCenter])
		channels=[chanLine]
		channels.append(chanLine.parallel_offset(chanSpacing, 'right', join_style=1, mitre_limit=chanSpacing))
		channels.append(chanLine.parallel_offset(chanSpacing, 'left', join_style=1, mitre_limit=chanSpacing))
		chanClip=optCavMask.union(dispCavMask)
		for line in channels:
			chanPoly = Polygon(line.buffer(chanWidth,cap_style=2).exterior).difference(chanClip)
			cavities.add(getHollowShapePoints(chanPoly, layer=DRIE1_layer))
		
		chanAreaMask = LineString([self.optCavCenter,self.dispCavCenter]).buffer(self.chanAreaWidth,cap_style=2)
		clippingArea = optCavMask.buffer(self.wallThickness)
		clippingArea = clippingArea.union(dispCavMask.buffer(self.wallThickness))
		clippingArea = clippingArea.union(chanAreaMask.buffer(self.wallThickness,cap_style=2))
		#Provide shape for making conform serpentine later
		self.refCell.heaterShape = clippingArea	
	
		clippingArea = Polygon(clippingArea.buffer(self.trenchThichness,join_style=2,mitre_limit=1000).exterior).difference(clippingArea)

		#Bridges
		maxBridgeLength = 2000
		for angle in self.botBridgeAngles:
			vector = np.array([maxBridgeLength*math.cos(angle*math.pi-math.pi/2),maxBridgeLength*math.sin(angle*math.pi-math.pi/2)])
			bridgeLine = LineString([self.dispCavCenter,self.dispCavCenter+vector]).buffer(self.bridgeWidth/2)
			clippingArea = clippingArea.difference(bridgeLine)
		for angle in self.topBridgeAngles:
			vector = np.array([maxBridgeLength*math.cos(angle*math.pi+math.pi/2),maxBridgeLength*math.sin(angle*math.pi+math.pi/2)])
			bridgeLine = LineString([self.optCavCenter,self.optCavCenter+vector]).buffer(self.bridgeWidth/2)
			clippingArea = clippingArea.difference(bridgeLine)

		if self.isolationCavs is True:
			cavities.add(getHollowShapePoints(clippingArea, layer=DRIE1_layer))
		return utils.translate(cavities,self.position), self.toElement()

	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("radius",str(self.radius))
		root.set("perimeter",str(self.outsiteLength))
		return root		

class SpecialExtendedCavs(OptCav):
	name="SpecialExtendedCavs"
	dispCavSize = 1600
	wallThickness = 150
	trenchThichness = 300
	boundaryThickness = 400
	optCavCenter = np.array([0,0])
	dispCavCenter = np.array([0,-2000])
	extensionPercentage = 0.0
	optCavArea=0
	optCavPerimeter=0
	dispCavArea=0
	dispCavPerimeter=0
	extensionCavArea=0
	extensionCavPerimeter=0
	
	def make(self):
		cavities = core.Elements()
		
		optCavShape = Point(self.optCavCenter).buffer(self.radius,resolution=50)
		optCavMask = Polygon(optCavShape.exterior)
		self.optCavArea = optCavShape.area
		self.optCavPerimeter = optCavShape.length
		if self.packing:
			optCavShape = optCavShape.difference(optCavShape.buffer(-self.packingWidth))
		self.outsiteLength = optCavShape.exterior.length
		cavities.add(getHollowShapePoints(optCavShape, layer=DRIE1_layer))
		
		dispCavShape = Point(self.dispCavCenter).buffer(self.dispCavSize/2,resolution=50)#.envelope
		dispCavMask = Polygon(dispCavShape.exterior)
		self.dispCavArea=dispCavShape.area
		self.dispCavPerimeter=dispCavShape.length
		if self.packing:
			dispCavShape = dispCavShape.difference(dispCavShape.buffer(-self.packingWidth))
		cavities.add(getHollowShapePoints(dispCavShape, layer=DRIE1_layer))
		
		chanSpacing = 400
		chanWidth = 5
		chanLine = LineString([self.optCavCenter,self.dispCavCenter])
		channels=[chanLine]
		channels.append(chanLine.parallel_offset(chanSpacing, 'right', join_style=1, mitre_limit=chanSpacing))
		channels.append(chanLine.parallel_offset(chanSpacing, 'left', join_style=1, mitre_limit=chanSpacing))
		chanClip=optCavMask.union(dispCavMask)
		for line in channels:
			chanPoly = Polygon(line.buffer(chanWidth,cap_style=2).exterior).difference(chanClip)
			cavities.add(getHollowShapePoints(chanPoly, layer=DRIE1_layer))
			
		chanAreaMask = LineString([self.optCavCenter,self.dispCavCenter]).buffer(500,cap_style=2)
		clippingArea = optCavMask.union(dispCavMask).union(chanAreaMask).exterior
		
		clippingArea = Polygon(clippingArea.union(Polygon(optCavShape.exterior)).buffer(self.wallThickness,join_style=2,mitre_limit=1000).exterior)
		
		#Provide shape for making conform serpentine later
		self.refCell.heaterShape = clippingArea	
		
		#Cell bounding box
		cellBox = shapely.geometry.box(0,-self.refCell.cellHeight,self.refCell.cellWidth,0)
		cellBox = shapely.affinity.translate(cellBox, -self.position[0],-self.position[1])
		pie = shapely.affinity.translate(uniteReflection(pieShape(3000,self.extensionPercentage*math.pi,20)),self.dispCavCenter[0],self.dispCavCenter[1])		
		#Keep only the bottom part up to the middle of the optical cavity
		extCavShape = shapely.geometry.box(-self.refCell.cellWidth/2,-self.refCell.cellHeight,self.refCell.cellWidth/2,0).intersection(cellBox)
		extCavShape = extCavShape.difference(pie)
		extCavShape = Polygon(cellBox.buffer(-self.boundaryThickness).exterior).difference(clippingArea).intersection(extCavShape)
		#Geodesic closing  (remove features whose size is smaller than 300)
		extCavShape = extCavShape.buffer(-300,join_style=1).buffer(300,join_style=1)
		self.extensionCavArea=extCavShape.area
		self.extensionCavPerimeter=extCavShape.length
		#Make packing donut shape
		extCavShape = extCavShape.difference(extCavShape.buffer(-300,join_style=1))

		if extCavShape.area > 0:
			cavities.add(getHollowShapePoints(extCavShape, layer=DRIE1_layer))
			
		return utils.translate(cavities,self.position), self.toElement()

	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("radius",str(self.radius))
		root.set("perimeter",str(self.outsiteLength))
		root.set("optCavArea",str(self.optCavArea))
		root.set("optCavPerimeter",str(self.optCavPerimeter))
		root.set("dispCavArea",str(self.dispCavArea))
		root.set("dispCavPerimeter",str(self.dispCavPerimeter))
		root.set("extensionCavArea",str(self.extensionCavArea))
		root.set("extensionCavPerimeter",str(self.extensionCavPerimeter))
		root.set("extensionPercentage",str(self.extensionPercentage))
		return root		

class VolVarCavs(OptCav):
	name="SpecialExtendedCavs2"
	radius = 900
	dispCavSize = 1600
	boundaryThickness = 350
	optCavCenter = np.array([0,-2000])
	dispCavCenter = np.array([0,-4600])
	optCavThick = 100
	dispCavThick = 200
	chanWidth = 5
	bottomCropLimit = 3000
	finNumber = 5
	botChan = False
	#Metrics
	optCavArea=0
	optCavPerimeter=0
	dispCavArea=0
	dispCavPerimeter=0
	extensionCavArea=0
	extensionCavPerimeter=0
	
	def make(self):
		self.position = np.array([self.refCell.cellWidth/2,0])
		cavities = core.Elements()
		
		#Cell bounding box
		cellBox = shapely.geometry.box(-self.refCell.cellWidth/2,-self.refCell.cellHeight,self.refCell.cellWidth/2,0)
		cropArea = shapely.geometry.box(-self.refCell.cellWidth/2,-self.refCell.cellHeight+self.bottomCropLimit,self.refCell.cellWidth/2,0)
		innerCellBox = cellBox.buffer(-self.boundaryThickness)

		#Disp cav
		dispCavShape = Point(self.dispCavCenter).buffer(self.dispCavSize/2,resolution=50)#.envelope
		self.dispCavArea=dispCavShape.area
		self.dispCavPerimeter=dispCavShape.length
		
		dispCavOuterWall = dispCavShape.buffer(self.dispCavThick,resolution=50)
		if self.packing:
			dispCavShape = dispCavShape.difference(dispCavShape.buffer(-self.packingWidth))
		cavities.add(getHollowShapePoints(dispCavShape, layer=DRIE1_layer))
		#Disp cav Bridge
		dispCavBridge = LineString([self.dispCavCenter,np.array([0,-self.refCell.cellHeight])]).buffer(250,cap_style=2)
		
		
		#Opt cav
		optCavShape = Point(self.optCavCenter).buffer(self.radius,resolution=50)
		optCavShape = optCavShape.difference(dispCavOuterWall).difference(dispCavBridge).intersection(innerCellBox).buffer(-149).buffer(149)
		self.optCavArea = optCavShape.area
		self.optCavPerimeter = optCavShape.length
		
		optCavOuterWall = optCavShape.buffer(self.optCavThick,resolution=50)
		if self.packing:
			optCavShape = optCavShape.difference(optCavShape.buffer(-self.packingWidth))
		cavities.add(getHollowShapePoints(optCavShape, layer=DRIE1_layer))
		#Opt cav Bridge
		optCavBridge = LineString([self.optCavCenter,np.array([0,0])]).buffer(250,cap_style=2)
		

		#Inter cav Bridge
		interCavBridge = LineString([self.dispCavCenter,self.optCavCenter]).buffer(200,cap_style=2)
		
		#Clipping area
		clippingArea = optCavOuterWall.union(optCavBridge).union(dispCavOuterWall).union(dispCavBridge).union(interCavBridge)
		
		#Fins
		finHalfThickness = 25
		firstLineYoffset = -self.boundaryThickness+finHalfThickness
		finLine =  LineString([np.array([-self.refCell.cellWidth/2,firstLineYoffset]),np.array([self.refCell.cellWidth/2,firstLineYoffset])])
		fins = finLine.buffer(finHalfThickness)
		for i in range(0,self.finNumber):
			fins =	fins.union(finLine.parallel_offset((300+finHalfThickness*2)*i, 'right', join_style=1, mitre_limit=400).buffer(finHalfThickness,cap_style=2))
		
		#Stripes
		stripes = innerCellBox.intersection(cropArea).difference(fins)
		
		#Ext cav
		extCavShape = stripes.difference(clippingArea).buffer(-149).buffer(149)
		self.extensionCavArea = extCavShape.area
		self.extensionCavPerimeter = extCavShape.length
		
		#Channels
		vertLinesXpos = 1400
		lineSet = MultiLineString()
		if extCavShape.area>0:
			linesYmin = extCavShape.bounds[1]+50
			linesYmax = extCavShape.bounds[3]-50
			leftLine = LineString([np.array([-vertLinesXpos,linesYmin]),np.array([-vertLinesXpos,linesYmax])])
			rightLine = LineString([np.array([vertLinesXpos,linesYmin]),np.array([vertLinesXpos,linesYmax])])
			topLine = LineString([np.array([-vertLinesXpos,linesYmax]),np.array([vertLinesXpos,linesYmax])])
			botLine = LineString([np.array([-vertLinesXpos,linesYmin]),np.array([vertLinesXpos,linesYmin])])
			topChanLine = LineString([np.array([0,linesYmax]),self.optCavCenter])
			botChanLine = LineString([np.array([0,linesYmin]),self.dispCavCenter])
			if self.botChan:
				lineSet = MultiLineString([leftLine,rightLine,topLine,botLine,topChanLine,botChanLine])
			else:
				lineSet = MultiLineString([leftLine,rightLine,topLine,topChanLine])
		if not self.botChan:
			midChanLine = LineString([self.optCavCenter,self.dispCavCenter])
			lineSet = lineSet.union(midChanLine)
		if not lineSet.is_empty:
			line = lineSet.buffer(self.chanWidth).difference(Polygon(dispCavShape.exterior).union(Polygon(optCavShape.exterior)))
			line = line.difference(extCavShape)
			cavities.add(getHollowShapePoints(line, layer=DRIE1_layer))
		packingArea = extCavShape.buffer(-300)
		if packingArea.area > 0:
			extCavShape = extCavShape.difference(packingArea)
		if extCavShape.area > 0:
			cavities.add(getHollowShapePoints(extCavShape, layer=DRIE1_layer))
			
		return utils.translate(cavities,self.position), self.toElement()

	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("radius",str(self.radius))
		# root.set("perimeter",str(self.outsiteLength))
		root.set("optCavArea",str(self.optCavArea))
		root.set("optCavPerimeter",str(self.optCavPerimeter))
		root.set("dispCavArea",str(self.dispCavArea))
		root.set("dispCavPerimeter",str(self.dispCavPerimeter))
		root.set("extensionCavArea",str(self.extensionCavArea))
		root.set("extensionCavPerimeter",str(self.extensionCavPerimeter))
		root.set("totalArea",str(self.extensionCavArea+self.dispCavArea+self.optCavArea))
		root.set("totalPerimeter",str(self.extensionCavPerimeter+self.dispCavPerimeter+self.optCavPerimeter))
		return root		

class BreakSealCavs(OptCav):
	name="BreakSealCavs"
	radius = 900
	dispCavSize = 1600
	boundaryThickness = 350
	optCavCenter = np.array([0,-2000])
	dispCavCenter = np.array([0,-4600])
	optCavThick = 100
	sealThickness = 250
	chanWidth = 5
	optCavCropLimit = 0
	northWidth = 300
	northHeight = 300
	northSeparation = 100
	southWidth = 300
	southHeight = 300
	breakSealThickness = 50
	southSeparation = 0
	dispCavYoffset = 0
	dispCavOpening = False
	#Metrics
	optCavArea=0
	dispCavArea=0
	
	def make(self):
		self.optCavArea=0
		self.dispCavArea=0
		self.position = np.array([self.refCell.cellWidth/2,0])
		cavities = core.Elements()
		clippingArea = MultiPolygon()
		
		#Cell bounding box
		cellBox = shapely.geometry.box(-self.refCell.cellWidth/2,-self.refCell.cellHeight,self.refCell.cellWidth/2,0)
		innerCellBox = cellBox.buffer(-self.boundaryThickness)
		
		#Opt cav
		optCavShape = Point(self.optCavCenter).buffer(self.radius,resolution=50)	
		optCavOuterWall = optCavShape.buffer(self.optCavThick,resolution=50)
		clippingArea = clippingArea.union(optCavOuterWall)
		self.optCavArea += optCavShape.area
		if self.packing:
			optCavShape = optCavShape.difference(optCavShape.buffer(-self.packingWidth))
		cavities.add(getHollowShapePoints(optCavShape, layer=DRIE1_layer))
		
		#Opt cav Bridge
		optCavBridge = LineString([self.optCavCenter,np.array([0,0])]).buffer(self.sealThickness,cap_style=2)
		clippingArea = clippingArea.union(optCavBridge)
		
		#North cav
		northCav = shapely.geometry.box(-self.northWidth/2,optCavShape.bounds[1]-self.northHeight-self.northSeparation,self.northWidth/2,optCavShape.bounds[1]-self.northSeparation)
		self.optCavArea += northCav.area
		cavities.add(getHollowShapePoints(northCav, layer=DRIE1_layer))
		
		
		#North Chan
		northChanLine = LineString([p2npa(northCav.centroid),self.optCavCenter])
		northChan = northChanLine.buffer(self.chanWidth)
		northChan = northChan.difference(Polygon(optCavShape.exterior))
		northChan = northChan.difference(northCav)
		cavities.add(getHollowShapePoints(northChan, layer=DRIE1_layer))	

		#South cav
		southCav = shapely.geometry.box(-self.southWidth/2,northCav.bounds[1]-self.southHeight-self.breakSealThickness,self.southWidth/2,northCav.bounds[1]-self.breakSealThickness)
		
		#Disp cav
		dispExlusionArea = Polygon(optCavShape.exterior).union(northChan).buffer(self.sealThickness).union(northCav.buffer(self.breakSealThickness)).union(optCavBridge)
		dispCavShape = shapely.geometry.box(-self.dispCavSize/2,southCav.bounds[1]-self.southSeparation-self.dispCavSize,self.dispCavSize/2,southCav.bounds[1]-self.southSeparation)#.envelope
		dispCavShape = shapely.affinity.translate(dispCavShape, yoff=self.dispCavYoffset)
		dispCavShape = dispCavShape.difference(dispExlusionArea)
		dispCavShape = dispCavShape.intersection(innerCellBox)
		if self.dispCavOpening: dispCavShape = dispCavShape.buffer(-149).buffer(149)
		self.dispCavArea += dispCavShape.union(southCav).area
		cavities.add(getHollowShapePoints(dispCavShape.union(southCav), layer=DRIE2_layer))
		
		#Opt cav extension
		optExclusionArea = dispCavShape.buffer(self.sealThickness).union(southCav.union(northCav).buffer(self.breakSealThickness))
		optCropArea = shapely.geometry.box(-self.refCell.cellWidth/2,0,self.refCell.cellWidth/2,-self.optCavCropLimit)
		optCavExt=innerCellBox.intersection(optCropArea).difference(clippingArea).difference(optExclusionArea).buffer(-149).buffer(149)
		self.optCavArea += optCavExt.area
		packingArea = optCavExt.buffer(-300)
		if packingArea.area > 0:
			optCavExt = optCavExt.difference(packingArea)
		if optCavExt.area > 0:
			cavities.add(getHollowShapePoints(optCavExt, layer=DRIE1_layer))
			
		#Channels
		vertLinesXpos = 1400
		lineSet = MultiLineString()
		if optCavExt.area>0:
			linesYmax = optCavExt.bounds[3]-150
			topLine = LineString([np.array([-vertLinesXpos,linesYmax]),np.array([vertLinesXpos,linesYmax])])
			topChanLine = LineString([np.array([0,linesYmax]),self.optCavCenter])
			lineSet = MultiLineString([topLine,topChanLine])
			line = lineSet.buffer(self.chanWidth)
			line = line.difference(Polygon(optCavShape.exterior))
			line = line.difference(optCavExt)
			cavities.add(getHollowShapePoints(line, layer=DRIE1_layer))	
		
		return utils.translate(cavities,self.position), self.toElement()

	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("radius",str(self.radius))
		root.set("optCavArea",str(self.optCavArea))
		root.set("dispCavArea",str(self.dispCavArea))
		root.set("areaRatioOptDisp",str(self.optCavArea/self.dispCavArea))
		return root		
		
class StarRing(OptCav):
	name="StarRingOpticalCavity"
	frac1=28
	frac2=1
	frac3=1000
	DR=300
	
	def make(self):
		angle1=math.pi
		points=[]
		for i in range(0,self.frac1):
			offset =math.pi+2*math.pi*i/self.frac1
			points.append(cartCoord(self.radius,-angle1/self.frac1+offset))
			points.append(cartCoord(self.radius,-angle1/(self.frac1*self.frac2)+offset))
			points.append(cartCoord(self.radius+self.DR,-angle1/(self.frac1*self.frac3)+offset))
			points.append(cartCoord(self.radius+self.DR,+angle1/(self.frac1*self.frac3)+offset))
			points.append(cartCoord(self.radius,+angle1/(self.frac1*self.frac2)+offset))

		polygon = Polygon(points).simplify(1,preserve_topology=True)
		self.outsiteLength = polygon.exterior.length
		if self.packing:
			polygon = polygon.difference(polygon.buffer(-self.packingWidth,join_style=1,mitre_limit=100))
			# polygon = Polygon(getHollowShapePoints(polygon)).buffer(0)
		self.shape = polygon
		optCav = getHollowShapePoints(polygon, layer=DRIE1_layer)
		return utils.translate(optCav,self.position), self.toElement()
		
	def toElement(self):
		root = super(self.__class__, self).toElement()
		root.set("DR",str(self.DR))
		root.set("Frac1",str(self.frac1))
		return root
		
class Trace(object):
	name="Trace"
	thickness = 100e-9
	material = Materials.platinum
	refCell=None
	wireLength = 0
	wireWidth = 10	  # wire width
	position = np.array([500,500])
	lineSpacing = 10
	connectorPoint = np.array([0,0])
	
	def __init__(self, refCell):
		self.refCell = refCell
		
	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("wireLength",str(self.wireLength))
		root.set("wireWidth",str(self.wireWidth))
		root.set("lineSpacing",str(self.lineSpacing))
		root.set("material",str(self.material.name))
		root.set("resistance",str(self.computeResistance()))
		root.set("thickness",str(self.thickness/1e-6))
		return root
	
	def computeResistance(self):
		return self.material.resistivity*self.wireLength/(self.thickness*self.wireWidth)

class Connector(Trace):
	destTraceObj = None
	padPointPos = None
	padSeparation = None

	def setParam(self, destTraceObj, padPointPos,padSeparation):
		self.destTraceObj = destTraceObj
		self.padPointPos = padPointPos
		self.padSeparation = padSeparation

	def make(self):
		sConnector = core.Elements()
		pointA = self.padPointPos
		pointB=self.destTraceObj.connectorPoint+self.destTraceObj.position
		lineSpacing=self.destTraceObj.lineSpacing
		wireWidth=self.destTraceObj.wireWidth
		midPoint1 = np.array([pointA[0],(pointA[1]+pointB[1])/2])
		midPoint2 = np.array([pointB[0],(pointA[1]+pointB[1])/2])
		points = np.array([pointA+np.array([0,-wireWidth/2]),midPoint1,midPoint2,pointB])
		line = LineString(points)
		
		rightLine = line.parallel_offset(lineSpacing, 'right', join_style=1, mitre_limit=lineSpacing)
		rightSegmentPt1 = pointA+np.array([0,-wireWidth/2])-np.array([-self.padSeparation/2,0])	
		rightLine = np.vstack((rightLine,rightSegmentPt1))	
		rightTrace=core.Path(rightLine,wireWidth, layer = SPUT2_layer)
		sConnector.add(rightTrace)
		
		leftLine = line.parallel_offset(lineSpacing, 'left', join_style=1, mitre_limit=lineSpacing)
		leftSegmentPt1 = pointA+np.array([0,-wireWidth/2])-np.array([+self.padSeparation/2,0])
		leftLine = np.vstack((leftSegmentPt1,leftLine))	
		leftTrace=core.Path(leftLine,wireWidth, layer = SPUT2_layer)
		sConnector.add(leftTrace)
		
		return sConnector

class SpecialConnector(Connector):
	align="right"
	wireWidth = 50
	lineSpacing = 45
	offsetAB = 0
	offset = 100
	reductorLength = 100
	overlapLength = 80
	
	def make(self):
		connector = core.Elements()
		#The portion with small thickness going down from the trace to the expander (A,B)
		pointA = self.destTraceObj.connectorPoint+self.destTraceObj.position
		pointB= pointA + np.array([0,-self.offsetAB])
		def makeDoubleLine(points,lineSpacing,wireWidth,layer):
			elements = core.Elements()
			line = LineString(points)
			if line.length>0:
				leftLine = line.parallel_offset(lineSpacing, 'left', join_style=1, mitre_limit=lineSpacing)
				rightLine = line.parallel_offset(lineSpacing, 'right', join_style=1, mitre_limit=lineSpacing)
				elements.add(core.Path(leftLine,wireWidth, layer))
				elements.add(core.Path(rightLine,wireWidth, layer))
			return elements
		connector.add(makeDoubleLine([pointA,pointB],self.destTraceObj.lineSpacing,self.destTraceObj.wireWidth,SPUT1_layer))
		#The large lines going from the pads to the expander (padPointPos,PointD,pointC)
		if self.align=="right": factor=1
		elif self.align=="left": factor=-1
		else: factor=0
		self.center = (self.destTraceObj.lineSpacing+self.destTraceObj.wireWidth/2-self.lineSpacing-self.wireWidth/2)*factor
		points=[]
		pointC = pointB + np.array([self.center,-self.reductorLength])
		pointD = pointC + np.array([0,-self.overlapLength])
		pointE = pointD + np.array([0,-self.offset])
		pointF = self.padPointPos+np.array([0,self.offset]) #F
		pointG = self.padPointPos #G
		points.append(pointG)
		points.append(pointF)
		points.append(pointE)
		points.append(pointD)
		points.append(pointC)
		connector.add(makeDoubleLine(points,self.lineSpacing,self.wireWidth,SPUT2_layer))

		#The expander (B,pointC)
		def getLineBorderPoints(point,lineSpacing,wireWidth):
			points=[]
			points.append(point+np.array([-lineSpacing-wireWidth/2,0]))
			points.append(point+np.array([-lineSpacing+wireWidth/2,0]))
			points.append(point+np.array([lineSpacing-wireWidth/2,0]))
			points.append(point+np.array([lineSpacing+wireWidth/2,0]))
			return points
		end1Pts = getLineBorderPoints(pointB,self.destTraceObj.lineSpacing,self.destTraceObj.wireWidth)
		end2Pts = getLineBorderPoints(pointC,self.lineSpacing,self.wireWidth)
		poly1, poly2 = Polygon([end1Pts[0],end1Pts[1],end2Pts[1],end2Pts[0]]),Polygon([end1Pts[2],end1Pts[3],end2Pts[3],end2Pts[2]])
		connector.add(core.Boundary(poly1.exterior.coords, layer = SPUT1_layer))
		connector.add(core.Boundary(poly2.exterior.coords, layer = SPUT1_layer))
		
		#Add the overlapping area
		connector.add(makeDoubleLine([pointC,pointD],self.lineSpacing,self.wireWidth,SPUT1_layer))
		
		return connector
		
class Serpentine(Trace):
	name = "Serpentine Trace"
	height = 5000  # device width
	width = 3000
	windNumber = 3		  # number of windings
	lineSpacing = 10
	
	def make(self):
		self.position = np.array([self.refCell.cellWidth/2-self.width/2,-self.refCell.cellHeight/2-self.height/2])
		serpentine = core.Elements()
		spacing = self.width/self.windNumber # spacing between windings
		endOffset = self.lineSpacing
		# startOffset = self.lineSpacing
		startOffset = 0
		unit = np.array([[0,0], [0, self.height], [spacing/2., self.height], [spacing/2., 0]])
		pts=unit
		for i in range(1,self.windNumber):
			next_unit = unit + i * np.array([spacing, 0])
			pts = np.vstack((pts, next_unit))
		pts=np.vstack(([0,-startOffset], pts, [spacing * self.windNumber, 0], [spacing * self.windNumber, self.height+endOffset]))
		line = LineString(pts)
		offset = line.parallel_offset(self.lineSpacing, 'right', join_style=1, mitre_limit=self.lineSpacing)
		offset2 = line.parallel_offset(self.lineSpacing, 'left', join_style=1, mitre_limit=self.lineSpacing)
		pts = np.vstack((offset2,offset))
		line = LineString(pts)
		self.wireLength = line.length
		trace=core.Path(pts, self.wireWidth, layer = SPUT1_layer)
		serpentine.add(trace)
		print >>f, ET.tostring(self.toElement())
		return utils.translate(serpentine,self.position)

class SquareSpiral(Trace):
	name = "SquareSpiral"
	height = 5000  # device width
	width = 3000
	windNumber = 4		  # number of windings
	windSeparation = 200
		
	def make(self):
		self.position = np.array([self.refCell.cellWidth/2-self.width/2,-self.refCell.cellHeight/2-self.height/2])
		spiral = core.Elements()
		points=np.array([0,0])
		for i in range(1,self.windNumber):
			offset = (i-1)*self.windSeparation
			point1=np.array([0+offset,self.height-offset])
			point2=np.array([self.width-offset,self.height-offset])
			point3=np.array([self.width-offset,0+offset])
			point4=np.array([0+offset+self.windSeparation,0+offset])
			points = np.vstack((points,point1,point2,point3,point4))
		line = LineString(points)	
		offset = line.parallel_offset(self.lineSpacing, 'right', join_style=1, mitre_limit=self.lineSpacing)
		offset2 = line.parallel_offset(self.lineSpacing, 'left', join_style=1, mitre_limit=self.lineSpacing)		
		points = np.vstack((offset2,offset))
		line = LineString(points)
		self.wireLength = line.length
		trace=core.Path(points, self.wireWidth, layer = SPUT1_layer)
		spiral.add(trace)
		return utils.translate(spiral,self.position), self.toElement()
	
class RoundSerpentine(Trace):
	name = "RoundSerpentine"
	width = 3200
	windNumber = 3		  # number of windings
	windSeparation = 300
	centralSep = 100
	separationXoffset = 0
	mirror = False
	autoPosition = True
		
	def make(self):
		# self.position = np.array([self.refCell.cellWidth/2-self.width/2,self.width/2-self.refCell.optCavYoffset])
		if self.autoPosition: self.position = np.array([self.refCell.cellWidth/2,-self.refCell.optCavYoffset])
		spiral = core.Elements()
		#Def points
		center = np.array([0,0])
		left = np.array([-self.width/2,0])
		bottom = np.array([0,-self.width/2])
		mask = LineString([bottom,center]).buffer(self.centralSep).envelope
		mask = shapely.affinity.translate(mask,self.separationXoffset,0)
		self.connectorPoint = np.array(bottom-[self.centralSep-self.separationXoffset,self.centralSep])
		points = [self.connectorPoint]
		for i in range(0,self.windNumber):
			circle = Point(center).buffer(self.width/2-i*self.windSeparation,resolution=36).buffer(0)
			line = LineString(circle.exterior.coords)
			line = line.difference(mask)
			arc = np.vstack((line[1].coords,line[0].coords))
			if i%2==0:
				points = np.vstack((points,arc))
			else:
				points = np.vstack((points,arc[::-1]))
		
		line = LineString(points)
		# Double the line
		offset = line.parallel_offset(self.lineSpacing, 'right', join_style=2, mitre_limit=self.lineSpacing)
		offset2 = line.parallel_offset(self.lineSpacing, 'left', join_style=2, mitre_limit=self.lineSpacing)		
		points = np.vstack((offset2,offset))
		line2 = LineString(points)
		if self.mirror:
			line2 = shapely.affinity.scale(line2,xfact=-1.0)
			self.connectorPoint = self.connectorPoint * np.array([-1,1])
		self.wireLength = line2.length
		spiral.add(getHollowShapePoints(line2.buffer(self.wireWidth/2,cap_style=2,join_style=2, mitre_limit=self.lineSpacing).simplify(0.1),layer=SPUT1_layer))
		#spiral.add(core.Path(line2, self.wireWidth, layer = SPUT1_layer))
		
		# print >>f, ET.tostring()
		return utils.translate(spiral,self.position), self.toElement()
		
class ConformSerpentine(RoundSerpentine):
	name = "ConformSepentine"
	windSeparation = 40
	windNumber = 1
	position = np.array([0,0])
	radialOffset = 50
	outConnLength = 150
	separationXoffset = 0

	def make(self):
		spiral = core.Elements()
		shape = self.refCell.heaterShape
		if shape is None: print "Warning: ConformSerpentine needs a shape to be defined"
		center = shape.centroid+np.array([0,0])
		bottom = center+np.array([0,-6000])
		mask = LineString([bottom,center]).buffer(self.centralSep).envelope
		mask = shapely.affinity.translate(mask,self.separationXoffset,0)
		localConnectorPoint = np.array([shape.centroid.x-self.centralSep+self.separationXoffset,shape.bounds[1]-self.outConnLength])
		points = [localConnectorPoint]
		self.connectorPoint = localConnectorPoint + self.refCell.optCavObj.position
		
		for i in range(0,self.windNumber):
			circle = Polygon(shape.buffer(self.radialOffset-i*self.windSeparation,resolution=50,join_style=1,mitre_limit=1000).exterior)
			line = LineString(circle.exterior)
			line = line.difference(mask)
			assert len(line)==2, "The line should be cut in 2 parts"
			arc = np.vstack((line[1].coords,line[0].coords))
			if i%2==0:
				points = np.vstack((points,arc))
			else:
				points = np.vstack((points,arc[::-1]))
		line = LineString(points)
		# Double the line
		offset = line.parallel_offset(self.lineSpacing, 'right', join_style=2, mitre_limit=1000).simplify(0.1, preserve_topology=False)
		offset2 = line.parallel_offset(self.lineSpacing, 'left', join_style=2, mitre_limit=1000).simplify(0.1, preserve_topology=False)
		points = np.vstack((offset2,offset))
		line2 = LineString(points)
		if self.mirror:
			line2 = shapely.affinity.scale(line2,xfact=-1.0)
			self.connectorPoint = localConnectorPoint* np.array([-1,1]) + self.refCell.optCavObj.position
		self.wireLength = line2.length
		# spiral.add(core.Path(line2, self.wireWidth, layer = SPUT1_layer))
		spiral.add(getHollowShapePoints(line2.simplify(0.8).buffer(self.wireWidth/2,cap_style=2,join_style=2, mitre_limit=self.lineSpacing),layer=SPUT1_layer))
		return utils.translate(spiral,self.refCell.optCavObj.position), self.toElement()
		
class Heater:
	name = "Heater"
	refCell = None
	padWidth = 400
	padSeparation = 40
	traceObj = None
	padPairPos = np.array([0,0])
	connectorObj= None

	def __init__(self, refCell):
		self.refCell = refCell
		# self.traceObj = Serpentine(refCell)
		# self.connectorObj = Connector(refCell)

	def make(self):
		bottomLeftPt = np.array([0,-self.refCell.cellHeight])
		heater = core.Elements()
		if self.traceObj:
			trace, traceElement = self.traceObj.make()
			heater.add(trace)
			pad1 = utils.translate(shapes.Rectangle((0,0),(self.padWidth,self.padWidth),layer=4),bottomLeftPt+self.padPairPos)
			pad2 = utils.translate(pad1,(self.padSeparation+self.padWidth,0))
			heater.add(pad1)
			heater.add(pad2)
			padConPoint = np.array([self.padWidth+self.padSeparation/2,self.padWidth])+self.padPairPos
			if self.connectorObj:
				self.connectorObj.destTraceObj = self.traceObj
				self.connectorObj.padPointPos = bottomLeftPt+padConPoint
				self.connectorObj.padSeparation = self.padSeparation
				heater.add(self.connectorObj.make())
		etNode = self.toElement()#.append(traceElement)
		if self.traceObj:
			etNode.append(traceElement)
		return heater, etNode

	def toElement(self):
		root = ET.Element(str(self.name))
		return root

class Sensor(Heater):
	name = "Sensor"
	padPairPos = np.array([2160,0])

class Interconnect(Heater):
	connectWidth = 250
	dicingWidth = 300
	clearWidth = 40
	leftLayer=SPUT2_layer
	rightLayer=SPUT2_layer
	extTop = False
	extBot = False
	extLeft = False
	extRight = False
	intRight = True
	intBottom = True
	topPad = False
	bottomPad = False

	def make(self):
		element = core.Elements()
		bottomLeftPt = np.array([0,-self.refCell.cellHeight])
		bottomRightPt = np.array([self.refCell.cellWidth,-self.refCell.cellHeight])
		leftPadMidLeft = bottomLeftPt + np.array([0,self.padWidth/2])
		rightPadBotMid = bottomLeftPt + np.array([self.padWidth+self.padSeparation+self.padWidth/2,0])
		if self.refCell.heaterObj is not None:
			leftPadMidLeft = bottomLeftPt + self.refCell.heaterObj.padPairPos + np.array([0,self.padWidth/2])
			rightPadBotMid = bottomLeftPt + self.refCell.heaterObj.padPairPos + np.array([self.padWidth+self.padSeparation+self.padWidth/2,0])
		leftConnect = core.Path([(self.connectWidth/2+self.clearWidth,self.dicingWidth+self.clearWidth),(self.connectWidth/2+self.clearWidth,-self.refCell.cellHeight+self.clearWidth)],self.connectWidth,layer=self.leftLayer)
		rightConnect = core.Path([(self.refCell.cellWidth-self.connectWidth/2-self.clearWidth,-self.clearWidth),(self.refCell.cellWidth-self.connectWidth/2-self.clearWidth,-self.refCell.cellHeight-self.dicingWidth-self.clearWidth)],self.connectWidth,layer=self.rightLayer)

		element.add(leftConnect)
		if self.intBottom and self.refCell.heaterObj is not None:
			leftPadConnect = core.Path([np.array([bottomLeftPt[0]+self.clearWidth,leftPadMidLeft[1]]),leftPadMidLeft],self.connectWidth,layer=self.leftLayer)
			bottomPadConnect = core.Path([rightPadBotMid,np.array([rightPadBotMid[0],bottomRightPt[1]-self.dicingWidth-self.clearWidth])],self.connectWidth,layer=self.rightLayer)
			element.add(bottomPadConnect)
			element.add(leftPadConnect)
		if self.extTop:
			extTopConnect = core.Path([(self.clearWidth,self.dicingWidth+self.connectWidth/2+self.clearWidth),(self.refCell.cellWidth+self.dicingWidth+self.clearWidth,self.dicingWidth+self.connectWidth/2+self.clearWidth)],self.connectWidth,layer=self.leftLayer)
			element.add(extTopConnect)
		else:
			if self.intRight: element.add(rightConnect)
		if self.extBot:
			extBotConnect = core.Path([(-self.dicingWidth-self.clearWidth,-self.refCell.cellHeight-self.dicingWidth-self.connectWidth/2-self.clearWidth),(self.refCell.cellWidth-self.clearWidth,-self.refCell.cellHeight-self.dicingWidth-self.connectWidth/2-self.clearWidth)],self.connectWidth,layer=self.rightLayer)
			element.add(extBotConnect)
		else:
			if self.intBottom:
				bottomConnect = core.Path([np.array([rightPadBotMid[0]-self.connectWidth/2,bottomRightPt[1]-self.dicingWidth-self.connectWidth/2-self.clearWidth]),bottomRightPt+np.array([-self.clearWidth,-self.dicingWidth-self.connectWidth/2-self.clearWidth])],self.connectWidth,layer=self.rightLayer)
				element.add(bottomConnect)
		if self.extLeft:
			extLeftConnect = core.Path([(-self.dicingWidth-self.connectWidth/2-self.clearWidth,-self.clearWidth),(-self.dicingWidth-self.connectWidth/2-self.clearWidth,-self.refCell.cellHeight-self.dicingWidth-self.connectWidth-self.clearWidth)],self.connectWidth,layer=self.rightLayer)
			element.add(extLeftConnect)
		if self.extRight:
			extRightConnect = core.Path([(self.refCell.cellWidth+self.dicingWidth+self.connectWidth/2+self.clearWidth,self.dicingWidth+self.connectWidth+self.clearWidth),(self.refCell.cellWidth+self.dicingWidth+self.connectWidth/2+self.clearWidth,-self.refCell.cellHeight+self.clearWidth)],self.connectWidth,layer=self.leftLayer)
			element.add(extRightConnect)
		if self.topPad:
			topPadShape = core.Path([np.array([self.refCell.cellWidth/2,self.dicingWidth+self.clearWidth]),np.array([self.refCell.cellWidth/2,self.refCell.cellWidth+self.dicingWidth+self.clearWidth])],self.refCell.cellWidth,layer=self.leftLayer)
			element.add(topPadShape)
		if self.bottomPad:
			bottomPadShape = core.Path([np.array([self.refCell.cellWidth/2,-self.refCell.cellHeight-self.dicingWidth-self.clearWidth]),np.array([self.refCell.cellWidth/2,-self.refCell.cellHeight-self.refCell.cellWidth-self.dicingWidth-self.clearWidth])],self.refCell.cellWidth,layer=self.leftLayer)
			element.add(bottomPadShape)
		return element

class Channels(object):
	name = "Channels"
	refCell = None
	chanNumber = 2
	chanSpacing = 200
	chanOffset = 500
	chanWidth = 12
	def __init__(self, refCell):
		self.refCell = refCell
		
	def make(self):
		cell = core.Cell(self.name)
		for i in range(0,self.chanNumber):
			assert self.refCell.optCavObj.shape is not None
			transOptCavShape = shapely.affinity.translate(self.refCell.optCavObj.shape,self.refCell.optCavObj.position[0],self.refCell.optCavObj.position[1])
			optCavCenter = p2npa(transOptCavShape.centroid)
			point1 = np.array([self.refCell.cellWidth/2-self.chanOffset-i*self.chanSpacing,-self.refCell.dispCavYoffset])
			delta = np.array([0,self.refCell.cavSeparation/2])
			point2 = point1 + delta
			channelLine = LineString([optCavCenter,point2,point1]).buffer(self.chanWidth/2,cap_style=2)
			channelLine = channelLine.difference(Polygon(transOptCavShape.exterior))
			channel = getHollowShapePoints(channelLine, layer=CHANETCH_layer)
			cell.add(channel)
			cell.add(utils.reflect(channel,'y',origin=(self.refCell.cellWidth/2,0)))
		# cell.add(channels)
		return 	cell, self.toElement()		
	
	def toElement(self):
		root = ET.Element(str(self.name))
		return root

class ChannelsSerp(Channels):
	name = "ChannelsSerp"
	overlap = 100
	dx = np.array([300,0])
	dy = np.array([0,100])
	
	def make(self):
		cell = core.Cell(self.name)
		dx = self.dx
		dy = self.dy
		transOptCavShape = shapely.affinity.translate(self.refCell.optCavObj.shape,self.refCell.optCavObj.position[0],self.refCell.optCavObj.position[1])
		transDispCavShape = shapely.affinity.translate(self.refCell.dispCavObj.shape,self.refCell.dispCavObj.position[0],self.refCell.dispCavObj.position[1])
		optCavCenter = p2npa(transOptCavShape.centroid)
		dispCavCenter = p2npa(transDispCavShape.centroid)
		point0 = np.array([self.refCell.cellWidth/2,-self.refCell.optCavYoffset-self.refCell.optCavObj.radius+self.overlap])
		point1 = np.array([self.refCell.cellWidth/2,-self.refCell.dispCavYoffset-self.overlap])
		pointA = (point0+point1)/2.
		
		channelLine = LineString([point0,pointA+dy,pointA+dx+dy,pointA+dx,pointA-dx,pointA-dx-dy,pointA-dy,point1]).buffer(self.chanWidth/2,cap_style=2)
		# channelLine = channelLine.difference(Polygon(transOptCavShape.exterior))
		channel = getHollowShapePoints(channelLine, layer=CHANETCH_layer)
		cell.add(channel)
		# cell.add(channels)
		return 	cell, self.toElement()		
	
	def toElement(self):
		root = ET.Element(str(self.name))
		return root		

class Cell:
	#Cell parameters
	name = "T-Cell"
	cell = None
	cellWidth = 4000
	cellHeight = 6000
	optCavYoffset = 2000
	dispCavWidth = 1600
	dispCavPacking = True
	cavSeparation = 500
	packingWidth = 300

	dicingCrossLength = 300
	dicingCrossLinewidth = 20
	dicingWidth = 300
	textBool = True
	text2Bool = True
	textSize = 600
	text2Size = 700
	textOffset = 140
	text2Offset = np.array([-1100,-100])
	textLinewidth = 10
	text2Linewidth = 30
	optCavObj = None
	dispCavObj = None
	chanObject = None
	heaterObj = None
	heaterShape = None
	sensorObj= None
	dispCavYoffset = None
	
	def __init__(self):
		self.optCavObj = OptCav(self)
		self.dispCavObj = DispCav(self)
		self.chanObject = ChannelsSerp(self)
	#Make Cell
	def make(self):
		cell = core.Cell(self.name)
		print "Making "+self.name
		#Computed values
		self.dispCavYoffset = self.cavSeparation + self.optCavYoffset + self.optCavObj.radius
		#Make surrounding box
		box = shapes.Rectangle((0,0), (self.cellWidth,-self.cellHeight),layer=LAYOUT_layer)
		cell.add(box)
		#Make optical cavity
		if self.optCavObj:
			self.optCavObj.position = np.array([self.cellWidth/2,-self.optCavYoffset])
			optCav, optCavElement = self.optCavObj.make()
			cell.add(optCav)
		#Make dispenser cavity
		if self.dispCavObj:
			self.dispCavObj.position = np.array([self.cellWidth/2,-self.dispCavYoffset-self.dispCavObj.dispCavSize/2])
			dispCav, dispCavElement = self.dispCavObj.make()
			cell.add(dispCav)
		#Make channels
		channels, channelsElement = self.chanObject.make()
		cell.add(channels)
		
		#Make heater
		if self.heaterObj:
			heater, heaterElement = self.heaterObj.make()
			cell.add(heater)
		if self.sensorObj:
			sensor, sensorElement = self.sensorObj.make()
			cell.add(sensor)
		self.cell = cell
		myShow(cell,self.name)
		#Make Element tree
		etNode = self.toElement()
		if self.optCavObj: etNode.append(optCavElement)
		if self.dispCavObj: etNode.append(dispCavElement)
		if self.heaterObj: etNode.append(heaterElement)
		if self.sensorObj: etNode.append(sensorElement)
		cellListNode.append(etNode)
		return cell, etNode

	def toElement(self):
		root = ET.Element(str(self.name))
		root.set("id",str(id(self)))
		return root

	#Make dicing mark
	def makeDicingCross(self):
		dicingPointBL = (self.dicingCrossLinewidth/2,-self.dicingWidth/2+self.dicingCrossLinewidth/2)
		dicingPointTR = (self.dicingCrossLength-self.dicingCrossLinewidth/2,self.dicingWidth/2-self.dicingCrossLinewidth/2)
		dicingPoints = [dicingPointTR,(dicingPointBL[0],dicingPointTR[1]),dicingPointBL,(dicingPointTR[0],dicingPointBL[1])]
		dicingMark1 = core.Path(dicingPoints, width=self.dicingCrossLinewidth,pathtype=2,layer=CHANETCH_layer)
		dicingMark1 = utils.translate(dicingMark1,(self.dicingWidth/2,0))
		dicingMark2 = utils.rotate(dicingMark1,90)
		dicingMark3 = utils.rotate(dicingMark1,180)
		dicingMark4 = utils.rotate(dicingMark1,270)
		dicingCross = core.Elements((dicingMark1,dicingMark2,dicingMark3,dicingMark4))
		return dicingCross	

	def makeCellDicingMark(self,marks="tblr"):
		dicingCross=self.makeDicingCross()
		dicingCross1 = utils.translate(dicingCross,(-self.dicingWidth/2,self.dicingWidth/2))
		dicingCross2 = utils.translate(dicingCross1,(self.cellWidth+self.dicingWidth,0))
		dicingCross3 = utils.translate(dicingCross1,(0,-self.cellHeight-self.dicingWidth))
		dicingCross4 = utils.translate(dicingCross1,(self.cellWidth+self.dicingWidth,-self.cellHeight-self.dicingWidth))
		dicingMark = core.Elements()
		if 't' in marks or 'l' in marks: dicingMark.add(dicingCross1)
		if 't' in marks or 'r' in marks: dicingMark.add(dicingCross2)
		if 'b' in marks or 'l' in marks: dicingMark.add(dicingCross3)
		if 'b' in marks or 'r' in marks: dicingMark.add(dicingCross4)
		return dicingMark

	def makeText(self,label):
		texts = core.Elements()
		bottomRight = np.array([self.cellWidth,-self.cellHeight])
		if self.textBool : texts.add(shapes.LineLabel(label, self.textSize, (self.textOffset*2,-self.textSize-self.textOffset),line_width=self.textLinewidth,layer=CHANETCH_layer))
		if self.text2Bool : texts.add(shapes.LineLabel(label, self.text2Size,bottomRight+self.text2Offset,line_width=self.text2Linewidth,layer=SPUT2_layer))
		return texts
		
	def adjustSizeFromSealingArea(self,sealingArea):
	#Will override cellWidth, cellHeight and optCavYoffset
		self.cellWidth = self.optCavRad*2 + sealingArea*2
		self.cellHeight = self.optCavRad*2 + sealingArea*2 + self.cavSeparation + self.dispCavWidth
		self.optCavYoffset = sealingArea + self.optCavRad
		return

	def getBoundingBox(self,position=np.array([0,0])):
		return MultiPoint([[0,0],[self.cellWidth,-self.cellHeight]]).envelope
		
class Array:
	#Array parameters
	waferRef = None
	defaultCellObj = None
	specialCells = []
	columnNumber = 19
	rowNumber = 13
	colStartInd = 0
	rowStartInd = 0
	leftDicMarkInd = []
	leftDicMarkExclInd = []
	topDicMarkInd = []
	topDicMarkExclInd = []
	rightDicMarkInd = []
	bottomDicMarkInd = []
	xOffset = 0
	yOffset = 0
	cellNumber = 0
	
	def __init__(self,wafer):
		self.waferRef = wafer

	def computeWH(self):
		arrayWidth = self.columnNumber * self.defaultCellObj.cellWidth + (self.columnNumber-1) * self.defaultCellObj.dicingWidth
		arrayHeight = self.rowNumber * self.defaultCellObj.cellHeight + (self.rowNumber-1) * self.defaultCellObj.dicingWidth
		return arrayWidth,arrayHeight
	
	def make(self):
		#Make cell array with labels
		array = core.Cell('ARRAY')
		arrayWidth,arrayHeight = self.computeWH()
		self.cellNumber=0
		
		posMatrix = np.zeros((self.columnNumber,self.rowNumber,2))
		impCodMatrix = np.zeros((self.rowNumber,self.columnNumber))
		for i in range(0,self.columnNumber):
			for j in range(0,self.rowNumber):
				#Compute cell position
				posX=i*(self.defaultCellObj.cellWidth+self.defaultCellObj.dicingWidth)-arrayWidth/2 + self.xOffset
				posY=-j*(self.defaultCellObj.cellHeight+self.defaultCellObj.dicingWidth)+arrayHeight/2 + self.yOffset
				#Cell is out of the implentation area
				impCode = 0
				ref = core.CellReference(self.defaultCellObj.cell, origin=(posX,posY))
				if shape_in_shape(self.defaultCellObj,np.array([posX,posY]),self.waferRef.implArea):
					if not shapeIntersect(self.defaultCellObj,np.array([posX,posY]),self.waferRef.alignMarks[0]) and not shapeIntersect(self.defaultCellObj,np.array([posX,posY]),self.waferRef.alignMarks[1]):
						#Cell is in implentation area and out of the alignment marks
						impCode = 1
					else:
						#Cell is on the alignement marks
						impCode = 2
				posMatrix[i][j]=(posX,posY)
				impCodMatrix[j][i]=impCode
		leftCumsumMat=np.cumsum(impCodMatrix,axis=0)
		rightCumsumMat=np.flipud(np.cumsum(np.flipud(impCodMatrix),axis=0))
		posList = ET.Element("PositionList")
		for i in range(0,self.columnNumber):
			for j in range(0,self.rowNumber):
				posX=posMatrix[i][j][0]
				posY=posMatrix[i][j][1]
				impCode=impCodMatrix[j][i]
				cellObj=self.defaultCellObj
				#Create a label from column and row numbers
				ref = core.CellReference(cellObj.cell, origin=(posX,posY))
				if impCode>0:
					#Search for special cell
					label=chr(ord('A')+j+self.rowStartInd)+chr(ord('A')+i+self.colStartInd)
					special = filter(lambda x: label in x[0], self.specialCells)
					if len(special)>0:
						cellObj = special[0][1]
						ref = core.CellReference(cellObj.cell, origin=(posX,posY))
					node = ET.Element("Position")
					node.set("label",label)
					node.set("id",str(id(cellObj)))
					posList.append(node)
					text = utils.translate(cellObj.makeText(label),(posX,posY))			
					if impCode==1:
						array.add(text)
						array.add(ref)
						#Add dicing marks
						if (i in self.leftDicMarkInd) and not (j in self.leftDicMarkExclInd):
							array.add(utils.translate(self.defaultCellObj.makeCellDicingMark(marks="l"),(posX,posY)))
						if (i in self.rightDicMarkInd) and not (j in self.leftDicMarkExclInd):
							array.add(utils.translate(self.defaultCellObj.makeCellDicingMark(marks="r"),(posX,posY)))
						if (j in self.topDicMarkInd) and not (i in self.topDicMarkExclInd):
							array.add(utils.translate(self.defaultCellObj.makeCellDicingMark(marks="t"),(posX,posY)))
						if (j in self.bottomDicMarkInd) and not (i in self.topDicMarkExclInd):
							array.add(utils.translate(self.defaultCellObj.makeCellDicingMark(marks="b"),(posX,posY)))
						self.cellNumber = self.cellNumber + 1
						
					#Look for neighbors to make proper connections
					interconnect = Interconnect(cellObj)
					if (j==0 or impCodMatrix[j-1][i]==0):
						if i!=self.columnNumber-1:
							interconnect.extTop = True
							if i==10: interconnect.topPad = True
						else:
							interconnect.intRight = False
					if (j==self.rowNumber-1 or impCodMatrix[j+1][i]==0) and i!=0:
						interconnect.extBot = True
						if i==10: interconnect.bottomPad = True
					if i!=0 and (i==0 or impCodMatrix[j][i-1]==0) and leftCumsumMat[j][i-1]!=0:
						interconnect.extLeft = True
					if i!=self.columnNumber-1 and (i==self.columnNumber-1 or impCodMatrix[j][i+1]==0) and rightCumsumMat[j][i+1]!=0:
						interconnect.extRight = True
					if impCodMatrix[j][i]!=1:
						interconnect.intBottom = False
					array.add(utils.translate(interconnect.make(),(posX,posY)))
					
		arrayNode = ET.Element("Array")
		arrayNode.set("cellNumber",str(self.cellNumber))
		arrayNode.append(posList)
		return array, arrayNode

class AlignMarks:
	inverse = False
	outterWidth = -4
	outterBoxWidth = 350
	crossLength = 150
	squareSize = 0
	crossWidth = 30
	crossSep = 350
	crossDist = 67150/2
	alignBoxesSize = []
	alignBoxesWidth = 10
	arrowsBool = True
	arrowWidth = 10
	arrowDistanceOffset = 0
	arrowSeparation = 100
	arrowNumber = 14
	arrowLength = 30
	arrowAngle = 0
	arrowStartInd = 2
	misProofMarkDist = 70
	misProofMarkSize = 20
	mirrorY = False
	# misProofMark = None

	def makeCross(self,layer,misProofPos=(-1,-1),mirrorX=False,mirrorY=False):
		square = Point(0,0).buffer(self.squareSize/2.).envelope
		line = LineString([(0,-self.crossLength/2),(0,self.crossLength/2)]).buffer(self.crossWidth/2,cap_style=2)
		line = line.union(shapely.affinity.rotate(line, 90, use_radians=False)).difference(square)
		misProofSquare = Point(self.misProofMarkDist*misProofPos[0],self.misProofMarkDist*misProofPos[1]).buffer(self.misProofMarkSize/2.).envelope
		line=line.union(misProofSquare)
		if self.inverse:
			bounds = Point(0,0).buffer(self.outterBoxWidth/2.).envelope
			line = bounds.difference(line.buffer(self.outterWidth,join_style=2,mitre_limit=1000))
		if mirrorX: line = shapely.affinity.scale(line, xfact=-1.0, origin=Point(0,0))
		if mirrorY: line = shapely.affinity.scale(line, yfact=-1.0, origin=Point(0,0))
		return getHollowShapePoints(line, layer)
		
	#Make alignment mark
	def make(self,layer,marks="lr"):
		cross1 = self.makeCross(layer,mirrorX=False,mirrorY=self.mirrorY)
		cross2 = self.makeCross(layer,mirrorX=True,mirrorY=self.mirrorY)
		cross1 = utils.translate(cross1,(self.crossSep/2,0))
		cross2 = utils.translate(cross2,(-self.crossSep/2,0))
		arrow = core.Path([(-self.arrowLength,0),(0,0),(0,self.arrowLength)],self.arrowWidth, layer=layer)
		arrowLine = core.Elements()
		for i in range(self.arrowStartInd, self.arrowNumber):
			arrowLine.add(utils.translate(arrow,(-self.arrowSeparation*i-self.arrowDistanceOffset,self.arrowSeparation*i+self.arrowDistanceOffset)))
		arrows = core.Elements(utils.rotate(arrowLine,0+self.arrowAngle))
		arrows.add(utils.rotate(arrowLine,90+self.arrowAngle))
		arrows.add(utils.rotate(arrowLine,180+self.arrowAngle))
		arrows.add(utils.rotate(arrowLine,270+self.arrowAngle))
		boxes = core.Elements()
		for boxWidth in self.alignBoxesSize:
			box = shapes.Box((-boxWidth,-boxWidth),(boxWidth,boxWidth),width=self.alignBoxesWidth, layer=layer)
			box.pathtype = 2
			boxes.add(box)
		mark1 = core.Elements()
		mark1.add(cross1)
		mark1.add(cross2)
		mark1.add(boxes)
		if self.arrowsBool: mark1.add(arrows)
		mark1 = utils.translate(mark1,(self.crossDist+self.crossSep/2,0))
		mark2 = utils.reflect(mark1,'y')
		mark = core.Elements()
		if "r" in marks:mark.add(mark1)
		if "l" in marks:mark.add(mark2)
		return mark

class Wafer:
	waferRadius = 101600/2
	implanAreaRad = waferRadius-4000
	box1Width = 1500
	box1Height = 20000
	box2Width = 3000
	box2Height = 36000
	box3Width = 2000
	box3Height = 6000
	centeringBoxSize = 500
	centeringBoxDistance = 107000/2
	flatLength = 32500
	implArea = None
	waferArea = None
	alignMarks = []
	featuresLayers = []

	def makeBoxes(self,layer):
		elements = core.Elements()
		def makeBox(width,height,waferRadius,angle,indentationSize,layer):
			box = shapely.geometry.box(-height/2,0,height/2,width)
			box = box.difference(Point(0,0).buffer(indentationSize,resolution=1))
			box = shapely.affinity.translate(box, yoff=self.waferRadius)
			box = shapely.affinity.rotate(box, angle, origin=Point(0,0))
			return getHollowShapePoints(box, layer=layer)

		def makeCenteringBoxes(centeringBoxSize,centeringBoxDistance,layer):
			distance = centeringBoxDistance + centeringBoxSize/2.
			box = shapely.geometry.box(-centeringBoxSize/2,-centeringBoxSize/2,centeringBoxSize/2,centeringBoxSize/2)
			box1 = getHollowShapePoints(shapely.affinity.translate(box, xoff=distance,yoff=distance), layer=layer)
			box2 = getHollowShapePoints(shapely.affinity.translate(box, xoff=-distance,yoff=distance), layer=layer)
			box3 = getHollowShapePoints(shapely.affinity.translate(box, xoff=-distance,yoff=-distance), layer=layer)
			box4 = getHollowShapePoints(shapely.affinity.translate(box, xoff=distance,yoff=-distance), layer=layer)
			return core.Elements((box1,box2,box3,box4))
			
		box1 = makeBox(self.box1Width,self.box1Height,self.waferRadius,-90,40,layer)
		box1b = makeBox(self.box1Width,self.box1Height,self.waferRadius,90,40,layer)
		box2 = makeBox(self.box2Width,self.box2Height,self.waferRadius,180,0,layer)
		box3 = makeBox(self.box3Width,self.box3Height,self.waferRadius,-135,0,layer)
		box3b = makeBox(self.box3Width,self.box3Height,self.waferRadius,135,0,layer)
		
		boxes = makeCenteringBoxes(self.centeringBoxSize,self.centeringBoxDistance,layer)
		elements.add(core.Elements((box1,box1b,box2,box3,box3b)))
		elements.add(boxes)
		return elements
		
	def make(self):
		waferElements = core.Elements()
		layer1marks = AlignMarks()
		self.alignMarks.append(layer1marks.make(DRIE1_layer,'l'))
		self.alignMarks.append(layer1marks.make(DRIE1_layer,'r'))
		layer2marks=AlignMarks()
		layer2marks.inverse = True
		layer2marks.arrowNumber = 10
		layer2marks.crossLength = 300
		waferElements.add(layer2marks.make(SPUT1_layer,'rl'))
		waferElements.add(layer2marks.make(SPUT2_layer,'rl'))
		layer3marks=AlignMarks()
		layer3marks.squareSize = 190
		layer3marks.arrowNumber = 9
		layer3marks.arrowStartInd = 3
		layer3marks.arrowAngle = 45
		layer3marks.crossLength = 300
		layer3marks.mirrorY = True
		waferElements.add(layer3marks.make(DRIE2_layer,'rl'))
		layer4marks=AlignMarks()
		layer4marks.inverse = True
		layer4marks.arrowNumber = 10
		layer4marks.crossLength = 300
		waferElements.add(layer4marks.make(OTS_MASK_layer,'rl'))
		
		layer5marks=AlignMarks()
		layer5marks.squareSize = 190
		layer5marks.arrowNumber = 9
		layer5marks.arrowStartInd = 3
		layer5marks.arrowAngle = 45
		layer5marks.crossLength = 300
		layer5marks.mirrorY = True
		waferElements.add(layer5marks.make(CHANETCH_layer,'rl'))
		
		for mark in self.alignMarks:
			waferElements.add(mark)
		self.waferArea = shapes.Disk((0,0), self.waferRadius, layer=LAYOUT_layer, number_of_points = 200)
		self.implArea = shapes.Disk((0,0),self.implanAreaRad, layer=LAYOUT_layer, number_of_points = 200)
		
		flatHeight = self.waferRadius-np.sqrt(self.waferRadius**2-(self.flatLength/2)**2)
		flatBox = shapes.Rectangle((-self.flatLength/2,-self.waferRadius), (self.flatLength/2,-self.waferRadius+flatHeight),layer=LAYOUT_layer)
		
		for layer in self.featuresLayers:
			waferElements.add(self.makeBoxes(layer))
	
		waferElements.add(self.waferArea)
		waferElements.add(self.implArea)
		waferElements.add(flatBox)
		
		return waferElements


wafer = Wafer()
wafer.featuresLayers = [DRIE1_layer,DRIE2_layer,SPUT1_layer,SPUT2_layer,OTS_MASK_layer,CHANETCH_layer]
waferElements = wafer.make()

####################Make default cells##################################
#Make cell 1
cellObj = Cell()
cellObj.name = "StdCell"
cellObj.cavSeparation = 600
cellObj.optCavObj.packing = False
cellObj.dispCavObj.packing = False
# cellObj.heaterObj = Heater(cellObj)
# cellObj.heaterObj.traceObj = RoundSerpentine(cellObj)
# cellObj.heaterObj.traceObj.thickness = 100e-9
# cellObj.heaterObj.traceObj.material = Materials.platinum
# cellObj.heaterObj.traceObj.width = 2690
# cellObj.heaterObj.traceObj.windNumber = 1
# cellObj.heaterObj.traceObj.windSeparation = 500
# cellObj.heaterObj.traceObj.centralSep = 80
# cellObj.heaterObj.traceObj.wireWidth = 40
# cellObj.heaterObj.traceObj.lineSpacing = 25
# cellObj.heaterObj.connectorObj = SpecialConnector(cellObj)
# cellObj.heaterObj.connectorObj.padSeparation = 10
# cellObj.heaterObj.padPairPos = np.array([1140,0])
# cellObj.sensorObj = Sensor(cellObj)
# cellObj.sensorObj.traceObj = RoundSerpentine(cellObj)
# cellObj.sensorObj.traceObj.thickness = 100e-9
# cellObj.sensorObj.traceObj.material = Materials.platinum
# cellObj.sensorObj.traceObj.width = 2530
# cellObj.sensorObj.traceObj.windNumber = 7
# cellObj.sensorObj.traceObj.windSeparation = 40
# cellObj.sensorObj.traceObj.centralSep = 20
# cellObj.sensorObj.traceObj.wireWidth = 10
# cellObj.sensorObj.traceObj.lineSpacing = 10
# cellObj.sensorObj.traceObj.mirror = True
# cellObj.sensorObj.connectorObj = SpecialConnector(cellObj)
# cellObj.sensorObj.connectorObj.align = "left"
# cellObj.sensorObj.connectorObj.offsetAB = 140
# cellObj.sensorObj.padPairPos = np.array([2020,0])
cellObj.make()

cellObj1 = copy.deepcopy(cellObj)
cellObj1.optCavObj.otsMaskOffset = 0
cellObj1.make()

cellObj5 = copy.deepcopy(cellObj)
cellObj5.optCavObj.radius = 1250
cellObj5.cavSeparation = 1600-1250
cellObj5.optCavObj.otsMaskOffset = 300
cellObj5.make()

cellObj2 = copy.deepcopy(cellObj)
cellObj2.optCavObj.radius = 750
cellObj2.cavSeparation = 1600-800
cellObj2.make()

cellObj3 = copy.deepcopy(cellObj)
cellObj3.optCavObj.radius = 500
cellObj3.cavSeparation = 1600-500
cellObj3.make()

cellObj4 = copy.deepcopy(cellObj)
cellObj4.optCavObj.radius = 1500
cellObj4.cavSeparation = 300
cellObj4.optCavObj.otsMaskOffset = 200
cellObj4.make()

####################Make alternative cells##############################

def makeStarCells():
	
	def makeStareCellX(frac1,frac2,frac3,DR):
		starCellX=copy.deepcopy(cellObj)
		starCellX.optCavObj = StarRing(starCellX)
		starCellX.name = "StarCell"+"F1_"+str(frac1)+"DR_"+str(DR)
		starCellX.optCavObj.DR=DR
		starCellX.optCavObj.frac1=frac1
		starCellX.optCavObj.frac2=frac2
		starCellX.optCavObj.frac3=frac3
		starCellX.cavSeparation = 500
		return starCellX
		
	starCells = []
	starCell1 =  makeStareCellX(60,1,1000,300)
	starCell1.optCavObj.packingWidth=200
	# starCell1.cavSeparation = 500
	starCell1.make()
	starCells.append(starCell1)

	starCell2=makeStareCellX(80,1,1000,100)
	starCell2.optCavObj.packingWidth=250
	starCell2.make()
	starCells.append(starCell2)

	starCell3=makeStareCellX(40,1,1000,100)
	starCell3.make()
	starCells.append(starCell3)

	starCell4=makeStareCellX(80,1,1000,0)
	starCell4.make()
	starCells.append(starCell4)

	starCell5=copy.deepcopy(starCell2)
	starCell5.optCavObj.radius=1100
	starCell5.optCavObj.DR=0
	starCell5.make()
	starCells.append(starCell5)
	
	return starCells

def makeThermCells():
	thermCells = []
	thermCell0=copy.deepcopy(cellObj)
	thermCell0.name = "ThermCell0"
	thermCell0.heaterObj.connectorObj.wireWidth = 80
	thermCell0.heaterObj.connectorObj.lineSpacing = 60
	thermCell0.heaterObj.traceObj.width = 3080
	thermCell0.heaterObj.traceObj.windNumber = 2
	thermCell0.heaterObj.traceObj.windSeparation = 440
	thermCell0.heaterObj.traceObj.wireWidth = 80
	thermCell0.heaterObj.traceObj.lineSpacing = 45
	thermCell0.heaterObj.traceObj.centralSep = 120
	thermCell0.sensorObj.connectorObj.offsetAB = 150
	thermCell0.sensorObj.traceObj.mirror = False
	thermCell0.sensorObj.traceObj.width = 2840
	thermCell0.sensorObj.traceObj.windNumber = 6
	thermCell0.sensorObj.traceObj.windSeparation = 40
	thermCell0.sensorObj.traceObj.centralSep = 120
	thermCell0.sensorObj.traceObj.separationXoffset = 120
	thermCell0.make()
	thermCells.append(thermCell0)

	thermCell1=copy.deepcopy(cellObj)
	thermCell1.name = "ThermCell1"
	thermCell1.optCavObj = ThermIsolCavs(thermCell1)
	thermCell1.optCavObj.wallThickness = 400
	thermCell1.optCavObj.dispCavCenter = np.array([0,-2200])
	thermCell1.chanNumber = 0
	thermCell1.dispCavObj = None
	thermCell1.textSize = 600
	thermCell1.heaterObj.traceObj = RoundSerpentine(thermCell1)
	thermCell1.heaterObj.traceObj.thickness = 100e-9
	thermCell1.heaterObj.traceObj.material = Materials.platinum
	thermCell1.heaterObj.traceObj.width = 2690
	thermCell1.heaterObj.traceObj.windNumber = 1
	thermCell1.heaterObj.traceObj.windSeparation = 500
	thermCell1.heaterObj.traceObj.centralSep = 80
	thermCell1.heaterObj.traceObj.wireWidth = 40
	thermCell1.heaterObj.traceObj.lineSpacing = 25
	thermCell1.heaterObj.connectorObj = SpecialConnector(thermCell1)
	thermCell1.heaterObj.connectorObj.padSeparation = 10
	thermCell1.heaterObj.padPairPos = np.array([1140,0])
	thermCell1.sensorObj = Sensor(thermCell1)
	thermCell1.sensorObj.traceObj = RoundSerpentine(thermCell1)
	thermCell1.sensorObj.traceObj.thickness = 100e-9
	thermCell1.sensorObj.traceObj.material = Materials.platinum
	thermCell1.sensorObj.traceObj.width = 2530
	thermCell1.sensorObj.traceObj.windNumber = 7
	thermCell1.sensorObj.traceObj.windSeparation = 40
	thermCell1.sensorObj.traceObj.centralSep = 20
	thermCell1.sensorObj.traceObj.wireWidth = 10
	thermCell1.sensorObj.traceObj.lineSpacing = 10
	thermCell1.sensorObj.traceObj.mirror = True
	thermCell1.sensorObj.connectorObj = SpecialConnector(thermCell1)
	thermCell1.sensorObj.connectorObj.align = "left"
	thermCell1.sensorObj.connectorObj.offsetAB = 140
	thermCell1.sensorObj.padPairPos = np.array([2020,0])
	thermCell1.make()
	thermCells.append(thermCell1)

	thermCell2=copy.deepcopy(thermCell1)
	thermCell2.name = "ThermCell2"
	thermCell2.optCavObj.bridgeWidth = 100
	thermCell2.make()
	thermCells.append(thermCell2)

	thermCell3=copy.deepcopy(thermCell2)
	thermCell3.name = "ThermCell3"
	thermCell3.optCavObj.topBridgeAngles= np.array([-4/6.,-2/6.,0,4/6.,2/6.])
	thermCell3.make()
	thermCells.append(thermCell3)

	thermCell4=copy.deepcopy(thermCell2)
	thermCell4.name = "ThermCell4"
	thermCell4.optCavObj.isolationCavs = False
	thermCell4.make()
	thermCells.append(thermCell4)

	thermCell5=copy.deepcopy(thermCell4)
	thermCell5.name = "ThermCell5"
	thermCell5.optCavObj.dispCavCenter = np.array([0,-2300])
	thermCell5.heaterObj.traceObj.lineSpacing = 60
	thermCell5.heaterObj.traceObj.width = 2760
	thermCell5.heaterObj.traceObj.centralSep = 115
	thermCell5.sensorObj.connectorObj.offsetAB = 210
	thermCell5.make()
	thermCells.append(thermCell5)

	thermCell6=copy.deepcopy(thermCell4)
	thermCell6.name = "ThermCell6"
	thermCell6.optCavObj.dispCavCenter = np.array([0,-2400])
	thermCell6.heaterObj.traceObj.lineSpacing = 100
	thermCell6.heaterObj.traceObj.width = 2840
	thermCell6.heaterObj.traceObj.centralSep = 150
	thermCell6.heaterObj.connectorObj.lineSpacing = 100
	thermCell6.sensorObj.connectorObj.offsetAB = 285
	thermCell6.make()
	thermCells.append(thermCell6)

	thermCell7=copy.deepcopy(thermCell4)
	thermCell7.name = "ThermCell7"
	thermCell7.heaterObj.connectorObj.lineSpacing = 70
	thermCell7.heaterObj.connectorObj.reductorLength = 245
	thermCell7.heaterObj.connectorObj.wireWidth = 100
	thermCell7.sensorObj.connectorObj.offsetAB = 285
	thermCell7.make()
	thermCells.append(thermCell7)

	thermCell8=copy.deepcopy(thermCell1)
	thermCell8.name = "ThermCell8"
	thermCell8.optCavObj.dispCavCenter = np.array([0,-2000])
	thermCell8.heaterObj.traceObj = ConformSerpentine(thermCell8)
	thermCell8.heaterObj.traceObj.wireWidth = 60
	thermCell8.heaterObj.traceObj.lineSpacing = 35
	thermCell8.heaterObj.traceObj.radialOffset = -95
	thermCell8.heaterObj.traceObj.outConnLength = 0
	thermCell8.sensorObj.traceObj = ConformSerpentine(thermCell8)
	thermCell8.sensorObj.traceObj.mirror = True
	thermCell8.sensorObj.traceObj.centralSep = 20
	thermCell8.sensorObj.traceObj.windNumber = 5
	thermCell8.sensorObj.traceObj.wireWidth = 10
	thermCell8.sensorObj.traceObj.lineSpacing = 10
	thermCell8.sensorObj.traceObj.radialOffset = -200
	thermCell8.sensorObj.traceObj.outConnLength = 0
	thermCell8.sensorObj.connectorObj.offsetAB = 0
	thermCell8.make()
	thermCells.append(thermCell8)

	thermCell9=copy.deepcopy(thermCell8)
	thermCell9.name = "ThermCell9"
	thermCell9.optCavObj.bridgeWidth = 100
	thermCell9.make()
	thermCells.append(thermCell9)

	thermCell10=copy.deepcopy(thermCell9)
	thermCell10.name = "ThermCell10"
	thermCell10.optCavObj.topBridgeAngles= np.array([-4/6.,-2/6.,0,4/6.,2/6.])
	thermCell10.make()
	thermCells.append(thermCell10)

	thermCell11=copy.deepcopy(thermCell9)
	thermCell11.name = "ThermCell11"
	thermCell11.optCavObj.isolationCavs = False
	thermCell11.make()
	thermCells.append(thermCell11)

	thermCell12=copy.deepcopy(thermCell11)
	thermCell12.name = "ThermCell12"
	thermCell12.optCavObj.dispCavCenter = np.array([0,-2000])
	thermCell12.heaterObj.traceObj.radialOffset = -300
	thermCell12.heaterObj.traceObj.outConnLength = 0
	thermCell12.heaterObj.traceObj.centralSep = 60
	thermCell12.sensorObj.traceObj.windNumber = 5
	thermCell12.sensorObj.traceObj.radialOffset = -45
	thermCell12.sensorObj.traceObj.outConnLength = 0
	thermCell12.sensorObj.traceObj.centralSep = 120
	thermCell12.sensorObj.traceObj.separationXoffset = 60
	thermCell12.sensorObj.connectorObj.offsetAB = 0
	thermCell12.make()
	thermCells.append(thermCell12)
	
	return thermCells

def makeVolVarCells():
	volVarCells=[]
	volVarCell1=copy.deepcopy(cellObj)
	volVarCell1.name = "volVarCell1"
	# volVarCell1.heaterObj = None #Temp
	# volVarCell1.sensorObj = None #Temp
	volVarCell1.textBool = False
	volVarCell1.optCavObj = VolVarCavs(volVarCell1)
	volVarCell1.chanNumber = 0
	volVarCell1.dispCavObj = None

	jiRange = [	(0,0),(0,4),(0,8),(0,11),(0,16),
				(1,0),(1,4),(1,9),(1,14),
				(2,0),(2,4),(2,8),(2,11),
				(2.7,0),(2.7,3),(2.8,5),(2.7,9),
				(3.4,0),(3.4,3),(3.4,7),
				(4.6,0),(4.5,4),
				(5.1,2)]

	def makeVolVarCellX(i,j):
		volVarCellX=copy.deepcopy(volVarCell1)
		volVarCellX.name = "volVarCellXi"+str(i)+"j"+str(j)
		volVarCellX.optCavObj.bottomCropLimit = 1000*j
		volVarCellX.optCavObj.finNumber = i
		return volVarCellX

	volVarCellX = makeVolVarCellX(16,0)
	volVarCellX.optCavObj.botChan = True
	volVarCellX.name = volVarCellX.name +"bis"
	volVarCellX.make()
	volVarCells.append(volVarCellX)

	volVarCellX = makeVolVarCellX(0,0)
	volVarCellX.optCavObj.botChan = True
	volVarCellX.name = volVarCellX.name +"bis"
	volVarCellX.make()
	volVarCells.append(volVarCellX)

	for (j,i) in jiRange:
		volVarCellX = makeVolVarCellX(i,j)
		volVarCellX.make()
		volVarCells.append(volVarCellX)

	for i in [0,2,4,6,8,10,12,16,24,36]:
		volVarCellY=copy.deepcopy(volVarCell1)
		volVarCellY.name = "volVarCellYi"+str(i)
		volVarCellY.optCavObj.radius = 500+i*100
		volVarCellY.optCavObj.bottomCropLimit = 6000
		volVarCellY.optCavObj.finNumber = 0
		volVarCellY.make()
		volVarCells.append(volVarCellY)

	volVarCellZ=copy.deepcopy(volVarCell1)
	volVarCellZ.name = "volVarCellZ"
	volVarCellZ.optCavObj.radius = 500
	volVarCellZ.optCavObj.bottomCropLimit = 0
	volVarCellZ.optCavObj.finNumber = 18
	volVarCellZ.make()
	volVarCells.append(volVarCellZ)
	return volVarCells

def makeBreakSealCells():
	breakSealCells=[]
	breakSealCell1=copy.deepcopy(cellObj)
	breakSealCell1.name = "BreakSealCell"
	breakSealCell1.optCavObj = BreakSealCavs(breakSealCell1)
	breakSealCell1.chanNumber = 0
	breakSealCell1.dispCavObj = None
	# breakSealCell1.make()
	
	# breakSealCell3=copy.deepcopy(breakSealCell1)
	# breakSealCell3.textBool = False
	# breakSealCell3.optCavObj.optCavCropLimit = 4000
	# breakSealCell3.optCavObj.southHeight = 600
	# breakSealCell3.make()
	def makeBreakSealCellX(thickness,nsWidth):
		breakSealCellX = copy.deepcopy(breakSealCell1)
		breakSealCellX.name = "BreakSealCellX"+str(i)
		breakSealCellX.optCavObj.breakSealThickness = thickness
		breakSealCellX.optCavObj.northWidth = nsWidth
		breakSealCellX.optCavObj.southWidth = nsWidth
		return breakSealCellX

	for i in [10,20,30,40,50,60,80,100,120]:
		breakSealCellX = makeBreakSealCellX(i,300)
		breakSealCellX.make()
		breakSealCells.append(breakSealCellX)

	for i in [30,40,50,60,80,100,120]:
		breakSealCellX = makeBreakSealCellX(i,600)
		breakSealCellX.make()
		breakSealCells.append(breakSealCellX)
		
	for i in [30,40,50,60,80,100,120]:
		breakSealCellX = makeBreakSealCellX(i,1200)
		breakSealCellX.make()
		breakSealCells.append(breakSealCellX)
		
	for i in [30,40,50,60,80,100,120]:
		breakSealCellX = makeBreakSealCellX(i,1200)
		breakSealCellX.optCavObj.northSeparation = 250
		breakSealCellX.optCavObj.dispCavSize = 2400
		breakSealCellX.optCavObj.dispCavYoffset = 600+i
		breakSealCellX.make()
		breakSealCells.append(breakSealCellX)

	for i in [30,40,50,60,80,100,120]:
		breakSealCellX = makeBreakSealCellX(i,300)
		breakSealCellX.optCavObj.sealThickness = 50
		breakSealCellX.optCavObj.northSeparation = 250
		breakSealCellX.optCavObj.dispCavSize = 2400
		breakSealCellX.optCavObj.dispCavYoffset = 600+i
		breakSealCellX.make()
		breakSealCells.append(breakSealCellX)

	for i in [30,40,50,60,80,100,120]:
		breakSealCellX = makeBreakSealCellX(i,300)
		breakSealCellX.optCavObj.sealThickness = 50
		breakSealCellX.optCavObj.northSeparation = 250
		breakSealCellX.optCavObj.dispCavSize = 2400
		breakSealCellX.optCavObj.dispCavYoffset = 450+i
		breakSealCellX.make()
		breakSealCells.append(breakSealCellX)
	return breakSealCells

#####################Define arrays######################################

#Make array 1
arrayObj = Array(wafer)
arrayObj.defaultCellObj = cellObj
arrayObj.columnNumber = 21
arrayObj.rowNumber = 13
arrayObj.leftDicMarkInd = [4]
arrayObj.rightDicMarkInd = [16]
arrayObj.topDicMarkInd = [2]
arrayObj.bottomDicMarkInd = [10]

cells = [cellObj,cellObj2,cellObj3,cellObj4,cellObj5]
labels = polyLabelRange([("BE","BK"),("AG","AK"),("LE","LK"),("MG","MK")])
arrayObj.specialCells.extend(fill(labels,cells,random=True,symetrise=True,colNum=arrayObj.columnNumber))

cells2 = [cellObj,cellObj1]
labels2 = polyLabelRange([("CC","CK")])
arrayObj.specialCells.extend(fill(labels2,cells2,random=True,symetrise=True,colNum=arrayObj.columnNumber))

# starCells = makeStarCells()
# labels = polyLabelRange([("AG","AO"),("MG","MO")])
# arrayObj.specialCells.extend(fill(labels,starCells,random=True))

# volVarCells = makeVolVarCells()
# labels = polyLabelRange([("DB","DT"),("JB","JT")])
# arrayObj.specialCells.extend(fill(labels,volVarCells,random=True))

# breakSealCells = makeBreakSealCells()
# labels = polyLabelRange([("FA","FU"),("GA","GB"),("GT","GU"),("HA","HU"),"GD","GR"])
# arrayObj.specialCells.extend(fill(labels,breakSealCells,random=True))

# thermCells = makeThermCells()
# labels = polyLabelRange([("BE","BK"),("LE","LK")])
# arrayObj.specialCells.extend(fill(labels,thermCells,random=True,symetrise=True,colNum=arrayObj.columnNumber))


#####################Make arrays############################
#Disable all special cells
# arrayObj.specialCells = []
array, arrayElement = arrayObj.make()

layoutNode = ET.Element("Layout")
layoutNode.append(cellListNode)
layoutNode.append(arrayElement)
tree = ET.ElementTree(layoutNode)
tree.write("output.xml", pretty_print=True)

#Make Top and Layout
top = core.Cell('TOP')
group = core.Cell('GROUP')
group.add(waferElements)
group.add(array)
top.add(group)
#Add wafer design to PDF
myShow(top,"TOP")
# amarks = core.GdsImport("getter-cavity-mask-centered2.gds",layers={3:GLASSETCH_layer})
# object = amarks.items()[0][1]
# top.add(object)

layout = core.Layout('LAYOUT')
layout.add(top)
filename='output.gds'
layout.save(filename)


os.startfile(filename)
pp.close()
print "Done. You can reload the layout."