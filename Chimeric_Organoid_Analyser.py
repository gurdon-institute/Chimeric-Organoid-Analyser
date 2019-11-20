# Segment organoids in DIC images, measure contained Tomato and GFP signal in each organoid.
# Minima are set to return no area for negative organoids with auto-thresholding.
#
# 	-by Richard Butler, Gurdon Institute Imaging Facility


import math as maths

from ij import IJ, ImagePlus
from ij.process import ImageProcessor, ByteProcessor, AutoThresholder, Blitter, FloodFiller
from ij.plugin.filter import ThresholdToSelection, RankFilters
from ij.measure import ResultsTable
from ij.gui import PolygonRoi, ShapeRoi, TextRoi, Overlay

from java.awt import Color, Font


FONT = Font(Font.SANS_SERIF, Font.BOLD, 42)

BFminA = 2500		#µm²
tomatoMinA = 100
GFPminA = 500


def fillHoles(ip):
	width = ip.getWidth()
	height = ip.getHeight()
	ff = FloodFiller(ip)
	ip.setColor(127)
	foreground = 127
	background = 0
	for y in range(height):
	    if ip.getPixel(0,y)==background:
	    	ff.fill(0, y)
	    if ip.getPixel(width-1,y)==background:
	    	ff.fill(width-1, y)
	for x in range(width):
	    if ip.getPixel(x,0)==background:
	    	ff.fill(x, 0)
	    if ip.getPixel(x,height-1)==background:
	    	ff.fill(x, height-1)
	n = width*height
	for i in range(n):
		if ip.get(i)==127:
		    ip.set(i, 0)
		else:
		    ip.set(i, 255)

def getMask(ip, method):
	W = ip.getWidth()
	H = ip.getHeight()
	
	hist = ip.getHistogram(256)
	stats = ip.getStatistics()
	
	thresh = AutoThresholder().getThreshold( method, hist )
	thresh = (thresh/float(255)) * stats.max

	mask = ByteProcessor(W,H)
	for i in range(W*H):
		value = ip.getf(i)
		bite = 255 if value>=thresh else 0
		mask.set(i, bite)

	fillHoles(mask)

	mask.erode()
	mask.dilate()
	
	return mask

def getRois(mask):
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	composite = ThresholdToSelection().convert(mask)
	rois = ShapeRoi(composite).getRois()
	return rois


imp = IJ.getImage()
cal = imp.getCalibration()
stack = imp.getStack()
z = imp.getSlice()

ipGFP = stack.getProcessor( imp.getStackIndex(1, z, 1) ).convertToFloatProcessor()
ipGFP.blurGaussian(4)
maskGFP = getMask(ipGFP, AutoThresholder.Method.Otsu)
roisGFP = getRois(maskGFP)

ipBF = stack.getProcessor( imp.getStackIndex(2, z, 1) ).convertToFloatProcessor()
RankFilters().rank(ipBF, 1.0, RankFilters.VARIANCE)
ipBF.blurGaussian(5)
maskBF = getMask(ipBF, AutoThresholder.Method.Triangle)
roisBF = getRois(maskBF)

ipTomato = stack.getProcessor( imp.getStackIndex(3, z, 1) ).duplicate()
sub = ipTomato.duplicate()
ipTomato.blurGaussian(5)
sub.blurGaussian(20)
ipTomato.copyBits(sub, 0,0, Blitter.SUBTRACT)
maskTomato = getMask(ipTomato, AutoThresholder.Method.MaxEntropy)
roisTomato = getRois(maskTomato)

bfMeasure = stack.getProcessor( imp.getStackIndex(2, z, 1) )
tomatoMeasure = stack.getProcessor( imp.getStackIndex(3, z, 1) )
gfpMeasure = stack.getProcessor( imp.getStackIndex(1, z, 1) )
ol = Overlay()
rt = ResultsTable()
rt.showRowNumbers(False)
for bf in roisBF:
	bfMeasure.setRoi(bf)
	roiStatsBF = bfMeasure.getStatistics()
	if roiStatsBF.area * cal.pixelWidth * cal.pixelHeight < BFminA:
		continue
	bounds = bf.getBounds()
	cXbf = int( bounds.x+(bounds.width/2.0) )
	cYbf = int( bounds.y+(bounds.height/2.0) )

	bf.setStrokeColor(Color.YELLOW)
	bf.setPosition(0,z,1)
	ol.add(bf)
	
	tomatoA = 0
	for tomato in roisTomato:
		tomatoMeasure.setRoi(tomato)
		roiStatsTomato = tomatoMeasure.getStatistics()
		if roiStatsTomato.area * cal.pixelWidth * cal.pixelHeight < tomatoMinA:
			continue
		bounds = tomato.getBounds()
		cXtomato = int( bounds.x+(bounds.width/2.0) )
		cYtomato = int( bounds.y+(bounds.height/2.0) )
		if bf.contains(cXtomato, cYtomato):
			tomatoA += roiStatsTomato.area
			tomato.setStrokeColor(Color.CYAN)
			tomato.setPosition(0,z,1)
			ol.add(tomato)
	
	gfpA = 0
	gfpSum = 0
	for gfp in roisGFP:
		gfpMeasure.setRoi(gfp)
		roiStatsGfp = gfpMeasure.getStatistics()
		if roiStatsGfp.area * cal.pixelWidth * cal.pixelHeight < GFPminA:
			continue
		bounds = gfp.getBounds()
		cXgfp = int( bounds.x+(bounds.width/2.0) )
		cYgfp = int( bounds.y+(bounds.height/2.0) )
		if bf.contains(cXgfp, cYgfp):
			gfpA += roiStatsGfp.area
			gfpSum += roiStatsGfp.mean * roiStatsGfp.area
			gfp.setStrokeColor(Color.MAGENTA)
			gfp.setPosition(0,z,1)
			ol.add(gfp)
		
	gfpMean = gfpSum / gfpA if gfpA>0 else 0

	row = rt.getCounter()

	uni = row%26
	cha = chr(ord('A')+uni)
	name = ""
	while len(name) < (row+1)/26.0:
		name += cha

	label = TextRoi(name, cXbf-20, cYbf+20, FONT)
	label.setStrokeColor(Color.GREEN)
	label.setPosition(0,z,1)
	ol.add(label)

	rt.setValue("Organoid", row, name)
	rt.setValue("X", row, cXbf*cal.pixelWidth)
	rt.setValue("Y", row, cYbf*cal.pixelHeight)
	rt.setValue("Area", row, roiStatsBF.area*cal.pixelWidth*cal.pixelHeight)
	rt.setValue("Contained Tomato Area", row, tomatoA*cal.pixelWidth*cal.pixelHeight)
	rt.setValue("Contained GFP Area", row, gfpA*cal.pixelWidth*cal.pixelHeight)
	rt.setValue("Contained GFP Mean", row, gfpMean)
	
imp.setOverlay(ol)
rt.show(imp.getTitle()+"_Z"+str(z)+"_Results")


