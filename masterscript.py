import os
from os import mkdir
from os.path import exists
from csbdeep.utils import normalize
import matplotlib.pyplot as plt
import skimage
from skimage import io
from skimage.measure import regionprops_table
from stardist.models import StarDist2D
import glob, shutil
import numpy as np
import time
import math
import random

basefname = 'F:/automation_v3tg/'
minarea = 5000          # minimum area of ROI
maxint = 50000          # maximum mean pixel 
minint = 0
extraOC = True          # Whether or not to take snapshot in the extraOC channel
defaultROIsize = [150,150]
resizeROI = True        # Whether or not to resize the ROI to fit tightly around the cell
selectionMode = "largest"   # options: largest - chooses the largest cell
                            #           random - chooses a random cell

def movefile(basefname,snapnum):    
    print(f'Moving snap{snapnum}.')
    print('Do not close this window while the macro is running.')

    if snapnum == 1:
        currfiles = glob.glob(basefname + 'output/snaps/*')
    else:
        currfiles = glob.glob(basefname + 'output/snaps%d/*' % snapnum)

    maxint = 0
    for f in currfiles:
        f = f.split('\\')
        f = f[-1]
        digits = f[:-4]
        maxint = max(maxint,int(digits))

    if snapnum == 1:
        shutil.move(basefname + 'temp/snap.tif', 
                    basefname + 'output/snaps/%d.tif' % (maxint + 1))
        os.remove(f'{basefname}temp/done.txt')
    else:
        shutil.move(basefname + 'temp/snap%d.tif' % snapnum, 
                    basefname + 'output/snaps%d/%d.tif' % (snapnum, maxint + 1))
        os.remove(f'{basefname}temp/done{snapnum}.txt')
        
def movesmt(basefname):
    currfiles = glob.glob(f'{basefname}/output/smt/*')
    sourcefile = glob.glob(f'{basefname}/temp/smt_temp/*nd2')

    # get last index
    maxint = 0
    for f in currfiles:
        f = f.split('\\')
        f = f[-1]
        digits = f[:-4]
        maxint = max(maxint,int(digits))

    shutil.move(sourcefile[0], f'{basefname}/output/smt/{maxint + 1}.nd2')
    os.remove(f'{basefname}temp/donesmt.txt')

def delfile(basefname,snapnum):
    # delete snap file rather than moving
    print(f'No cell found. Deleting snap {snapnum}.')
    print('Do not close this window while the macro is running.')

    if snapnum == 1:
        os.remove(f'{basefname}temp/snap.tif')
        os.remove(f'{basefname}temp/done.txt')
    else:
        os.remove(f'{basefname}temp/snap{snapnum}.tif')
        os.remove(f'{basefname}temp/done{snapnum}.txt')

def locatecell(basefname,maxint,minarea,verbose=True,selectionMode=selectionMode,minint=minint):

    movestring = "StgMoveXY(%f, %f, 1);\n"
    
    # below are the pieces of the "movethere" macro
    movethere_part1 = f"""StgMoveXY(%f, %f, 1);
counter = counter + 1;
"counter";
counter;    
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{basefname}/temp/snap2.tif",18,0);
CloseAllDocuments(0);
WriteFile("{basefname}/temp/done2.txt","DONE",8);  

//WAIT FOR THE PYTHON SCRIPT TO PROCESS SNAP2
while(ExistFile("{basefname}/temp/snap2.tif")) {{
"Waiting for snap2 processing.";
}}

"""

    extraOCmac = f"""// TAKE AN IMAGE WITH AN EXTRA OC
SelectOptConf(extraOC);
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{basefname}/temp/snap4.tif",18,0);
WriteFile("{basefname}/temp/done4.txt","DONE",8);  
CloseAllDocuments(0);

"""

    movethere_part2 = f"""// SHRINK ROI TO CELL AND TAKE A THIRD IMAGE
SelectOptConf(cellFindingOC); // Set optical configuration here for finding cells.
RunMacro("{basefname}/temp/roiupdate.mac");
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{basefname}/temp/snap3.tif",18,0); 
WriteFile("{basefname}/temp/done3.txt","DONE",8);  
CloseAllDocuments(0);
    
// PRE-BLEACH
SelectOptConf(preBleachingOC); // Set optical configuration for pre-bleaching
CameraSet_Exposure(1, 7.000000);
LiveSync();
Wait(bleachTime);    
Freeze();
SelectOptConf(ISOC); // Set optical configuration for illumination sequence
CameraSet_Exposure(1, 7.000000);

// RUN ILLUMINATION SEQUENCE
IlluminationSequence_Run(ISname);    
CloseAllDocuments(0);
WriteFile("{basefname}/temp/donesmt.txt","DONE",8);  

// MOVE BACK TO GRID POINT CENTER
StgMoveXY(%f, %f, 1);"""

    # creates a pretrained model
    model = StarDist2D.from_pretrained('2D_versatile_fluo')

    #I = plt.imread('../output/snaps/3.tif')
    I = io.imread(f"{basefname}/temp/snap.tif", plugin='pil')
    labels, details = model.predict_instances(normalize(I), prob_thresh=0.5) 

    regiondata = skimage.measure.regionprops(labels)

    # get a bunch of properties of the ROIs
    rp = regionprops_table(labels, I,
                  properties=('centroid',
                             'orientation',
                             'axis_major_length',
                             'axis_minor_length',
                              'area',
                              'area_filled',
                              'eccentricity',
                              'intensity_mean',
                             )
                      )

    background = I[labels==0].mean()
    
    sel = ((rp['intensity_mean']-background<maxint) & (rp['area'] > minarea) & (rp['intensity_mean']-background>minint))

    if verbose:
        print('intensity_mean\t\tarea')
        for i in range(len(regiondata)):
            print(f"{rp['intensity_mean'][i]}\t\t{rp['area'][i]} ",end="")
            if sel[i]:
                print("***")
            else:
                print()

    large_regions = [];
    for i,region in enumerate(regiondata):    
      if sel[i]:
        large_regions.append(region)
        
    # if verbose:
        # print(large_regions)

    if not bool(large_regions):
      shiftx = 0
      shifty = 0
      with open(f'{basefname}/temp/movethere.mac', 'w') as fh:
        fh.write(movestring % (shiftx,-shifty))     
      return False # return value indicating that a suitable cell was not found

    else:

      # "largest_circular_region" is a bit of a misnomer, because I've added other possible selection modes
      if selectionMode == "largest":
        largest_circular_region = max(large_regions, key=lambda item: item.area)
      elif selectionMode == "random":
        largest_circular_region = max(large_regions, key=lambda item: random.random())
      else: # defaults to largest nucleus
        print("selectionMode not recognized! Defaulting to largest ROI.")
        largest_circular_region = max(large_regions, key=lambda item: item.area)

      centroid = largest_circular_region.centroid

      # remember that we used "valid" mode in the convolution, so indices are already shifted by half of the box width
      micronsperpx = 0.16 
      xcenter = centroid[1]
      ycenter = centroid[0]
      shiftx = micronsperpx * (xcenter - 256)
      shifty = micronsperpx * (ycenter - 256)
      
      # if extraOC:
        # mactemp = f'{basefname}/scripts/movethere_template_extra.mac'
      # else:
        # mactemp = f'{basefname}/scripts/movethere_template.mac'
        
      # OLD VERSION: Read in a macro and replaced the first and last lines
      # with open(mactemp, 'r') as fh:
        # macro_data = fh.readlines()
      # macro_data[0] = (movestring % (shiftx,shifty))
      # macro_data[-1] = (movestring % (-shiftx,-shifty))
      
      # New version: Use the strings above
      with open(f'{basefname}/temp/movethere.mac', 'w') as fh:
        if extraOC:
            fh.writelines(movethere_part1 % (shiftx,shifty) + extraOCmac + movethere_part2 % (-shiftx,-shifty))
        else:
            fh.writelines(movethere_part1 % (shiftx,shifty) + movethere_part2 % (-shiftx,-shifty))
      return True # return value indicating that a suitable cell was found

# def relocatecell(basefname,maxint,minarea,verbose=True,
        # resizeROI=resizeROI,selectionMode=selectionMode,defaultROIsize=defaultROIsize,
        # minint=minint):
def relocatecell(basefname,verbose=True,
        resizeROI=resizeROI,selectionMode=selectionMode,defaultROIsize=defaultROIsize):
    # creates a pretrained model
    model = StarDist2D.from_pretrained('2D_versatile_fluo')

    I = io.imread(f"{basefname}/temp/snap2.tif", plugin='pil')
    labels, details = model.predict_instances(normalize(I), prob_thresh=0.5) 

    # get a bunch of properties of the ROIs
    rp = regionprops_table(labels, I,
                  properties=('centroid',
                             'orientation',
                             'axis_major_length',
                             'axis_minor_length',
                              'area',
                              'area_filled',
                              'eccentricity',
                              'intensity_mean',
                             )
                      )

    background = I[labels==0].mean()

    # sel = ((rp['intensity_mean']-background<maxint) 
                # & (rp['area'] > minarea)
                # & (rp['intensity_mean']-background>minint)
                # & (rp['centroid-0'] > 199)
                # & (rp['centroid-0'] < 312)
                # & (rp['centroid-1'] > 199)
                # & (rp['centroid-1'] < 312))

    sel = ((rp['centroid-0'] > 199)
                & (rp['centroid-0'] < 312)
                & (rp['centroid-1'] > 199)
                & (rp['centroid-1'] < 312)
                & (rp['area'] > minarea))

    regiondata = skimage.measure.regionprops(labels)

    large_regions = []
    for i,region in enumerate(regiondata):    
      if sel[i]:
        large_regions.append(region)    
    
    if verbose:
        print(large_regions)

    if bool(large_regions):

        # "largest_circular_region" is a bit of a misnomer, because I've added other possible selection modes
        if selectionMode == "largest":
            largest_circular_region = max(large_regions, key=lambda item: item.area)
        elif selectionMode == "random":
            largest_circular_region = max(large_regions, key=lambda item: random.random())
        else: # defaults to largest nucleus
            print("selectionMode not recognized! Defaulting to largest ROI.")
            largest_circular_region = max(large_regions, key=lambda item: item.area)

        bbox = largest_circular_region.bbox
        centroid = largest_circular_region.centroid

        if resizeROI:
            minx = bbox[1]
            maxx = bbox[3]
            miny = bbox[0]
            maxy = bbox[2]
        else:
            meanx = bbox[1] + bbox[3]
            meany = bbox[0] + bbox[2]
            minx = meanx - math.floor(defaultROIsize[0]/2)
            maxx = meanx + math.ceil(defaultROIsize[0]/2)
            miny = meany - math.floor(defaultROIsize[1]/2)
            maxy = meany + math.ceil(defaultROIsize[1]/2)

    else:  # if a cell is not relocated, set the ROI by default to be a default sized ROI in the center of the FOV
        minx = 256 - math.floor(defaultROIsize[0]/2)
        maxx = 256 + math.ceil(defaultROIsize[0]/2)
        miny = 256 - math.floor(defaultROIsize[1]/2)
        maxy = 256 + math.ceil(defaultROIsize[1]/2)
        
        if verbose:
            print('Cell not relocated. Defaulting to centered 150 x 150 px ROI.')
        
    ROImacrostring = """ROISet(%d,%d,%d,%d);
    CameraFormatSet(1, "FMT 1x1 (X-7244) 16-bit");
    CameraFormatSet(2, "FMT 1x1 (X-7244) 16-bit");
    ROIEnable(1); """ % (minx,miny,maxx,maxy) 

    with open(f'{basefname}/output/rois.txt', 'a') as fh:
        fh.write('%d,%d,%d,%d\n' % (minx,miny,maxx,maxy))

    with open(f'{basefname}/temp/roiupdate.mac', 'w') as fh:
        fh.write(ROImacrostring)


if __name__ == "__main__":

    # make output directories

    outdirs = [f'{basefname}/{j}/' for j in ['output',
        'output/snaps',
        'output/snaps2',
        'output/snaps3',
        'output/smt',
        'temp',
        'temp/smt_temp']]

    for outdir in outdirs:
        if not exists(outdir):
            mkdir(outdir)
        
    # make snaps4 folder only if extraOC is selected
    if extraOC:
        outdir = f'{basefname}/output/snaps4/'
        if not exists(outdir):
            mkdir(outdir)
    
    # loop that will keep running to detect, process, and move files as they are created by NIS Elements 
    while True:
        if os.path.isfile(f'{basefname}/temp/done.txt'):
            # attempt to find a cell in the first snapshot. If it is found, move the snapshot and write out movethere.mac. 
            # Otherwise, delete the snapshot and write out a placeholder movethere.mac
            if locatecell(basefname,maxint,minarea):
                movefile(basefname,1)
            else:
                delfile(basefname,1)
        if os.path.isfile(f'{basefname}/temp/done2.txt'):
            # relocate cell in second snapshot and move second snapshot to snaps2 folder
            relocatecell(basefname)
            movefile(basefname,2)
        if os.path.isfile(f'{basefname}/temp/done3.txt'):
            try:
                movefile(basefname,3)
            except Exception as e:
                print(e)
        if os.path.isfile(f'{basefname}/temp/done4.txt'):
            try:
                movefile(basefname,4)
            except Exception as e:
                print(e)
        if os.path.isfile(f'{basefname}/temp/donesmt.txt'):
            movesmt(basefname)
    
# TODO: 
# Add a stop signal from the macro to break out of the while loop, clean up temp files, and exit python program














