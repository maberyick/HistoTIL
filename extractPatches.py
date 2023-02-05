import os, glob, sys
from openslide import open_slide, deepzoom
import cv2
from PIL import Image, ImageStat
import concurrent.futures
import time

def startPatchExtraction(ind, files):
    maskFname = files[ind]
    baseFname = os.path.splitext(os.path.basename(maskFname))[0]
    print(ind, len(files), outPath, baseFname, dataPath)
    if os.path.exists(os.path.join(outPath, baseFname)):
        print("Path exists")
        return
    else:
        os.makedirs(os.path.join(outPath, baseFname))

    slideFname = os.path.join(dataPath, baseFname + ext)
    slide = open_slide(slideFname)
    dz = deepzoom.DeepZoomGenerator(slide, tile_size=tileSize, overlap=0, limit_bounds=True)
    maskLevel = 3
    maskTileSize = tileSize*slide.level_downsamples[imageLevel]//slide.level_downsamples[maskLevel]
    dzLevel = dz.level_count-imageLevel-1

    mask = cv2.imread(maskFname)
    mask = Image.fromarray(mask).convert("L")
    mask = mask.resize(slide.level_dimensions[maskLevel])
    fn = lambda x : 1 if x > 0 else 0
    mask = mask.point(fn, mode='1')

    counter = 0
    for i in range(dz.level_tiles[dzLevel][0]):
        for j in range(dz.level_tiles[dzLevel][1]):
            coord = dz.get_tile_coordinates(dzLevel, (i, j))            
            if coord[2] != (tileSize, tileSize):
                continue
            else:
                coord = coord[0]
            cenX = (coord[0]+tileSize*slide.level_downsamples[imageLevel]//2)//slide.level_downsamples[maskLevel]
            cenY = (coord[1]+tileSize*slide.level_downsamples[imageLevel]//2)//slide.level_downsamples[maskLevel]
            maskRegion = mask.crop((cenX-(maskTileSize//2), cenY-(maskTileSize//2), cenX+(maskTileSize//2), cenY+(maskTileSize//2)))
            if ImageStat.Stat(maskRegion).mean[0] > minCellularRegionDensityPerTile:
                counter += 1
                tile = dz.get_tile(dzLevel, (i, j)).convert("RGB")
                tile.save(os.path.join(outPath,baseFname, baseFname+'_r'+str(coord[1])+'_c'+str(coord[0])+patchExt))
    
    slide.close()
    print(slideFname, baseFname, maskFname, counter)


def startPatchExtractionWithoutMask(slides, ind):
    slideFname = slides[ind]
    baseFname = os.path.basename(slideFname)
    baseFname = baseFname[:-len(ext)]
    print(ind, len(files), baseFname)
    if os.path.exists(os.path.join(outPath, baseFname)) or os.path.exists(os.path.join(outPath, baseFname.replace('-', '_'))):
        print("Path already exists!. (" + outPath + "/" + baseFname + ")") 
        return
    else:
        #print("Path not  exists! Just created! (" + outPath + "/" + baseFname + ")")
        os.makedirs(os.path.join(outPath, baseFname))

    print(slideFname)
    slide = open_slide(slideFname)
    dz = deepzoom.DeepZoomGenerator(slide, tile_size=tileSize, overlap=0, limit_bounds=True)
    dzLevel = dz.level_count-imageLevel-1

    maskLevel = 3
    maskTileSize = tileSize*slide.level_downsamples[imageLevel]//slide.level_downsamples[maskLevel]

    mask = slide.read_region((0,0), maskLevel, slide.level_dimensions[maskLevel]).convert("L")
    fn1 = lambda x : 1 if x > 200 else 0
    fn2 = lambda x : 1 if x < 50 else 0
    mask1 = mask.point(fn1, mode='1')
    mask2 = mask.point(fn2, mode='1')

    counter = 0
    for i in range(dz.level_tiles[dzLevel][0]):
        for j in range(dz.level_tiles[dzLevel][1]):
            coord = dz.get_tile_coordinates(dzLevel, (i, j))            
            if coord[2] != (tileSize, tileSize):
                continue
            else:
                coord = coord[0]

            cenX = (coord[0]+tileSize*slide.level_downsamples[imageLevel]//2)//slide.level_downsamples[maskLevel]
            cenY = (coord[1]+tileSize*slide.level_downsamples[imageLevel]//2)//slide.level_downsamples[maskLevel]

            maskRegion1 = mask1.crop((cenX-(maskTileSize//2), cenY-(maskTileSize//2), cenX+(maskTileSize//2), cenY+(maskTileSize//2)))
            maskRegion2 = mask2.crop((cenX-(maskTileSize//2), cenY-(maskTileSize//2), cenX+(maskTileSize//2), cenY+(maskTileSize//2)))
            if ImageStat.Stat(maskRegion1).mean[0] < minCellularRegionDensityPerTile and ImageStat.Stat(maskRegion2).mean[0] < minCellularRegionDensityPerTile:
                counter += 1
                tile = dz.get_tile(dzLevel, (i, j)).convert("RGB")
                tile.save(os.path.join(outPath,baseFname,baseFname+'_r'+str(coord[1])+'_c'+str(coord[0])+patchExt))

    slide.close()
    print(slideFname, baseFname, counter)


if __name__ == "__main__":
    minCellularRegionDensityPerTile = 0.4
    max_workers = 16

    tileSize = int(sys.argv[1])
    imageLevel = int(sys.argv[2]) #0 -> 40x, 1 -> 20x, 2 -> 10x, 3 -> 5x
    dataPath = sys.argv[3]
    ext = sys.argv[4]
    outPath = sys.argv[5]
    patchExt = sys.argv[6]
    startInd = int(sys.argv[7])
    endInd = int(sys.argv[8])

    start = time.time()
    if len(sys.argv) == 10:
        maskPath = sys.argv[9]
        print(os.path.join(maskPath, '*.png'))
        files = sorted(glob.glob(os.path.join(maskPath, '*.png')))
        print(f'{len(files)} mask files have been read.', files)
        if endInd == -1:
            endInd = len(files)-1
        startPatchExtraction(0,files)
        #with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        #    {executor.submit(startPatchExtraction, ind, files): ind for ind in range(startInd, endInd+1)}
    else:
        files = sorted(glob.glob(os.path.join(dataPath, '*'+ext)))
        print(f'{len(files)} files have been read.')
        if endInd == -1:
            endInd = len(files)-1
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            {executor.submit(startPatchExtractionWithoutMask, files, ind): ind for ind in range(startInd, endInd+1)}

    finish = time.time()
    print("Patch extraction done!")
    print(f'Elapsed time: {finish-start}')


