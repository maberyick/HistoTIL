import subprocess
from timeit import default_timer as timer
import argparse
import os, glob, sys
from openslide import open_slide, deepzoom
import cv2
from PIL import Image, ImageStat
import concurrent.futures
import time
from alive_progress import alive_bar
from PIL import Image, ImageFont, ImageDraw
import pandas as pd
import numpy as np

def main(args):
    columns = ['File']
    wsi_file_list = pd.read_csv(args.wsi_file_list, delimiter='\t', header=None, names=columns)
    files = wsi_file_list['File'].tolist()
    print(f'{len(files)} files have been read.')
    if args.end_ind == -1:
        args.end_ind = len(files)-1
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_worker) as executor:
        {executor.submit(startPatchExtraction, ind, files, args): ind for ind in range(args.start_ind, args.end_ind+1)}

def startPatchExtraction(ind, files, args):
    font = ImageFont.truetype(font="./pylib/Arca-Heavy.ttf", size=14)
    fill_color = (255, 255, 255)
    stroke_color = (0, 0, 0)
    maskFname = os.path.basename(files[ind])[:-4] + '.png'
    print(files[ind])
    maskPath = os.path.join(args.mask_path, maskFname)
    baseFname = os.path.splitext(os.path.basename(maskFname))[0]
    print('Image number: %s, Out of: %s, Image ID: %s' % (ind+1, len(files), baseFname))
    #print('Image number: %s, Out of: %s, Saving path: %s, Image ID: %s, From the path: %s' % (ind+1, len(files), args.output_path, baseFname, files[ind]))

    if args.save_summary_only == 0:
        if os.path.exists(os.path.join(args.output_path, baseFname)):
            print("Path exists")
            return
        else:
            os.makedirs(os.path.join(args.output_path, baseFname))
    #slideFname = os.path.join(args.input_path, baseFname + '.' + args.ext)
    #slideFname = os.path.join(files[ind])
    slideFname = files[ind]
    slide = open_slide(slideFname)
    dz = deepzoom.DeepZoomGenerator(slide, tile_size=args.tile_size, overlap=0, limit_bounds=True)
    # Size of the smallest size from the WSI to be used for tracking the extraction of patches
    maskLevel = slide.level_downsamples[-args.visual_size]

    mask_tile_size = args.tile_size*slide.level_downsamples[args.mag_level]//maskLevel
    dzLevel = dz.level_count-args.mag_level-1
    
    if not os.path.exists(maskPath):
        print(f"Mask file '{maskPath}' does not exist. Skipping to the next file.")
        return
    mask = cv2.imread(maskPath)
    mask = Image.fromarray(mask).convert("L")
    #mask_org_size = mask.size
    # size of the chosen mask
    slide_siz = slide.level_dimensions[-args.visual_size]
    mask = mask.resize(slide_siz)
    fn = lambda x : 1 if x > 0 else 0
    mask = mask.point(fn, mode='1')
    
    # get slide thumbnail
    img_patch = slide.get_thumbnail(slide_siz)
    img_patches = ImageDraw.Draw(img_patch)
    
    # counter ids
    counter = 0
    total_counter = 0
    
    # store the image IDs and basic info
    counter_id = []
    wsi_id = []
    wsi_tiss_per = []
    wsi_coord_x = []; wsi_coord_y = []
    wsi_coord_r = []; wsi_coord_c = []

    # use the tracking bar
    with alive_bar(dz.level_tiles[dzLevel][0]*dz.level_tiles[dzLevel][1], title=f'Extracting patches from {baseFname}') as bar:
        # loop through the i and j (positions) of the WSI
        for i in range(dz.level_tiles[dzLevel][0]):
            for j in range(dz.level_tiles[dzLevel][1]):
                bar()
                coord = dz.get_tile_coordinates(dzLevel, (i, j))     
                # if the patch is not square (given dimension), continue
                if coord[2] != (args.tile_size, args.tile_size):
                    total_counter += 1
                    continue
                else:
                    coord = coord[0]
                cenX = (coord[0]+args.tile_size*slide.level_downsamples[args.mag_level]//2)//maskLevel
                cenY = (coord[1]+args.tile_size*slide.level_downsamples[args.mag_level]//2)//maskLevel
                maskRegion = mask.crop((cenX-(mask_tile_size//2), cenY-(mask_tile_size//2), cenX+(mask_tile_size//2), cenY+(mask_tile_size//2)))
                # include the text
                label = '#'+str(counter+1)+'  '+'R'+str(j+1)+' - '+'C'+str(i+1)+'\n'+'X'+str(coord[1]+1)+'\n'+'Y'+str(coord[0]+1)+'\n'+str(round(100*ImageStat.Stat(maskRegion).mean[0]))+'%'
                # the coordinates
                coord_rect = [(cenX-(mask_tile_size//2), cenY-(mask_tile_size//2)), (cenX+(mask_tile_size//2), cenY+(mask_tile_size//2))]
                coord_text1 = [cenX-(mask_tile_size//2)+3, cenY-(mask_tile_size//2)+3]
                coord_text2 = [cenX-(mask_tile_size//2)+2, cenY-(mask_tile_size//2)+2]
                if ImageStat.Stat(maskRegion).mean[0] > args.min_cellular_density:
                    if args.save_summary_only == 0:
                        tile = dz.get_tile(dzLevel, (i, j)).convert("RGB")
                        if args.img_resize == 1:
                            tile.resize((args.img_resize_value,args.img_resize_value)).save(os.path.join(args.output_path,baseFname, baseFname+'_r'+str(coord[1])+'_c'+str(coord[0])+'.'+args.patch_ext))
                        else:
                            tile.resize((args.tile_size,args.tile_size)).save(os.path.join(args.output_path,baseFname, baseFname+'_r'+str(coord[1])+'_c'+str(coord[0])+'.'+args.patch_ext))
                    counter += 1
                    wsi_tiss_per.append(ImageStat.Stat(maskRegion).mean[0])
                    wsi_id.append(baseFname+'_r'+str(coord[1])+'_c'+str(coord[0]))
                    counter_id.append(counter)
                    wsi_coord_r.append(j+1)
                    wsi_coord_c.append(i+1)
                    wsi_coord_x.append(coord[1]+1)
                    wsi_coord_y.append(coord[0]+1)
                    if ImageStat.Stat(maskRegion).mean[0] >= 0.8:
                        img_patches.rectangle(coord_rect, fill=None, outline = (127,255,0), width=2)
                        img_patches.text(coord_text1, label, fill=stroke_color, stroke_fill=stroke_color, font=font, spacing=2)
                        img_patches.text(coord_text2, label, fill=fill_color, stroke_with=0, font=font, spacing=2)
                    else:
                        img_patches.rectangle(coord_rect, fill=None, outline = (255,215,0), width=2)
                        img_patches.text(coord_text1, label, fill=stroke_color, stroke_fill=stroke_color, font=font, spacing=2)
                        img_patches.text(coord_text2, label, fill=fill_color, stroke_with=0, font=font, spacing=2)
                else:
                    if ImageStat.Stat(maskRegion).mean[0] > args.min_cellular_density/2:
                        label = '# NA'+'  '+'R'+str(j+1)+' - '+'C'+str(i+1)+'\n'+'X'+str(coord[1]+1)+'\n'+'Y'+str(coord[0]+1)+'\n'+str(round(100*ImageStat.Stat(maskRegion).mean[0]))+'%'
                        img_patches.rectangle(coord_rect, fill=None, outline = (255,48,48), width=2)
                        img_patches.text(coord_text1, label, fill=stroke_color, stroke_fill=stroke_color, font=font, spacing=2)
                        img_patches.text(coord_text2, label, fill=fill_color, stroke_with=0, font=font, spacing=2)
        img_patch.save(args.output_path+'/patch_summary/'+baseFname+'.png', format='PNG')
        slide.close()
        df = pd.DataFrame(np.column_stack([counter_id,wsi_id,wsi_tiss_per,wsi_coord_x,wsi_coord_y,wsi_coord_r,wsi_coord_c]), columns=["Number","Patch name","Tissue percentage","WSI coord (x)", "WSI coord (y)","WSI position (row)", "WSI position (column)"])
        df.to_csv(args.output_path+'/patch_summary/'+baseFname+'_patch_list.csv', index=False)
        print('Image path: %s, Image ID: %s, Image mask path: %s, Number of patches saved: %s' % (slideFname, baseFname, maskFname, counter+1))

### Working paths
parser = argparse.ArgumentParser(description='configuration for running the histoQC')
parser.add_argument('--output_path', type=str, default='./testimage/patch_output/', help='output directory to save the results of the patches')
parser.add_argument('--ext', type=str, default='svs', help='extension of the wsi image (e.g. svs,tif,tiff)')
parser.add_argument('--patch_ext', type=str, default='png', help='extension of the patch image (png,jpeg)')
### Stats parameters
parser.add_argument('--min_cellular_density', type=float, default=0.4, help='Minimum value of celular density (quantity of cells (%) on a patch) required to save an image')
parser.add_argument('--max_worker', type=int, default=16, help='Number of cores to be used')
parser.add_argument('--tile_size', type=int, default=2048, help='Size of the patch in pixels (e.g. 20048x2048)')
parser.add_argument('--mag_level', type=int, default=0, help='Mag. level of the image (For 40x:0, For 20x:1, For 10x:2, For 5x:3)')
parser.add_argument('--visual_size', type=int, default=1, help='Relative size of the visual patch summary (1 = small, 2 = medium, 3 = big, 4 = large')
parser.add_argument('--start_ind', type=int, default=0, help='Index to start the process (0 start from the first wsi image)')
parser.add_argument('--end_ind', type=int, default=-1, help='Index to end the process (-1 ends at the last wsi image)')
parser.add_argument('--img_resize', type=int, default=-1, help='Resize the final patch image. Useful if image is 40x and wants to be 20x)')
parser.add_argument('--img_resize_value', type=int, default=2048, help='Value of the final patch (resized)')
### Added for including a tsv files
parser.add_argument('--wsi_file_list', type=str, default='./testimage/wsi_file_list.tsv', help='path to the TSV file containing the list of wsi file paths and names')
parser.add_argument('--mask_path', type=str, default='./testimage/histoqc_mask/', help='directory of the mask generated by histoqc')
parser.add_argument('--save_summary_only', type=int, default=0, help='If 1, just save the summary information of the patches, 0 if not and do the patches')


args = parser.parse_args()

if not os.path.isdir(args.output_path):
    os.mkdir(args.output_path)
if not os.path.isdir(args.output_path+'/patch_summary/'):
    os.mkdir(args.output_path+'/patch_summary/')

if __name__ == "__main__":
    start = timer()
    results = main(args)
    end = timer()
    print('Script Time: %f secons' % (end - start))