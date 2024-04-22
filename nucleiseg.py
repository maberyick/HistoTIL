import os
import tensorflow as tf2
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import pylib.utils as utils
from PIL import Image
import numpy as np
from timeit import default_timer as timer
import argparse
import pandas as pd
import random

from histomicstk.preprocessing.color_normalization.\
    deconvolution_based_normalization import deconvolution_based_normalization

def main(args):
    tf.flags.DEFINE_string('gpu_index', '1', 'gpu index, default: 0')
    FLAGS = tf.flags.FLAGS
    flags = FLAGS
    run_config = tf.ConfigProto()
    sess = tf.Session(config=run_config)
    model = utils.net_arch(sess, flags, (2000,2000))
    saver = tf.train.Saver()
    sess.run(tf.global_variables_initializer())
    model_out_dir = "./models/model_nuclei"  
    # define color tcga image
    # TCGA-A2-A3XS-DX1_xmin21421_ymin37486_.png, Amgad et al, 2019)
    # for macenco (obtained using rgb_separate_stains_macenko_pca()
    # and reordered such that columns are the order:
    # Hamtoxylin, Eosin, Null
    W_target = np.array([
        [0.5807549,  0.08314027,  0.08213795],
        [0.71681094,  0.90081588,  0.41999816],
        [0.38588316,  0.42616716, -0.90380025]
    ])
    if utils.load_model(sess, saver, model_out_dir):
        best_auc_sum = sess.run(model.best_auc_sum)
        print('====================================================')
        print(' Best auc_sum: {:.3}'.format(best_auc_sum))
        print('=============================================')    
        print(' [*] Load Success!\n')
    if args.cohort_type == 'patch':
        names = utils.all_files_under(args.input_path,'.png')
        for name in names:
            print(name)
            print(os.path.basename(name))
            savingName = args.output_path + os.path.basename(name)[:-4]
            if os.path.isfile(savingName):
                print('file existed, continue to next')
                continue
            image = Image.open(name)
            print(np.shape(image))
            # ToDo: Improve the detection of cells, image scaling and color normalization
            stain_unmixing_routine_params = {
                'stains': ['hematoxylin', 'eosin'],
                'stain_unmixing_method': 'macenko_pca',}
            tissue_rgb_normalized = deconvolution_based_normalization(
                np.asarray(image), W_target=W_target,
                stain_unmixing_routine_params=stain_unmixing_routine_params)
            image = Image.fromarray(np.uint8(tissue_rgb_normalized))
            # Split image into 4 tiles
            #image_numpydata = asarray(image)
            #M = image_numpydata.shape[0]//2
            #N = image_numpydata.shape[1]//2
            #tiles = [image_numpydata[x:x+M,y:y+N] for x in range(0,image_numpydata.shape[0],M) for y in range(0,image_numpydata.shape[1],N)]
            # get image back
            #pilImage = Image.fromarray(tiles[0])
            #
            image = image.resize((2000,2000), Image.ANTIALIAS)
            print('predicting on image ', name)
            image = np.expand_dims(image, axis=0)
            print(np.shape(image))
            samples = sess.run(tf.get_default_graph().get_tensor_by_name('g_/Sigmoid:0'), \
                    feed_dict={tf.get_default_graph().get_tensor_by_name('image:0'): image})
            samples = np.squeeze(samples*255).astype(np.uint8)   
            # save as png file
            image = Image.fromarray(samples)
            image.save(savingName+'.png', format='PNG')
    elif args.cohort_type == 'folder':
        # list of the folders (cohorts)
        folder_list = utils.all_files_under(args.input_wsi)
        for folder in folder_list:
            folder_dir = os.path.splitext(folder)[0]
            folder_name = os.path.basename(folder_dir)
            if not os.path.isdir(args.input_path + folder_name +'/'):
                print("folder empty")
                continue
            print(folder_name)
            if not os.path.isdir(args.output_path + folder_name):
                os.mkdir(args.output_path + folder_name)
            # loop through the files of each folder
            names = utils.all_files_under(args.input_path+folder_name,'.png')
            for name in names:
                print(name)
                savingName = args.output_path + folder_name +'/'+ os.path.basename(name)[:-4]+'.png'
                if os.path.isfile(savingName):
                    print('file existed, continue to next')
                    continue
                image = Image.open(name)
                print(np.shape(image))
                image = image.resize((2000,2000), Image.ANTIALIAS) 
                print('predicting on image ', name)
                image = np.expand_dims(image, axis=0)
                print(np.shape(image))
                samples = sess.run(tf.get_default_graph().get_tensor_by_name('g_/Sigmoid:0'), \
                        feed_dict={tf.get_default_graph().get_tensor_by_name('image:0'): image})
                samples = np.squeeze(samples*255.0).astype(np.uint8)
                # save as png file
                image = Image.fromarray(samples)
                image.save(savingName, format='PNG')
    elif args.cohort_type == 'tsv_list':
        # list of the folders (cohorts)
        columns = ['File']
        wsi_file_list = pd.read_csv(args.wsi_file_list, delimiter='\t', header=None, names=columns)
        files = wsi_file_list['File'].tolist()
        print(f'{len(files)} files have been read.')
        if args.processing_order == 'random':
            random.shuffle(files)
        if args.processing_order == 'from end':
            files = reversed(files)    
        #folder_list = utils.all_files_under(args.input_wsi)
        for file in files:
            folder_dir = os.path.splitext(file)[0]
            folder_name = os.path.basename(folder_dir)
            #print(os.path.join(args.input_path, folder_name))
            if not os.path.isdir(os.path.join(args.input_path, folder_name)):
                #print("folder empty")
                continue
            #print(folder_name)
            if not os.path.isdir(os.path.join(args.output_path, folder_name)):
                os.mkdir(os.path.join(args.output_path, folder_name))
            # loop through the files of each folder
            names = utils.all_files_under(os.path.join(args.input_path,folder_name),'.png')
            for name in names:
                #print(name)
                savingName = os.path.join(args.output_path, folder_name, os.path.basename(name)[:-4]+'.png')
                if os.path.isfile(savingName):
                    #print('file existed, continue to next')
                    continue
                image = Image.open(name)
                print(np.shape(image))
                image = image.resize((2000,2000), Image.ANTIALIAS) 
                print('predicting on image ', name)
                image = np.expand_dims(image, axis=0)
                print(np.shape(image))
                samples = sess.run(tf.get_default_graph().get_tensor_by_name('g_/Sigmoid:0'), \
                        feed_dict={tf.get_default_graph().get_tensor_by_name('image:0'): image})
                samples = np.squeeze(samples*255.0).astype(np.uint8)
                # save as png file
                image = Image.fromarray(samples)
                image.save(savingName, format='PNG')
                
### Working paths
parser = argparse.ArgumentParser(description='configuration for running the nuclei segmentation')
parser.add_argument('--input_wsi', type=str, default='./testimage/wsi/', help='input directory of the wsi (folder option)')
parser.add_argument('--input_path', type=str, default='./testimage/patch/', help='input directory of the images (patch option)')
parser.add_argument('--output_path', type=str, default='./testimage/nucleimask/', help='output directory for nuclei mask')
parser.add_argument('--cohort_type', type=str, default='patch', help='Choose between patch or folder. running on a single folder (/cohort/*.png) of png images or a directory of folders (e.g. /cohorts/cohort1/*.png')
### Added for including a tsv files
parser.add_argument('--processing_order', type=str, default='from start', help='Choose between , from start, from end, random. To deal with the images to process')
parser.add_argument('--wsi_file_list', type=str, default='./testimage/wsi_file_list.tsv', help='path to the TSV file containing the list of wsi file paths and names')

args = parser.parse_args()

if not os.path.isdir(args.output_path):
    os.mkdir(args.output_path)

if __name__ == "__main__":
    start = timer()
    results = main(args)
    end = timer()
    print('Script Time: %f secons' % (end - start))