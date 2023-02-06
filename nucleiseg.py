import os
import tensorflow as tf2
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import pylib.utils as utils
from PIL import Image
import numpy as np
from timeit import default_timer as timer
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def main(args):
    tf.flags.DEFINE_string('gpu_index', '0', 'gpu index, default: 0')
    FLAGS = tf.flags.FLAGS
    flags = FLAGS
    run_config = tf.ConfigProto()
    sess = tf.Session(config=run_config)
    model = utils.net_arch(sess, flags, (2000,2000))
    saver = tf.train.Saver()
    sess.run(tf.global_variables_initializer())
    model_out_dir = "./models/model_nuclei"  
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
            savingName = args.output_path + os.path.basename(name)[:-4]
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
            samples = np.squeeze(samples*255).astype(np.uint8)   
            # save as png file
            image = Image.fromarray(samples)
            image.save(savingName+'.png', format='PNG')
    elif args.cohort_type == 'folder':
        # list of the folders (cohorts)
        folder_list = utils.all_files_under(args.input_wsi)
        for folder in folder_list:
            folder_dir = os.path.splitext(folder)[0]
            folder_name = os.path.splitext(os.path.basename(folder_dir))[0]
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
                #img_norm = samples
                #img_norm = (samples - np.amin(samples))*255.0/(np.amax(samples) - np.amin(samples))
                #img_norm = np.uint8(img_norm)
                #print(samples.shape)
                #print(np.amin(samples))
                #print(np.amax(samples))
                # save as png file
                image = Image.fromarray(samples)
                image.save(savingName, format='PNG')
                #print(np.amin(img_norm))
                #rgb_img = np.dstack((img_norm,img_norm,img_norm))
                #plt.imsave(savingName, rgb_img, cmap=cm.gray)    

### Working paths
parser = argparse.ArgumentParser(description='configuration for running the nuclei segmentation')
parser.add_argument('--input_wsi', type=str, default='./testimage/wsi/', help='input directory of the wsi (folder option)')
parser.add_argument('--input_path', type=str, default='./testimage/patch/', help='input directory of the images (patch option)')
parser.add_argument('--output_path', type=str, default='./testimage/nucleimask/', help='output directory for nuclei mask')
parser.add_argument('--cohort_type', type=str, default='patch', help='Choose between patch or folder. running on a single folder (/cohort/*.png) of png images or a directory of folders (e.g. /cohorts/cohort1/*.png')

args = parser.parse_args()

if not os.path.isdir(args.output_path):
    os.mkdir(args.output_path)

if __name__ == "__main__":
    start = timer()
    results = main(args)
    end = timer()
    print('Script Time: %f secons' % (end - start))