import os
import tensorflow as tf2
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import pylib.utils as utils
from PIL import Image
import numpy as np
# use only if results are to be saved as h5 files
#import h5py
from timeit import default_timer as timer
import argparse

def main(args):
    tf.flags.DEFINE_string('gpu_index', '0', 'gpu index, default: 0')
    FLAGS = tf.flags.FLAGS
    flags = FLAGS
    #os.environ['CUDA_VISIBLE_DEVICES'] = FLAGS.gpu_index
    run_config = tf.ConfigProto()
    #run_config.gpu_options.allow_growth = True
    sess = tf.Session(config=run_config)
    model = utils.net_arch(sess, flags, (2000,2000))
    saver = tf.train.Saver()
    sess.run(tf.global_variables_initializer())
    model_out_dir = "./models/model_epistroma"  
    if utils.load_model(sess, saver, model_out_dir):
        best_auc_sum = sess.run(model.best_auc_sum)
        print('====================================================')
        print(' Best auc_sum: {:.3}'.format(best_auc_sum))
        print('=============================================')    
        print(' [*] Load Success!\n')
    names = utils.all_files_under(args.input_path,'.png')
    for name in names:
        print(name)
        savingName = args.output_path + os.path.basename(name)[:-4]
        if os.path.isfile(savingName) or os.path.getsize(name)>>20 < 2:
            print('size less than 2MB or file existed , continue to next')
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
        image.save(savingName, format='PNG')
        # save as hdf5 file (.mat)
        #hdf5 = h5py.File(savingName,mode='w')
        #hdf5.create_dataset('mask', data = samples)
        #hdf5.close()

### Working paths
parser = argparse.ArgumentParser(description='configuration for running the Epithelium-Stroma segmentation')
parser.add_argument('--input_path', type=str, default='./testimage/patch/', help='input directory')
parser.add_argument('--output_path', type=str, default='./testimage/epimask/', help='output directory for epistroma mask')

args = parser.parse_args()

if not os.path.isdir(args.output_path):
    os.mkdir(args.output_path)

if __name__ == "__main__":
    start = timer()
    results = main(args)
    end = timer()
    print('Script Time: %f secons' % (end - start))