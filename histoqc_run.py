import subprocess
from timeit import default_timer as timer
import argparse

def main(args):
    path_in = args.input_path #'/home/ubuntu/wsi/bms/'
    path_out = args. output_path #'/mnt/data01/histoqc/bms/'
    path_mask = args.mask_path #'/mnt/data01/histoqc/bms/mask_use/'
    ext = 'svs'
    subprocess.run(["python","-m","histoqc","-c","v2.1","-n","2",path_in+'*.'+ext,"--outdir",path_out,"--maskuse_outdir",path_mask], cwd="./nucleus_segmentation_wsi_hovernet/HistoQC/")

### Working paths
parser = argparse.ArgumentParser(description='configuration for running the histoQC')
parser.add_argument('--input_path', type=str, default='./testimage/wsi/', help='input directory of the wsi images')
parser.add_argument('--output_path', type=str, default='./testimage/nucleimask/', help='output directory to save the results of histoQC')
parser.add_argument('--mask_path', type=str, default='./testimage/nucleimask/', help='output directory to save (move) the tissue masks of the wsi')

args = parser.parse_args()

if not os.path.isdir(args.output_path):
    os.mkdir(args.output_path)

if __name__ == "__main__":
    start = timer()
    results = main(args)
    end = timer()
    print('Script Time: %f secons' % (end - start))