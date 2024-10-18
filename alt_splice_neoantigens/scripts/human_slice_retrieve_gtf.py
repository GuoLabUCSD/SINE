import os
import argparse
import numpy as np
import re

def slice_gtf(junction, gtf_path, save_dir='/tmp'):


        if os.path.isfile(gtf_path):

                # no colon in junction name for easier processing
                junction = junction.replace(':','-')
                savepath = os.path.join(save_dir, '{}.gtf'.format(junction))

                if os.path.isfile(savepath):
                        return savepath


                # make junction bed
                bed_path = os.path.join(save_dir, '{}.bed'.format(junction))
                np.savetxt(bed_path, [junction.replace('-','\t')], delimiter='\n', fmt='%s')

                os.system('bedtools intersect -wa -a {} -b {} > {}'.format(gtf_path, bed_path, savepath))
                return savepath

        else:
                print('{} does not exist.'.format(gtf_path))
                return None

def retrieve(new_gtf_path):
        if os.path.isfile(new_gtf_path):
                enst_list = []
                with open(new_gtf_path, 'r') as f:
                    for line in f.readlines():
                        enst_list.extend(re.findall('ENST[0-9]*',line))
                return list(set(enst_list))
        else:
                print('{} does not exist.'.format(new_gtf_path))
                return None
        

if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        path_args = parser.add_argument_group("Input/output options:")

        path_args.add_argument('--gtf_path', type=str, required=True) 
        path_args.add_argument('--junction', type=str, required=True) 
        path_args.add_argument('--save_dir', type=str, default='/tmp/') 
        args = parser.parse_args()

        savepath = slice(args.gtf_path, args.junction, args.save_dir)
        retrieve(savepath)
