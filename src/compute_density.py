import argparse
import numpy as np
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encodings', type=str, required=True)
    parser.add_argument('--chrm', type=str, required=True)
    parser.add_argument('--out', type=str, required=True)

    return parser.parse_args()

def compute_segment_density(segment_file, density_file):
    '''
    write sample and segment density out to file
    '''
    try:
        with open(segment_file, 'r') as sf:
            df = open(density_file, 'w')
            for line in sf:
                L = line.strip().split()
                sample_ID = L[0]
                encoding = [int(i) for i in L[1:]]
                density = np.sum(encoding)
                df.write(f'{sample_ID}\t{density}\n')
        df.close()
    except:
        print('cannot open...' + segment_file)
        

def main():
    args = get_args()
    
    encodings_dir = args.encodings
    chrm = args.chrm
    out_dir = args.out

    segments = range(0,170)
    for seg in segments:
        seg = str(seg)
        gt_file = encodings_dir + 'chrm'+chrm+'.segment'+seg+'.gt'
        density_file = out_dir + 'chrm'+chrm+'.segment'+seg+'.density'
    
        print(gt_file)
        compute_segment_density(gt_file, density_file)    

if __name__ == '__main__':
    main()
