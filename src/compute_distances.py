import argparse
import numpy as np
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encodings', type=str, required=True)
    parser.add_argument('--embeddings', type=str, required=True)
    parser.add_argument('--chrm', type=str, required=True)
    parser.add_argument('--out', type=str, required=True)

    return parser.parse_args()

def compute_segment_distances(gt_file,
                            emb_file):
    '''
    compute pairwise sample encoding and embedding distances
    '''
    try:
        gtf = open(gt_file, 'r')
        gt_dict = dict()
        for enc in gtf:
            L_enc = enc.strip().split()
            sample_ID = L_enc[0]
            encoding = [int(i) for i in L_enc[1:]]
            gt_dict[sample_ID] = encoding
        gtf.close()

        embf = open(emb_file, 'r')
        emb_dict = dict()
        for emb in embf:
            L_emb = emb.strip().split()
            sample_ID = L_emb[0]
            print('emb: ' + sample_ID)
            embedding = [float(i) for i in L_emb[1:]]
            emb_dict[sample_ID] = embedding
        embf.close()

    except:
        print('cannot open...' + gt_file + ' or ' + emb_file)
    
    return gt_dict, emb_dict

def compute_sample_dist(gt_dict,
                        emb_dict, 
                        out_file):
    '''
    compute distances and write out to file
    '''
    df = open(out_file, 'w')
    for sample_A in gt_dict:
        for sample_B in gt_dict:
            if sample_A == sample_B:
                continue
            else:
                enc_dist = euclidean_dist(gt_dict[sample_A], gt_dict[sample_B])
                emb_dist = euclidean_dist(emb_dict[sample_A], emb_dict[sample_B])
                df.write(sample_A + '\t' + sample_B + '\t' + str(enc_dist) + '\t' + str(emb_dist) + '\n')
                
    df.close()

def euclidean_dist(vector1, vector2): 
    '''
    compute euclidean distance between two vectors
    '''
    return np.linalg.norm(np.array(vector1) - np.array(vector2))
    
def main():
    args = get_args()
    
    encodings_dir = args.encodings
    embeddings_dir = args.embeddings
    chrm = args.chrm
    out_dir = args.out

    segments = [81]
    for seg in segments:
        seg = str(seg)
        gt_file = encodings_dir + 'chrm'+chrm+'.segment'+seg+'.gt'
        emb_file = embeddings_dir + 'chrm'+chrm+'.segment'+seg+'.emb'
        distance_file = out_dir + 'chrm'+chrm+'.segment'+seg+'.dist'
    
        print(gt_file)
        gt_dict, emb_dict = compute_segment_distances(gt_file,
                                                        emb_file)

        compute_sample_dist(gt_dict, emb_dict,
                            distance_file)

if __name__ == '__main__':
    main()
