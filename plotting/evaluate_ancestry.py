import argparse
import ancestry_helpers as ah
import violin_plot as vp

def parse_args():
    parser=argparse.ArgumentParser(description="evaluate ancestry for 1KG")
    parser.add_argument("--pop", help="path to population file")
    parser.add_argument("--knn", help="path to knn results file")
    parser.add_argument("--png", help="directory to write png files")
    return parser.parse_args()

def main():
    args=parse_args()
    population_file = args.pop
    knn_results = args.knn
    png_dir = args.png

    # get population mappings
    sample_subpopulations, sub_to_super, super_to_sub = ah.get_population_maps(population_file)

    # get KNN results
    sample_knn = ah.get_knn_results(knn_results)

    # plot violin plots for each subpopulation
    vp.plot_subpopulation(sample_subpopulations,
                          sub_to_super,
                          super_to_sub,
                          sample_knn,
                          png_dir)


if __name__ == '__main__':
    main()