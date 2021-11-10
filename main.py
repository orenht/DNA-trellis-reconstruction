import sys

from consts import *
import encoder
import trellis_graph
import algorithms.multi_trace
import algorithms.trellis_bma
import random

import logging
import argparse
import time
from typing import NamedTuple
from collections import Counter
import matplotlib.pyplot as plt


class ResultEntry(NamedTuple):
    index: int
    original: str
    estimate: str
    hamming: int
    levenstein: int


def read_centers_and_clusters():
    centers = []
    with open("Centers.txt") as f:
        centers = f.read().splitlines()

    with open("Clusters.txt") as f:
        clusters = [[c for c in cluster.splitlines() if c != ""]
                    for cluster in f.read().split("===============================")
                    if cluster != ""]

    return centers, clusters


def print_results(results):
    for entry in results:
        print(f"{entry.index}:")
        print(f"hamming: {entry.hamming}, levenstein: {entry.levenstein}")
        print(f"original: {entry.original}")
        print(f"estimate: {entry.estimate}")
        print("----")


if __name__ == '__main__':
    logging.basicConfig(level=LOGGING_LEVEL)
    parser = argparse.ArgumentParser(description='Reconstruct DNA string from noisy traces using trellis.')
    main_options = parser.add_mutually_exclusive_group()
    main_options.add_argument("-p", "--parse", action="store_true", help="parse the input file instead of performing reconstruction")
    main_options.add_argument("-r", "--reconstruct", action="store_true", help="reconstruct DNA string from traces")

    parser.add_argument("--trace-num", type=int, help="Number of traces to use for each cluster. "
                                                      "if a cluster has less than this number, it will be ignored")
    parser.add_argument("--from-idx", type=int, help="index of clusters in dataset to start from")
    parser.add_argument("--to-idx", type=int, help="index of clusters in dataset to finish (non inclusive)")
    parser.add_argument("--results-file", type=str, help="write the reconstruction results to this file")
    parser.add_argument("--input-results-file", type=str, help="when --parse is used, parse this results file")
    parser.add_argument("-wh", "--worst-n-hamming", type=int, help="output the worst N reconstructions by hamming distance")
    parser.add_argument("-wl", "--worst-n-levenstein", type=int, help="output the worst N reconstructions by levenstein distance")
    parser.add_argument("-eh", "--error-histogram", action="store_true", help="plot hamming and levenstein error histograms")
    args = parser.parse_args()

    if not args.results_file:
        current_datetime = time.strftime("%Y%m%d-%H_%M_%S")
        args.results_file = f"results_{args.trace_num}_traces_{current_datetime}.txt"

    if args.parse:
        if not args.input_results_file:
            print("missing input results file!")
            sys.exit(1)
        with open(args.input_results_file) as f:
            data = f.readlines()
        stripped_data = [[s.strip() for s in line.split(",")] for line in data]
        results = [ResultEntry(int(r[0]), r[4], r[3], int(r[1]), int(r[2])) for r in stripped_data]

        # general info
        hammings = [r.hamming for r in results]
        levensteins = [r.levenstein for r in results]
        avg_hamming = sum(hammings) / len(hammings)
        avg_levenstein = sum(levensteins) / len(levensteins)
        normalized_hamming = sum([r.hamming/len(r.original) for r in results]) / len(results)
        normalized_levenstein = sum([r.levenstein / len(r.original) for r in results]) / len(results)

        print(f"read {len(results)} samples:")
        print(f"Normalized hamming: {normalized_hamming}, avg hamming: {avg_hamming}")
        print(f"Normalized levenstein: {normalized_levenstein}, avg levenstein: {avg_levenstein}")

        if args.error_histogram:
            hammings_counter = Counter(hammings)
            levensteins_counter = Counter(levensteins)

            fig, (sub1, sub2) = plt.subplots(2)

            sub1.hist(hammings, bins=max(hammings))
            sub1.set_ylabel("count")
            sub1.set_xlabel("hamming distance")
            sub2.hist(levensteins, bins=max(levensteins))
            sub2.set_ylabel("count")
            sub2.set_xlabel("levenstein distance")
            #plt.hist(levensteins)
            #plt.bar(levensteins_counter.keys(), levensteins_counter.values())
            #fig.show()
            plt.tight_layout()
            plt.show()
        if args.worst_n_hamming:
            results.sort(key=lambda res: res.hamming, reverse=True)
            print(f"worst {args.worst_n_hamming} cases by hamming distance:")
            print_results(results[:args.worst_n_hamming])
        if args.worst_n_levenstein:
            results.sort(key=lambda res: res.levenstein, reverse=True)
            print(f"worst {args.worst_n_levenstein} cases by levenstein distance:")
            print_results(results[:args.worst_n_levenstein])
    elif args.reconstruct:
        if USE_NANOPORE_DATA_FROM_FILE:
            centers, clusters = read_centers_and_clusters()
            results = []
            total_processed = 0
            with open(args.results_file, "w") as f:
                for idx, (original, traces) in enumerate(zip(centers[args.from_idx: args.to_idx], clusters[args.from_idx: args.to_idx])):
                    if len(traces) < args.trace_num:
                        logging.info(f"skipping cluster {args.from_idx+idx}, has {len(traces)} traces and {args.trace_num} required")
                        continue
                    chosen_traces = random.sample(traces, args.trace_num)
                    logging.info(f"original: {original}")
                    logging.info("traces:")
                    logging.info(",\n".join(chosen_traces))
                    total_processed += 1

                    result = algorithms.trellis_bma.compute_trellis_bma_estimation(chosen_traces, original)
                    results.append(result)
                    original, estimate, hamm, levenstein = result
                    f.write(f"{args.from_idx+idx}, {hamm}, {levenstein}, {estimate}, {original}\n")
                    f.flush()

        else:
            if USE_CUSTOM_TRACES:
                original = CUSTOM_ORIGINAL
                traces = CUSTOM_TRACES
            else:
                original = ORIGINAL
                traces = [encoder.create_noisy_trace(ORIGINAL) for i in range(TRACE_NUM)]
            logging.info(traces)
            #trellis_graph = trellis_graph.build_new(original, traces)
            #algorithms.multi_trace.compute_marginal_prob_for_each_vertex(trellis_graph, traces, original)
            algorithms.trellis_bma.compute_trellis_bma_estimation(traces, original)
    else:
        print("must choose args or reconstruct!")
