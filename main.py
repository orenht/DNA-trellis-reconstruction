from consts import *
import encoder
import trellis_graph
import algorithms.multi_trace
import algorithms.trellis_bma
import random

import logging
import argparse
import time

def read_centers_and_clusters():
    centers = []
    with open("Centers.txt") as f:
        centers = f.read().splitlines()

    with open("Clusters.txt") as f:
        clusters = [[c for c in cluster.splitlines() if c != ""]
                    for cluster in f.read().split("===============================")
                    if cluster != ""]

    return centers, clusters

if __name__ == '__main__':
    logging.basicConfig(level=LOGGING_LEVEL)
    parser = argparse.ArgumentParser(description='Reconstruct DNA string from noisy traces using trellis.')
    parser.add_argument("--trace-num", type=int, help="Number of traces for each string")
    parser.add_argument("--from-idx", type=int, help="index of clusters in dataset to start from")
    parser.add_argument("--to-idx", type=int, help="index of clusters in dataset to finish (non inclusive)")
    parser.add_argument("--results-file", type=str, help="write the reconstruction results to this file")
    args = parser.parse_args()

    if not args.results_file:
        current_datetime = time.strftime("%Y%m%d-%H_%M_%S")
        args.results_file = f"results_{args.trace_num}_traces_{current_datetime}.txt"

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
