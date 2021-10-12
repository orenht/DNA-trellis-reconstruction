from consts import *
import encoder
import trellis_graph
import algorithms.multi_trace
import algorithms.trellis_bma
import random

import logging


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
    if USE_NANOPORE_DATA_FROM_FILE:
        centers, clusters = read_centers_and_clusters()
        from_idx = 2500
        to_idx = 2600
        trace_num = 4
        results = []
        for original, traces in zip(centers[from_idx:to_idx], clusters[from_idx:to_idx]):
            chosen_traces = random.sample(traces, trace_num)
            logging.info(f"original: {original}")
            logging.info("traces:")
            logging.info(",\n".join(chosen_traces))

            results.append(algorithms.trellis_bma.compute_trellis_bma_estimation(chosen_traces, original))

        with open(f"results_{trace_num}_traces.txt", "w") as f:
            f.writelines((f"{hamm}, {leven}, {original}\n" for original, hamm, leven in results))
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
