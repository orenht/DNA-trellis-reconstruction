import trellis_graph
import algorithms.common
from consts import *

from typing import NamedTuple, List, Dict, DefaultDict
import networkx as nx
import numpy as np
from collections import defaultdict
import math

beta_b = 0.1
beta_i = 0
beta_e = 1


class TrellisMetadata(NamedTuple):
    trellis: nx.DiGraph
    topological_ordering: List[trellis_graph.TrellisVertex]
    F_values: Dict[trellis_graph.TrellisVertex, float]
    B_values: Dict[trellis_graph.TrellisVertex, float]
    Sm_per_stage_per_symbol: Dict[int, Dict[str, List[trellis_graph.TrellisVertex]]]
    vertices_by_stage_by_substage_sorted_topologically: Dict[int, Dict[int, List[trellis_graph.TrellisVertex]]]
    Vk_estimations_per_stage_per_symbol: DefaultDict[int, Dict[str, float]] = {}


def compute_trellis_bma_estimation(traces, original):
    if ESTIMATE_SECOND_HALF_REVERSED:
        first_half_estimates = build_trellis_and_estimate(traces, original, len(original) // 2)
        reversed_traces = [t[::-1] for t in traces]
        reversed_original = original[::-1]
        second_half_length = len(original) - (len(original)//2)
        second_half_estimates = list(reversed(build_trellis_and_estimate(reversed_traces, reversed_original, second_half_length)))
        final_estimate = first_half_estimates + second_half_estimates
    else:
        final_estimate = build_trellis_and_estimate(traces, original, len(original))

    #for stage in reversed(range(len(original) // 2, len(original))):
    #    sync_probabilities_and_estimate_stage(stage, trellises_metadata, traces,
    #                                          V_per_index_per_symbol[stage], is_first_half=False)
    #
    #    second_half_estimates.insert(0, max(Alphabet, key=lambda c: V_per_index_per_symbol[stage][c]))

    final_estimate_str = "".join(final_estimate)
    print("estimate:")
    print(final_estimate_str)
    print("original:")
    print(original)
    hamm = algorithms.common.hamming_distance(original, final_estimate_str)
    print(f"hamming distance from original: {hamm}")
    levenstein = algorithms.common.levenshtein_distance(original, final_estimate_str)
    print(f"levenstein distance: {levenstein}")
    print(f"traces: {traces}")

    return original, final_estimate_str, hamm, levenstein


def build_trellis_and_estimate(traces, original, estimation_len):
    trellises_metadata = []

    for trace in traces:
        trellis = trellis_graph.build_new(original, [trace])
        topological_ordering = list(nx.topological_sort(trellis))
        F_values = algorithms.common.compute_Fs_for_all_nodes(trellis, topological_ordering)
        B_values = algorithms.common.compute_Bs_for_all_nodes(trellis, topological_ordering, [trace])
        Sm_per_stage = algorithms.common.compute_Sm_foreach_stage_and_symbol(trellis, [trace])
        vertices_by_stage_by_substage = algorithms.common.get_vertices_by_stage_by_substage_sorted_topologically(
            topological_ordering)
        vk_estimations = defaultdict(dict)
        trellises_metadata.append(TrellisMetadata(trellis, topological_ordering,
                                                  F_values, B_values, Sm_per_stage,
                                                  vertices_by_stage_by_substage, vk_estimations))

    # Vk_per_trace_per_stage = defaultdict(lambda: defaultdict(float))
    V_per_index_per_symbol = defaultdict(dict)
    estimates = []

    # first half of stages
    # for stage in range(len(original) // 2):
    for stage in range(estimation_len):
        # # compute V^k for all trellises
        # for trace_idx, trellis_metadata in enumerate(trellises_metadata):
        #     for symbol in Alphabet:
        #         stage_final_vertices = trellis_metadata.Sm_per_stage_per_symbol[stage][symbol]
        #         #for v in stage_final_vertices:
        #         #    print(f"F[{v}]={trellis_metadata.F_values[v]}")
        #         if USE_LOG_PROB:
        #             #values = [(v, trellis_metadata.F_values[v] + beta_b * trellis_metadata.B_values[v])
        #             #          if beta_b != 0 or trellis_metadata.B_values[v] != -math.inf
        #             #          else (v, -math.inf)
        #             #          for v in stage_final_vertices]
        #             f_values = [trellis_metadata.F_values[v] for v in stage_final_vertices]
        #             #b_values = [beta_b * trellis_metadata.B_values[v] for v in stage_final_vertices]
        #             print(f_values)
        #             #print(b_values)
        #             #print(values)
        #             #estimated_v = np.logaddexp.reduce([x[1] for x in values])
        #             estimated_v = np.logaddexp.reduce([trellis_metadata.F_values[v] + beta_b * trellis_metadata.B_values[v]
        #                                               if beta_b != 0 or trellis_metadata.B_values[v] != -math.inf
        #                                               else -math.inf
        #                                               for v in stage_final_vertices])
        #         else:
        #             estimated_v = sum(trellis_metadata.F_values[v] * (trellis_metadata.B_values[v] ** beta_b)
        #                               if not math.isnan(trellis_metadata.B_values[v] ** beta_b)
        #                               else 0
        #                               for v in stage_final_vertices)
        #         trellis_metadata.Vk_estimations_per_stage_per_symbol[stage][symbol] = estimated_v
        #         print(f"[trace{trace_idx}[{stage}]:{symbol}] = {estimated_v}")
        #
        # # Use V^k to update forward values for each trellis
        # for trace_idx, trellis_metadata in enumerate(trellises_metadata):
        #     # Compute gamma^k(m) for all trellises
        #     for symbol in Alphabet:
        #         vk_estimation = trellis_metadata.Vk_estimations_per_stage_per_symbol[stage][symbol]
        #         if USE_LOG_PROB:
        #             internal_factor = vk_estimation * beta_i
        #             external_factor = sum(trellis_j.Vk_estimations_per_stage_per_symbol[stage][symbol]
        #                                   for trellis_j in trellises_metadata)
        #             external_factor -= vk_estimation
        #             external_factor = external_factor * beta_e
        #             gamma_k = internal_factor + external_factor
        #
        #         else:
        #             internal_factor = vk_estimation ** beta_i
        #             external_factor = math.prod(trellis_j.Vk_estimations_per_stage_per_symbol[stage][symbol]
        #                                         for trellis_j in trellises_metadata)
        #             external_factor /= vk_estimation
        #             external_factor = external_factor ** beta_e
        #             gamma_k = internal_factor * external_factor
        #
        #         # update relevant F^k values
        #         print(gamma_k)
        #         for v in trellis_metadata.Sm_per_stage_per_symbol[stage][symbol]:
        #             if USE_LOG_PROB:
        #                 print(f"{v} before: {trellis_metadata.F_values[v]}")
        #                 trellis_metadata.F_values[v] += gamma_k
        #                 print(f"{v} after: {trellis_metadata.F_values[v]}")
        #             else:
        #                 trellis_metadata.F_values[v] *= gamma_k
        #
        #     # forward pass F values
        #     # update all substages in next stage, including -1
        #     for substage in [-1] + list(range(len(traces))):
        #         for v in trellis_metadata.vertices_by_stage_by_substage_sorted_topologically[stage + 1][substage]:
        #             new_f = algorithms.common.compute_F_value_for_single_node(trellis_metadata.trellis, v,
        #                                                                       trellis_metadata.F_values)
        #             trellis_metadata.F_values[v] = new_f
        #
        # # combine V^k to compute V(M_l = m)
        # for symbol in Alphabet:
        #     if USE_LOG_PROB:
        #         V_per_index_per_symbol[stage][symbol] = sum(metadata.Vk_estimations_per_stage_per_symbol[stage][symbol]
        #                                                     for metadata in trellises_metadata)
        #     else:
        #         V_per_index_per_symbol[stage][symbol] = math.prod(metadata.Vk_estimations_per_stage_per_symbol[stage][symbol]
        #                                                           for metadata in trellises_metadata)
        print(f"estimating stage: {stage}")
        sync_probabilities_and_estimate_stage(stage, trellises_metadata, traces,
                                              V_per_index_per_symbol[stage], is_first_half=True)

        estimates.append(max(Alphabet, key=lambda c: V_per_index_per_symbol[stage][c]))

    return estimates


def sync_probabilities_and_estimate_stage(stage, trellises_metadata, traces, V_per_symbol, is_first_half=True):
    # compute V^k for all trellises
    for trace_idx, trellis_metadata in enumerate(trellises_metadata):
        compute_vk_for_trellis(stage, trace_idx, trellis_metadata)
        # if is_first_half:
        #     compute_vk_for_trellis(stage, trace_idx, trellis_metadata,
        #                            trellis_metadata.F_values, trellis_metadata.B_values)
        # else:
        #     compute_vk_for_trellis(stage, trace_idx, trellis_metadata,
        #                            trellis_metadata.B_values, trellis_metadata.F_values)

    # Use V^k to update forward values for each trellis
    for trace_idx, trellis_metadata in enumerate(trellises_metadata):
        # Compute gamma^k(m) for all trellises
        symbol_gammas = {c: compute_gamma_coeff_for_trellis_and_symbol(stage, c, trellis_metadata, trellises_metadata)
                         for c in Alphabet}
        # normalize gammas to avoid drift
        if USE_LOG_PROB:
            # normalize by subtracting most dominant probability
            factor = max(symbol_gammas.values())
            normalized_gammas = {c: symbol_gammas[c] - factor for c in Alphabet}
        else:
            factor = 1 / sum(symbol_gammas.values())
            normalized_gammas = {c: symbol_gammas[c] * factor for c in Alphabet}
        #if stage > 0:
            #print(normalized_gammas)
            #print(f"trace{trace_idx}[{stage}]: {symbol_gammas}")
            #print(f"trace{trace_idx}[{stage}]: {normalized_gammas}")
        for symbol in Alphabet:
            #gamma_k = compute_gamma_coeff_for_trellis_and_symbol(stage, symbol, trellis_metadata, trellises_metadata)

            # update relevant F^k values
            #print(normalized_gammas[symbol])
            for v in trellis_metadata.Sm_per_stage_per_symbol[stage][symbol]:
                if is_first_half:
                    if USE_LOG_PROB:
                        #print(f"{v} before: {trellis_metadata.F_values[v]}")
                        trellis_metadata.F_values[v] += normalized_gammas[symbol]
                        #trellis_metadata.F_values[v] += symbol_gammas[symbol]
                        #print(f"{v} after: {trellis_metadata.F_values[v]}")
                    else:
                        trellis_metadata.F_values[v] *= normalized_gammas[symbol]
                else:
                    if USE_LOG_PROB:
                        #print(f"{v} before: {trellis_metadata.B_values[v]}")
                        trellis_metadata.B_values[v] += normalized_gammas[symbol]
                        #trellis_metadata.F_values[v] += symbol_gammas[symbol]
                        #print(f"{v} after: {trellis_metadata.B_values[v]}")
                    else:
                        trellis_metadata.B_values[v] *= normalized_gammas[symbol]

        if is_first_half:
        #if True:
            # forward pass F values
            # update all substages in next stage, including -1
            for substage in [-1] + list(range(len(traces))):
                for v in trellis_metadata.vertices_by_stage_by_substage_sorted_topologically[stage + 1][substage]:
                    new_f = algorithms.common.compute_F_value_for_single_node(trellis_metadata.trellis, v,
                                                                              trellis_metadata.F_values)
                    trellis_metadata.F_values[v] = new_f
        else:
            # backward pass B values
            # update all substages in current stage, including -1, in reverse order
            for substage in reversed([-1] + list(range(len(traces)))):
                for v in trellis_metadata.vertices_by_stage_by_substage_sorted_topologically[stage][substage]:
                    new_b = algorithms.common.compute_B_value_for_single_node(trellis_metadata.trellis, v, traces,
                                                                              trellis_metadata.B_values)
                    trellis_metadata.B_values[v] = new_b

    # combine V^k to compute V(M_l = m)
    for symbol in Alphabet:
        if USE_LOG_PROB:
            V_per_symbol[symbol] = sum(metadata.Vk_estimations_per_stage_per_symbol[stage][symbol]
                                       for metadata in trellises_metadata)
        else:
            V_per_symbol[symbol] = math.prod(metadata.Vk_estimations_per_stage_per_symbol[stage][symbol]
                                             for metadata in trellises_metadata)


def compute_vk_for_trellis(stage: int, trace_idx: int, trellis_metadata: TrellisMetadata):#, dominant_prob_dict, secondary_prob_dict):
    for symbol in Alphabet:
        stage_final_vertices = trellis_metadata.Sm_per_stage_per_symbol[stage][symbol]
        # for v in stage_final_vertices:
        #    print(f"F[{v}]={dominant_prob_dict[v]}")
        if USE_LOG_PROB:
            # values = [(v, dominant_prob_dict[v] + beta_b * secondary_prob_dict[v])
            #          if beta_b != 0 or secondary_prob_dict[v] != -math.inf
            #          else (v, -math.inf)
            #          for v in stage_final_vertices]
            f_values = [trellis_metadata.F_values[v] for v in stage_final_vertices]
            b_values = [beta_b * trellis_metadata.B_values[v]
                        if beta_b != 0 or trellis_metadata.B_values[v] != -math.inf
                        else -math.inf
                        for v in stage_final_vertices]
            # print(values)
            # estimated_v = np.logaddexp.reduce([x[1] for x in values])
            if stage > 1145:
                print(f_values)
                print(b_values)

            estimated_v = np.logaddexp.reduce([trellis_metadata.F_values[v] + beta_b * trellis_metadata.B_values[v]
                                               if beta_b != 0 or trellis_metadata.B_values[v] != -math.inf
                                               else -math.inf
                                               for v in stage_final_vertices])
        else:
            estimated_v = sum(trellis_metadata.F_values[v] * (trellis_metadata.B_values[v] ** beta_b)
                              if not math.isnan(trellis_metadata.B_values[v] ** beta_b)
                              else 0
                              for v in stage_final_vertices)
        trellis_metadata.Vk_estimations_per_stage_per_symbol[stage][symbol] = estimated_v
        #if stage >= 9990:
        #    print(f"[trace{trace_idx}[{stage}]:{symbol}] = {estimated_v}")


def compute_gamma_coeff_for_trellis_and_symbol(stage, symbol, trellis_metadata, trellises_metadata):
    vk_estimation = trellis_metadata.Vk_estimations_per_stage_per_symbol[stage][symbol]
    if USE_LOG_PROB:
        internal_factor = vk_estimation * beta_i
        external_factor = sum(trellis_j.Vk_estimations_per_stage_per_symbol[stage][symbol]
                              for trellis_j in trellises_metadata)
        external_factor -= vk_estimation
        external_factor = external_factor * beta_e
        gamma = internal_factor + external_factor

    else:
        internal_factor = vk_estimation ** beta_i
        external_factor = math.prod(trellis_j.Vk_estimations_per_stage_per_symbol[stage][symbol]
                                    for trellis_j in trellises_metadata)
        external_factor /= vk_estimation
        external_factor = external_factor ** beta_e
        gamma = internal_factor * external_factor

    return gamma
