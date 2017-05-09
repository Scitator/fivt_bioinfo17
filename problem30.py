import click
import numpy as np
import collections


def chromosome2cycle(chromosome):
    nodes = []
    weight = [-1, 1]
    bias = [[0, -1], [-1, 0]]
    for i in chromosome:
        nodes.append(weight[i > 0] * 2 * i + bias[i > 0][0])
        nodes.append(weight[i > 0] * 2 * i + bias[i > 0][1])

    return nodes


def color_edges(genome):
    edges = []
    for chrome in genome:
        cycle = chromosome2cycle(chrome)
        for j in range(len(cycle) // 2):
            edge = (cycle[1 + 2 * j], cycle[(2 + 2 * j) % len(cycle)])
            edges.append(edge)
    return edges


def graph2genome(graph):
    adj = np.zeros(len(graph) * 2, dtype=np.int)
    for edge in graph:
        adj[edge[0] - 1] = edge[1] - 1
        adj[edge[1] - 1] = edge[0] - 1

    genome = []
    closed = []
    for edge in graph:
        start_v = edge[0]
        if start_v in closed:
            continue
        closed.append(start_v)

        closing = start_v - 1 if start_v % 2 == 0 else start_v + 1

        block = []
        while True:
            next_vert = start_v / 2 if start_v % 2 == 0 else -(start_v + 1) / 2
            block.append(next_vert)

            finish_v = adj[start_v - 1] + 1
            closed.append(finish_v)

            if finish_v == closing:
                genome.append(block)
                break

            start_v = finish_v - 1 if finish_v % 2 == 0 else finish_v + 1
            closed.append(start_v)

    return genome


def two_break_on_genome_graph(graph, i, j, k, l):
    to_del = ((i, j), (j, i), (k, l), (l, k))
    break_genome = [t for t in graph if t not in to_del]
    break_genome.append((i, k))
    break_genome.append((j, l))

    return break_genome


def two_break_on_genome(genome, i, j, k, l):
    graph = color_edges(genome)
    graph = two_break_on_genome_graph(graph, i, j, k, l)
    genome = graph2genome(graph)
    return genome


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem30_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()[:2]
    P = data[0].strip().replace("\n", "").replace("+", "")[1:-1].split(")(")
    Q = list(map(int, data[1].split(", ")))
    P = list(list(map(int, line.split())) for line in P)
    result = two_break_on_genome(P, *Q)

    pretty_value_fn = lambda value: ['-', '+'][value > 0] + str(abs(int(value)))
    result = " ".join(
        ['(' + ' '.join([pretty_value_fn(value) for value in block]) + ')'
         for block in result])
    # result = list(map(lambda block: "(", result))
    print(result)


if __name__ == '__main__':
    main()
