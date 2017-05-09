import click
import numpy as np
import collections


def two_break_distance(P, Q):
    # build the graph
    graph = collections.defaultdict(list)
    for perm_cycle in P + Q:
        cycle_len = len(perm_cycle)
        for i in range(cycle_len):
            graph[perm_cycle[i]].append(-1 * perm_cycle[(i + 1) % cycle_len])
            graph[-1 * perm_cycle[(i + 1) % cycle_len]].append(perm_cycle[i])

    # count it all!
    components_num = 0
    curr_keys = set(graph.keys())
    while len(curr_keys) > 0:
        components_num += 1
        queue_keys = [curr_keys.pop()]
        while queue_keys:
            key = queue_keys.pop()
            queue_keys += list(filter(lambda node: node in curr_keys, graph.get(key, [])))
            curr_keys -= set(queue_keys)

    result = sum(map(len, P)) - components_num
    return result


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem28_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()[:2]
    P, Q = [line.strip().replace("\n", "").replace("+", "")[1:-1].split(")(") for line in data]
    P = list(list(map(int, line.split())) for line in P)
    Q = list(list(map(int, line.split())) for line in Q)
    result = two_break_distance(P, Q)
    print(result)


if __name__ == '__main__':
    main()
