import pandas as pd
import click
import collections


def kmer_suffix(kmer):
    return kmer[1:]


def kmer_prefix(kmer):
    return kmer[:-1]


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def build_graph(kmers):
    graph = collections.defaultdict(list)
    for kmer in kmers:
        prefix = kmer_prefix(kmer)
        suffix = kmer_suffix(kmer)
        graph[prefix].append(suffix)
    return graph


def find_start_vertex(graph):
    counter = collections.defaultdict(lambda: 0)
    for key, value in graph.items():
        counter[key] += 0
        if len(value) == 0:
            return key
        for node in value:
            counter[node] += 1
    counter_sort = sorted(counter.items(), key=lambda x: x[1])
    return counter_sort[0][0]


def find_eulerian_tour(graph):
    """
    stack St;
    в St кладём любую вершину (стартовая вершина);
    пока St не пустой
        пусть V - значение на вершине St;
        если степень(V) = 0, то
            добавляем V к ответу;
            снимаем V с вершины St;
        иначе
            находим любое ребро, выходящее из V;
            удаляем его из графа;
            второй конец этого ребра кладём в St;
    """
    ans = []
    stack = [find_start_vertex(graph)]
    while stack:
        curr_v = stack[-1]
        if len(graph[curr_v]) == 0:
            ans.append(curr_v)
            stack.pop()
        else:
            next_v = graph[curr_v].pop()
            stack.append(next_v)
    return list(reversed(ans))


def dna_reconstruction(k, dna):
    kmers = [x for x in chunks(dna, k)]
    graph = build_graph(kmers)
    path = find_eulerian_tour(graph)
    result = [x[0] for x in path] + [path[-1][1:]]
    return "".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem11_input.tsv")
def main(fin):
    df = pd.read_csv(fin, sep="\t")
    assert all(x in df.columns.values.tolist() for x in ["k", "dna"])
    for i, row in df.iterrows():
        print(dna_reconstruction(row["k"], row["dna"]))


if __name__ == '__main__':
    main()
