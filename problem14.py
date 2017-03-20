import numpy as np
import click
import collections


def kmer_suffix(kmer):
    if isinstance(kmer, (tuple)):
        return tuple([x[1:] for x in kmer])
    elif isinstance(kmer, str):
        return kmer[1:]
    else:
        raise Exception()


def kmer_prefix(kmer):
    if isinstance(kmer, (list, tuple)):
        return tuple([x[:-1] for x in kmer])
    elif isinstance(kmer, str):
        return kmer[:-1]
    else:
        raise Exception()


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


def count_in_out(graph):
    counter = collections.defaultdict(lambda: (0, 0))
    for key, value in graph.items():
        counter[key] = (counter[key][0], counter[key][1] + len(value))
        for node in value:
            counter[node] = (counter[node][0] + 1, counter[node][1])
    return counter


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


def node_contig(graph, node, in_out):
    contigs = []
    for value in graph[node]:
        path = [node, value]
        in_, out_ = in_out[value]
        while in_ == 1 and out_ == 1:
            node = value
            value = graph[node][0]
            path.append(value)
            in_, out_ = in_out[value]
        contigs.append(path)
    return contigs


def graph2contig(graph):
    in_out = count_in_out(graph)
    contigs = []
    for key, (in_, out_) in in_out.items():
        if out_ > 0 and not (out_ == 1 and in_ == 1):
            contigs.extend(node_contig(graph, key, in_out))
    return contigs


def contig_generation(kmers):
    graph = build_graph(kmers)
    path = graph2contig(graph)
    contig = list(map(
        lambda ys: "".join([x[0] for x in ys] + [ys[-1][1:]]), path))
    return "\t".join(contig)
    

@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem14_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    kmers = list(map(lambda x: x.replace("\n", "").strip(), data))
    print(contig_generation(kmers))


if __name__ == '__main__':
    main()
