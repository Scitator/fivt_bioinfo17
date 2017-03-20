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


def path_reconstruction(k, d, kmers):
    path1 = [x[0] for x in kmers]
    path2 = [x[1] for x in kmers]
    result1 = "".join([x[0] for x in path1] + [path1[-1][1:]])
    result2 = "".join([x[0] for x in path2] + [path2[-1][1:]])
    result = result1[:k+d] + result2
    return result


def dna_reconstruction(k, d, kmers):
    graph = build_graph(kmers)
    path = find_eulerian_tour(graph)
    result = path_reconstruction(k, d, path)
    return "".join(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem13_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    k, d = tuple(map(int, data[0].split()))
    kmers = list(map(lambda x: tuple(x.replace("\n", "").split("|")), data[1:]))
    print(dna_reconstruction(k, d, kmers))


if __name__ == '__main__':
    main()
