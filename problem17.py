import click
char2mass = dict(G=57, A=71, S=87, P=97, V=99, T=101, C=103, I=113, L=113, N=114, D=115, K=128,
                 Q=128, E=129, M=131, H=137, F=147, R=156, Y=163, W=186)


def factorial(n):
    num = 1
    while n >= 1:
        num *= n
        n -= 1
    return num


def permutation_count(str):
    str = sorted(str)
    result = factorial(len(str))
    curr = str[0]
    count = 0

    for c in str:
        if c != curr:
            result /= factorial(count)
            curr = c
            count = 0

        count += 1

    result /= factorial(count)
    return result


def problem2solve(mass):
    result = 0

    del char2mass['Q']
    del char2mass['L']

    keys = list(char2mass.keys())
    curr_state = [(i, key, value) for i, (key, value) in enumerate(char2mass.items())]
    next_state = []

    while curr_state:
        for (i, key, value) in curr_state:
            for idx in range(i, len(keys)):
                next_key = keys[idx]
                curr_mass = value + char2mass[next_key]
                curr_key = key + next_key

                if curr_mass > mass:
                    continue
                if curr_mass == mass:
                    result += permutation_count(curr_key)
                    continue

                next_state.append((idx, curr_key, curr_mass))

        curr_state = next_state
        next_state = []

    return int(result)


@click.command()
@click.option(
    "--fin",
    type=str,
    default="problem17_input.txt")
def main(fin):
    with open(fin) as fin:
        data = fin.readlines()
    m = tuple(map(int, data[0].split()))[0]
    print(problem2solve(m))


if __name__ == '__main__':
    main()
