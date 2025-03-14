import argparse
import utils

def get_args():
    parser = argparse.ArgumentParser(description='Suffix Array')

    parser.add_argument('--reference',
                        help='Reference sequence file',
                        type=str)

    parser.add_argument('--string',
                        help='Reference sequence',
                        type=str)

    parser.add_argument('--query',
                        help='Query sequences',
                        nargs='+',
                        type=str)

    return parser.parse_args()

def build_suffix_array(T):
    suffixes = [(T[i:], i) for i in range(len(T))]
    sorted_suffixes = sorted(suffixes)
    suffix_array = [suffix[1] for suffix in sorted_suffixes]
    return suffix_array

def longest_common_prefix(s1, s2):
    lcp_length = 0
    for c1, c2 in zip(s1, s2):
        if c1 == c2:
            lcp_length += 1
        else:
            break
    return lcp_length

def search_array(T, suffix_array, q):
    lo = 0
    hi = len(suffix_array)
    matches = []
    prefix_lengths = []
    max_lcp_length = 0

    while lo < hi:
        mid = (lo + hi) // 2
        suffix = T[suffix_array[mid]:]
        lcp_length = longest_common_prefix(suffix, q)
        max_lcp_length = max(max_lcp_length, lcp_length)
        if suffix.startswith(q):
            # Find the range of matches
            left = mid
            while left >= 0 and T[suffix_array[left]:].startswith(q):
                matches.append(suffix_array[left])
                prefix_lengths.append(len(q))
                left -= 1
            right = mid + 1
            while right < len(suffix_array) and T[suffix_array[right]:].startswith(q):
                matches.append(suffix_array[right])
                prefix_lengths.append(len(q))
                right += 1
            break
        elif suffix < q:
            lo = mid + 1
        else:
            hi = mid

    matches.sort()
    if matches:
        return len(matches), matches, prefix_lengths
    else:
        return 0, [], [max_lcp_length]

def main():
    args = get_args()

    T = None

    if args.string:
        T = args.string
    elif args.reference:
        reference = utils.read_fasta(args.reference)
        T = reference[0][1]

    array = build_suffix_array(T)

    if args.query:
        for query in args.query:
            print(f"Searching for query: {query}")
            match_count, match_indices, prefix_lengths = search_array(T, array, query)
            if match_count > 0:
                print(f'{query} found {match_count} times at indices {match_indices} with prefix lengths {prefix_lengths}')
            else:
                print(f'{query} not found, longest prefix overlap length: {prefix_lengths[0]}')

if __name__ == '__main__':
    main()
