import argparse
import utils

def get_args():
    parser = argparse.ArgumentParser(description='Suffix Trie')

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

def build_suffix_trie(s):
    root = {}
    for i in range(len(s)):
        current_node = root
        for char in s[i:]:
            if char not in current_node:
                current_node[char] = {}
            current_node = current_node[char]
        current_node['$'] = True  # Mark the end of a suffix
        if i % 1000 == 0:  # Print progress every 1000 iterations
            print(f"Inserted suffix starting at index {i}")
    return root

def search_trie(trie, pattern):
    current_node = trie
    match_length = 0
    for char in pattern:
        if char not in current_node:
            break
        current_node = current_node[char]
        match_length += 1
    return match_length

def main():
    args = get_args()

    T = None

    if args.string:
        T = args.string
    elif args.reference:
        reference = utils.read_fasta(args.reference)
        T = reference[0][1]

    print(f"Building suffix trie for sequence of length {len(T)}")
    trie = build_suffix_trie(T)
    print("Suffix trie built successfully")

    if args.query:
        for query in args.query:
            match_len = search_trie(trie, query)
            print(f'{query} : {match_len}')

if __name__ == '__main__':
    main()
