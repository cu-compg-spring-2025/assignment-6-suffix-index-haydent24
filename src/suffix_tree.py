import argparse
import utils

SUB = 0
CHILDREN = 1

def get_args():
    parser = argparse.ArgumentParser(description='Suffix Tree')

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

def add_suffix(nodes, suf):
    n = 0
    i = 0
    while i < len(suf):
        b = suf[i]
        children = nodes[n][CHILDREN]
        if b not in children:
            n2 = len(nodes)
            nodes.append([suf[i:], {}])
            nodes[n][CHILDREN][b] = n2
            return
        else:
            n2 = children[b]

        sub2 = nodes[n2][SUB]
        j = 0
        while j < len(sub2) and i + j < len(suf) and suf[i + j] == sub2[j]:
            j += 1

        if j < len(sub2):
            n3 = n2
            n2 = len(nodes)
            nodes.append([sub2[:j], {sub2[j]: n3}])
            nodes[n3][SUB] = sub2[j:]
            nodes[n][CHILDREN][b] = n2

        i += j
        n = n2

def build_suffix_tree(text):
    text += "$"

    nodes = [['', {}]]

    for i in range(len(text)):
        add_suffix(nodes, text[i:])
        if i % 1000 == 0:  # Print progress every 1000 iterations
            print(f"Added suffix starting at index {i}")

    return nodes

def search_tree(suffix_tree, P, text_length):
    n = 0
    i = 0
    while i < len(P):
        b = P[i]
        children = suffix_tree[n][CHILDREN]
        if b not in children:
            return 0, []
        n2 = children[b]
        sub2 = suffix_tree[n2][SUB]
        j = 0
        while j < len(sub2) and i + j < len(P) and P[i + j] == sub2[j]:
            j += 1
        if j < len(sub2):
            return 0, []
        i += j
        n = n2

    # Now we need to collect all suffixes in the subtree
    match_indexes = []
    
    def collect_suffixes(node, path_length):
        if not suffix_tree[node][CHILDREN]:  # Leaf node
            # Calculate the original position
            suffix_start = text_length - path_length
            match_indexes.append(suffix_start)
        else:
            for child_char, child_node in suffix_tree[node][CHILDREN].items():
                collect_suffixes(child_node, path_length + len(suffix_tree[child_node][SUB]))
    
    # Start collecting from the node where we found the pattern
    collect_suffixes(n, len(P))
    
    return len(match_indexes), match_indexes

def main():
    args = get_args()

    T = None

    if args.string:
        T = args.string
    elif args.reference:
        reference = utils.read_fasta(args.reference)
        T = reference[0][1]

    print(f"Building suffix tree for: {T[:50]}... (truncated)")
    tree = build_suffix_tree(T)
    print(f"Suffix tree built with {len(tree)} nodes")
        
    if args.query:
        for query in args.query:
            print(f"Searching for query: {query}")
            match_len, match_indexes = search_tree(tree, query, len(T))
            print(f'{query} : {match_len} matches at indexes {match_indexes}')

if __name__ == '__main__':
    main()
