import argparse
import time
import psutil
import matplotlib.pyplot as plt
import numpy as np
import os
import gc
import tracemalloc
from memory_profiler import memory_usage

# Import your implementations
import suffix_tree
import suffix_array
import utils

def get_args():
    parser = argparse.ArgumentParser(description='Benchmark Suffix Structures')

    parser.add_argument('--reference',
                        help='Reference sequence file',
                        type=str,
                        required=True)

    parser.add_argument('--queries',
                        help='File containing queries, one per line',
                        type=str,
                        required=True)

    parser.add_argument('--output',
                        help='Output directory for plots',
                        type=str,
                        default='benchmark_results')

    return parser.parse_args()

def read_queries(query_file):
    with open(query_file, 'r') as f:
        return [line.strip() for line in f]

def categorize_queries(queries):
    """Group queries by length"""
    query_groups = {}
    for query in queries:
        length = len(query)
        if length not in query_groups:
            query_groups[length] = []
        query_groups[length].append(query)
    return query_groups

def measure_build_time_and_memory(build_func, text):
    """Measure build time and memory usage"""
    # Force garbage collection before measurement
    gc.collect()
    
    # Start memory tracking
    tracemalloc.start()
    
    # Measure time
    start_time = time.time()
    structure = build_func(text)
    build_time = time.time() - start_time
    
    # Get memory stats
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    # Convert to MB
    peak_mb = peak / (1024 * 1024)
    
    return structure, build_time, peak_mb

def measure_search_time_and_memory(search_func, structure, query, text_length=None, text=None):
    """Measure search time and memory usage using tracemalloc"""
    # Force garbage collection
    gc.collect()
    
    # Start memory tracking
    tracemalloc.start()
    
    # Add error handling
    try:
        # Measure time
        start_time = time.time()
        
        # Call appropriate search function based on parameters
        if text_length is not None and text is None:
            # For suffix tree
            result = search_func(structure, query, text_length)
        elif text is not None and search_func.__module__ == 'suffix_array':
            # For suffix array - note the order of parameters
            result = search_func(text, structure, query)
        else:
            # Generic case
            result = search_func(structure, query)
            
        search_time = time.time() - start_time
        
        # Get memory stats
        current, peak = tracemalloc.get_traced_memory()
        
        # Convert to MB
        peak_mb = peak / (1024 * 1024)
        
        # Stop tracking
        tracemalloc.stop()
        
        return result, search_time, peak_mb
    except Exception as e:
        # Stop tracking if there was an error
        tracemalloc.stop()
        print(f"Error during search: {e}")
        return (0, []), 0, 0  # Return default values on error

def benchmark_structures(reference_file, queries):
    # Read reference sequence
    reference = utils.read_fasta(reference_file)
    text = reference[0][1]  # Assuming read_fasta returns a list of (header, sequence) tuples
    
    # Debug the reference structure
    print(f"Reference type: {type(reference)}")
    print(f"Reference length: {len(reference)}")
    
    # Ensure text is a string
    if not isinstance(text, str):
        text = ''.join(map(str, text))
    
    print(f"Successfully read reference sequence of length {len(text)}")
    print(f"Text type: {type(text)}")
    print(f"First 20 characters: {text[:20]}")
    
    # Group queries by length
    query_groups = categorize_queries(queries)
    
    # Results will be stored here
    results = {
        'tree': {'build': {}, 'search': {}},
        'array': {'build': {}, 'search': {}}
    }
    
    # Build structures and measure performance
    print("Building suffix structures...")
    
    # Suffix Tree
    print("Building suffix tree...")
    tree, tree_build_time, tree_build_mem = measure_build_time_and_memory(
        suffix_tree.build_suffix_tree, text
    )
    results['tree']['build'] = {'time': tree_build_time, 'memory': tree_build_mem}
    
    # Suffix Array
    print("Building suffix array...")
    array, array_build_time, array_build_mem = measure_build_time_and_memory(
        suffix_array.build_suffix_array, text
    )
    results['array']['build'] = {'time': array_build_time, 'memory': array_build_mem}
    
    print("Measuring search performance...")
    
    # For each query length
    for length, queries_of_length in query_groups.items():
        print(f"Processing queries of length {length}...")
        results['tree']['search'][length] = {'time': [], 'memory': []}
        results['array']['search'][length] = {'time': [], 'memory': []}
        
        # For each query in this length group
        for i, query in enumerate(queries_of_length):
            if i % 10 == 0:
                print(f"  Processing query {i+1}/{len(queries_of_length)}")
                
            # Suffix Tree
            tree_result, tree_search_time, tree_search_mem = measure_search_time_and_memory(
                suffix_tree.search_tree, tree, query, len(text)
            )
            results['tree']['search'][length]['time'].append(tree_search_time)
            results['tree']['search'][length]['memory'].append(tree_search_mem)
            
            # Suffix Array - handle the different return format
            array_result, array_search_time, array_search_mem = measure_search_time_and_memory(
                suffix_array.search_array, array, query, text=text
            )
            
            # Handle the case where array_result might have 3 values
            if isinstance(array_result, tuple) and len(array_result) >= 2:
                # Only take the first two values (count and matches)
                array_result = (array_result[0], array_result[1])
                
            results['array']['search'][length]['time'].append(array_search_time)
            results['array']['search'][length]['memory'].append(array_search_mem)
    
    return results

def plot_results(results, output_dir):
    """Generate plots from benchmark results"""
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Plot build times
    plt.figure(figsize=(10, 6))
    structures = ['tree', 'array']
    build_times = [results[s]['build']['time'] for s in structures]
    build_mem = [results[s]['build']['memory'] for s in structures]
    
    plt.subplot(1, 2, 1)
    bars = plt.bar(structures, build_times)
    plt.title('Build Time')
    plt.ylabel('Time (seconds)')
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                 f'{height:.2f}s',
                 ha='center', va='bottom')
    
    plt.subplot(1, 2, 2)
    bars = plt.bar(structures, build_mem)
    plt.title('Build Memory Usage')
    plt.ylabel('Memory (MB)')
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                 f'{height:.2f}MB',
                 ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'build_performance.png'))
    plt.close()
    
    # Plot search times by query length
    plt.figure(figsize=(12, 8))
    
    # Get all unique query lengths
    query_lengths = sorted(results['tree']['search'].keys())
    
    # Average search times for each length
    avg_search_times = {
        'tree': [np.mean(results['tree']['search'][l]['time']) for l in query_lengths],
        'array': [np.mean(results['array']['search'][l]['time']) for l in query_lengths]
    }
    
    # Average memory usage for each length
    avg_search_mem = {
        'tree': [np.mean(results['tree']['search'][l]['memory']) for l in query_lengths],
        'array': [np.mean(results['array']['search'][l]['memory']) for l in query_lengths]
    }
    
    plt.subplot(2, 1, 1)
    plt.plot(query_lengths, avg_search_times['tree'], 's-', label='Suffix Tree')
    plt.plot(query_lengths, avg_search_times['array'], '^-', label='Suffix Array')
    plt.title('Average Search Time by Query Length')
    plt.xlabel('Query Length')
    plt.ylabel('Time (seconds)')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    plt.plot(query_lengths, avg_search_mem['tree'], 's-', label='Suffix Tree')
    plt.plot(query_lengths, avg_search_mem['array'], '^-', label='Suffix Array')
    plt.title('Average Search Memory Usage by Query Length')
    plt.xlabel('Query Length')
    plt.ylabel('Memory (MB)')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'search_performance.png'))
    plt.close()
    
    # Create search performance tables with detailed statistics
    for metric in ['time', 'memory']:
        plt.figure(figsize=(12, 6))
        
        # Table data
        labels = []
        tree_values = []
        array_values = []
        
        for length in query_lengths:
            labels.append(f"Length {length}")
            tree_data = results['tree']['search'][length][metric]
            array_data = results['array']['search'][length][metric]
            
            tree_values.append([
                np.mean(tree_data),
                np.min(tree_data),
                np.max(tree_data),
                np.std(tree_data)
            ])
            
            array_values.append([
                np.mean(array_data),
                np.min(array_data),
                np.max(array_data),
                np.std(array_data)
            ])
        
        # Create figure for the table
        plt.figure(figsize=(12, len(query_lengths) + 2))
        table_data = [
            ['Query Length', 'Structure', 'Mean', 'Min', 'Max', 'Std Dev']
        ]
        
        for i, length in enumerate(query_lengths):
            table_data.append([
                length, 'Tree', 
                f"{tree_values[i][0]:.4f}", 
                f"{tree_values[i][1]:.4f}", 
                f"{tree_values[i][2]:.4f}", 
                f"{tree_values[i][3]:.4f}"
            ])
            table_data.append([
                '', 'Array', 
                f"{array_values[i][0]:.4f}", 
                f"{array_values[i][1]:.4f}", 
                f"{array_values[i][2]:.4f}", 
                f"{array_values[i][3]:.4f}"
            ])
        
        plt.axis('off')
        plt.table(
            cellText=table_data[1:],
            colLabels=table_data[0],
            loc='center',
            cellLoc='center',
            colWidths=[0.1, 0.1, 0.2, 0.2, 0.2, 0.2]
        )
        
        plt.title(f'Search {metric.capitalize()} Statistics')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'search_{metric}_stats.png'))
        plt.close()
    
    # Save raw data
    with open(os.path.join(output_dir, 'benchmark_results.txt'), 'w') as f:
        f.write('Build Performance:\n')
        for structure in structures:
            f.write(f"{structure} - Time: {results[structure]['build']['time']:.4f}s, Memory: {results[structure]['build']['memory']:.4f}MB\n")
        
        f.write('\nSearch Performance:\n')
        for length in query_lengths:
            f.write(f"\nQuery Length: {length}\n")
            for structure in structures:
                times = results[structure]['search'][length]['time']
                mems = results[structure]['search'][length]['memory']
                
                f.write(f"{structure} - Time (s): Mean={np.mean(times):.4f}, Min={np.min(times):.4f}, "
                        f"Max={np.max(times):.4f}, StdDev={np.std(times):.4f}\n")
                        
                f.write(f"{structure} - Memory (MB): Mean={np.mean(mems):.4f}, Min={np.min(mems):.4f}, "
                        f"Max={np.max(mems):.4f}, StdDev={np.std(mems):.4f}\n")

def main():
    args = get_args()
    
    # Read queries
    queries = read_queries(args.queries)
    
    # Run benchmarks
    results = benchmark_structures(args.reference, queries)
    
    # Plot results
    plot_results(results, args.output)
    
    print(f"Benchmark completed. Results saved to {args.output}")

if __name__ == '__main__':
    main()