def longest_common_substring(seq1, seq2):
    """
    Computes the Longest Common Substring between two DNA sequences.
    
    This algorithm uses dynamic programming to find the longest sequence
    of nucleotides that appear in the same order in both sequences consecutively.
    
    Args:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        
    Returns:
        int: Length of the Longest Common Substring
    """
    # Convert sequences to uppercase for case-insensitive comparison
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    m, n = len(seq1), len(seq2)
    # Create DP table with one extra row and column
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    max_length = 0
    
    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                # If characters match, extend the previous substring
                dp[i][j] = dp[i-1][j-1] + 1
                max_length = max(max_length, dp[i][j])
            else:
                # If characters don't match, reset the substring length
                dp[i][j] = 0
    
    return max_length


def longest_common_subsequence(seq1, seq2):
    """
    Computes the Longest Common Subsequence between two DNA sequences.
    
    This algorithm uses dynamic programming to find the longest sequence
    of nucleotides that appear in the same order in both sequences,
    not necessarily consecutively.
    
    Args:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        
    Returns:
        int: Length of the Longest Common Subsequence
    """
    # Convert sequences to uppercase for case-insensitive comparison
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    m, n = len(seq1), len(seq2)
    
    # Create a 2D DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Fill the DP table
    for i in range(m):
        for j in range(n):
            if seq1[i] == seq2[j]:
                # If characters match, extend the previous LCS
                dp[i + 1][j + 1] = dp[i][j] + 1
            else:
                # Take the maximum from either skipping a character in seq1 or seq2
                dp[i + 1][j + 1] = max(dp[i + 1][j], dp[i][j + 1])
    
    # The bottom-right cell contains the length of the LCS
    return dp[m][n]


def needleman_wunsch(first_str, second_str):
    """
    Computes the Needleman-Wunsch algorithm between two DNA sequences.
    
    This algorithm uses dynamic programming to find the best possible
    DNA sequence alignment. It does this by inserting gaps into the DNA
    sequences, and calculates the best score by minimizing the number of gaps,
    while maximizing the number of DNA sequence pair matchings. 
    
    Args:
        first_str (str): First DNA sequence
        second_str (str): Second DNA sequence
        
    Returns:
        int: Alignment score
    """
    # initialize rows and cols for dp array
    rows = len(first_str) + 1
    cols = len(second_str) + 1

    # initialize dp table to save all values
    dp_arr = [[0] * cols for i in range(rows)]

    ## fill in gap points
    for i in range(1, cols):
        dp_arr[0][i] = -1 * i
    for i in range(1, rows):
        dp_arr[i][0] = -1 * i

    # iterate through every row, col and compute the values from the top, left and diagonal
    # top represents inserting a gap in sequence 2 and left represents inserting a gap in sequence 1
    # diagonal represents a character alignment and we compare if it is a match or mismatch
    for row in range(1, rows):
        for col in range(1, cols):
            top = dp_arr[row - 1][col] - 1
            left = dp_arr[row][col - 1] - 1
            diagonal = dp_arr[row - 1][col - 1] + 1 if first_str[row - 1] == second_str[col - 1] else dp_arr[row - 1][col - 1] - 1
            dp_arr[row][col] = max(top, left, diagonal)

    return dp_arr[rows - 1][cols - 1]


def editDistance(query, sequence):
    '''
    Computes the Edit Distance Algorithm to compare two sequences.

    This algorithm uses dynamic programming to find the smallest number of
    operations possible to transform one sequence into another, using insertion,
    dleetion, or substitution.

        Args:
        query (str): First DNA sequence
        sequence (str): Second DNA sequence
        
    Returns:
        int: number of operations performed

    '''

    # Convert sequences to uppercase for case-insensitive comparison
    query = query.upper()
    sequence = sequence.upper()

    # counter for the operations performed
    fixes = []

    m, n = len(query), len(sequence)
    
    # Create a 2D DP table
    dp = [[0] * (n + 1) for i in range(m + 1)]
    
    # Fill the DP table
    for i in range(m + 1):
        for j in range(n + 1):
            current = []
            # fill in base cases
            if i == 0:
               dp[i][j] = j
            elif j == 0:
               dp[i][j] = i

            else:
                current.append(dp[i - 1][j] + 1)
                current.append(dp[i][j - 1] + 1)
                
                # do the characters match?
                if query[i - 1] != sequence[j - 1]:
                    current.append(dp[i-1][j-1] + 2)
                else:
                    current.append(dp[i-1][j-1])

                # choose lowest-cost operation
                dp[i][j] = min(current)

    # return number of operations you did on the sequence
    return dp[m][n]


def read_fasta(file_path):
    """
    Reads a FASTA format file and returns a dictionary mapping sequence headers to sequences.
    
    Args:
        file_path (str): Path to the FASTA file
        
    Returns:
        dict: Dictionary with headers as keys and sequences as values
    """
    sequences = {}
    current_header = None
    current_sequence = ""
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                # If line is a header
                if line.startswith('>'):
                    # Save the previous sequence if it exists
                    if current_header is not None:
                        sequences[current_header] = current_sequence
                    # Start a new sequence
                    current_header = line
                    current_sequence = ""
                # If line is a sequence
                elif line:
                    current_sequence += line
            
            # Add the last sequence
            if current_header is not None:
                sequences[current_header] = current_sequence
                
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found")
        return {}
        
    return sequences


# Main execution
if __name__ == "__main__":
    import os
    import argparse
    import sys
    
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='DNA Sequence Comparison Tool')
    parser.add_argument('--query', type=str, help='Path to query sequence file')
    parser.add_argument('--database', type=str, help='Path to database sequences file')
    parser.add_argument('--algorithm', type=str, choices=['lcs', 'substring', 'needleman', 'editdistance', 'all'],
                        help='Algorithm to use for comparison')
    
    args = parser.parse_args()
    
    # Look for files in current directory by default
    # Use the specific file paths
    query_file = "DNA_query.txt"
    database_file = "DNA_sequences.txt"
    
    print(f"Attempting to read files")
    print(f"Query file: {query_file}")
    print(f"Database file: {database_file}")
    
    # Read the query sequence
    query_sequences = read_fasta(query_file)
    if not query_sequences:
        print("No query sequence found. Checking for raw query text...")
        # Try to read the query file directly as raw text (not FASTA format)
        try:
            with open(query_file, 'r') as f:
                raw_query = f.read().strip()
                query_sequences = {">Query": raw_query}
                print("Read query as raw text")
        except:
            print("Could not read query file as raw text either.")
            exit(1)
    
    # Just take the first sequence from the query file
    query_header, query_sequence = next(iter(query_sequences.items()))
    print(f"Query sequence: {query_header}")
    print(f"Query length: {len(query_sequence)} nucleotides")
    
    # Read the database of sequences
    database_sequences = read_fasta(database_file)
    if not database_sequences:
        print("No database sequences found.")
        exit(1)
    
    print(f"Found {len(database_sequences)} sequences in the database.")
    
    # Choose algorithm
    selected_algorithm = args.algorithm
    
    # If no algorithm was specified, ask user to choose
    if not selected_algorithm:
        print("\nAvailable algorithms:")
        print("1. Longest Common Subsequence (lcs)")
        print("2. Longest Common Substring (substring)")
        print("3. Needleman-Wunsch (needleman)")
        print("4. Edit Distance (editdistance)")
        print("5. Run all algorithms (all)")
        
        try:
            choice = input("\nEnter your choice (1-5): ")
            if choice == '1':
                selected_algorithm = 'lcs'
            elif choice == '2':
                selected_algorithm = 'substring'
            elif choice == '3':
                selected_algorithm = 'needleman'
            elif choice == '4':
                selected_algorithm = 'editdistance'
            else:
                selected_algorithm = 'all'
        except:
            # If input fails, default to all
            print("Input error. Running all algorithms.")
            selected_algorithm = 'all'
    
    # Run the selected algorithm(s)
    results = {}
    
    # Longest Common Subsequence
    if selected_algorithm in ['lcs', 'all']:
        best_header_lcs = None
        best_score_lcs = -1
        
        print("\nComputing LCSubsequence scores...")
        for header, sequence in database_sequences.items():
            score = longest_common_subsequence(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = ""
            if len(header.split()) > 0:
                for i in range(1, len(header.split())):
                    short_header += header.split()[i] + " "
            
            print(f"LCSubsequence Score for {short_header}: {score}")
            
            if score > best_score_lcs:
                best_score_lcs = score
                best_header_lcs = header
        
        results['lcs'] = (best_header_lcs, best_score_lcs)
    
    # Longest Common Substring
    if selected_algorithm in ['substring', 'all']:
        best_header_substring = None
        best_score_substring = -1
        
        print("\nComputing LCSubstring scores...")
        for header, sequence in database_sequences.items():
            score = longest_common_substring(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = ""
            if len(header.split()) > 0:
                for i in range(1, len(header.split())):
                    short_header += header.split()[i] + " "
            
            print(f"LCSubstring Score for {short_header}: {score}")
            
            if score > best_score_substring:
                best_score_substring = score
                best_header_substring = header
        
        results['substring'] = (best_header_substring, best_score_substring)
    
    # Needleman-Wunsch
    if selected_algorithm in ['needleman', 'all']:
        best_header_nw = None
        best_score_nw = -float('inf')
        
        print("\nComputing Needleman-Wunsch scores...")
        for header, sequence in database_sequences.items():
            score = needleman_wunsch(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = ""
            if len(header.split()) > 0:
                for i in range(1, len(header.split())):
                    short_header += header.split()[i] + " "
            
            print(f"N-W Score for {short_header}: {score}")
            
            if score > best_score_nw:
                best_score_nw = score
                best_header_nw = header
        
        results['needleman'] = (best_header_nw, best_score_nw)
    
    # Edit Distance
    if selected_algorithm in ['editdistance', 'all']:
        best_header_ed = None
        best_score_ed = float('inf')  # Lower is better for edit distance
        
        print("\nComputing Edit Distance scores...")
        for header, sequence in database_sequences.items():
            score = editDistance(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = ""
            if len(header.split()) > 0:
                for i in range(1, len(header.split())):
                    short_header += header.split()[i] + " "
            
            print(f"Edit Distance Score for {short_header}: {score}")
            
            # For edit distance, lower is better
            if score < best_score_ed:
                best_score_ed = score
                best_header_ed = header
        
        results['editdistance'] = (best_header_ed, best_score_ed)
    
    # Print summary of results
    print("\n=== RESULTS SUMMARY ===")
    
    for alg, (header, score) in results.items():
        algorithm_name = {
            'lcs': 'Longest Common Subsequence',
            'substring': 'Longest Common Substring',
            'needleman': 'Needleman-Wunsch',
            'editdistance': 'Edit Distance'
        }.get(alg)
        
        print(f"\nMost similar sequence using {algorithm_name}:")
        
        short_header = ""
        if len(header.split()) > 0:
            for i in range(1, len(header.split())):
                short_header += header.split()[i] + " "

        print(f"{short_header}")
        print(f"Score: {score}")