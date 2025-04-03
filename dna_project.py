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
    # Look for files directly in the current directory
    query_file = "/Users/autumnle/Downloads/DNA_query.txt"
    database_file = "/Users/autumnle/Downloads/DNA_sequences.txt"
    
    print(f"Attempting to read files in current directory")
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
    
    # Find the most similar sequence using LCS
    best_score = -1
    best_header = None
    
    print("\nComputing LCS scores...")
    for header, sequence in database_sequences.items():
        score = longest_common_subsequence(query_sequence, sequence)
        
        # Print a short version of the header for easier reading
        short_header = header.split()[0] if len(header.split()) > 0 else header
        print(f"LCS Score for {short_header}: {score}")
        
        if score > best_score:
            best_score = score
            best_header = header
    
    print("\nMost similar sequence:")
    print(best_header)
    print(f"LCS Score: {best_score}")
    