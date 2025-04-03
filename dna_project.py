def editDistance(query, sequence):
    # counter for the operations performed
    fixes = 0
    
    # check if they're the same length
    if len(sequence) != len(query):
        # number of insertions/deletes that need to happen
        l = len(sequence) - len(query)
        fixes += abs(l)
        
        # if the query is bigger add the missing characters to the sequence 
        if l < 0:
            extra = query[abs(l):]
            sequence += extra
        
        # if sequence is bigger than query delete the extra
        else:
            sequence = sequence[0:len(query)]

    # check if a substitution needs to be made
    for idx in range(len(query)):
        if sequence[idx] != query[idx]:
            fixes += 1

    # return number of operations you did on the sequence
    return fixes


def longest_common_substring(seq1,seq2):
  """
    Computes the Longest Common Substring between two DNA sequences.
    
    This algorithm uses a table to find the longest sequence
    of nucleotides that appear in the same order in both sequences consecutively.
    
    Args:
        seq1 (str): First DNA sequence
        seq2 (str): Second DNA sequence
        
    Returns:
        int: Length of the Longest Common Subsequence
    """
  
  # Convert sequences to uppercase for case-insensitive comparison
  seq1 = seq1.upper()
  seq2 = seq2.upper()
    
  # create table to store string matches
  dp_array = [[0] * len(seq1) for i in range(len(seq2))]

  # save rows and cols
  rows, cols = len(dp_array), len(dp_array[0])
  longest_str = ""

  # check every character match in both sequences and mark with a 1 if there is a match
  for i in range(len(seq2)):
    for j in range(len(seq1)):
      if seq1[j] == seq2[i]:
        dp_array[i][j] = 1
  
  # check diagonals from bottom up 
  for i in range(rows - 1, -1, -1):
    current_str = ""
    row, col = i, 0
    
    # if the current diagonal position contains a match, add the string 
    while row < rows and col < cols:
        # longest substring will be contained in the longest consecutive diagonals of 1s
        if dp_array[row][col] == 1:
          current_str += seq2[row]
          
        # check for longest length seen so far
        if len(current_str) > len(longest_str):
          longest_str = current_str

        # move to next cell in diagonal
        row += 1
        col += 1

  # repeat and check for the remaining diagonals
  for i in range(1, cols):
    current_str = ""
    row, col = 0, i

    # check current diagonal and if there is a character match, add it to current string
    while row < rows and col < cols:
        # longest substring will be contained in the longest consecutive diagonals of 1s
        if dp_array[row][col] == 1:
          current_str += seq2[row]

        # check for longest string seen so far
        if len(current_str) > len(longest_str):
          longest_str = current_str

        # move to next cell in diagonal
        row += 1
        col += 1

  return len(longest_str)



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
    Computes the needleman wunsch algorithm between two DNA sequences.
    
    This algorithm uses dynamic programming to find the best possible
    dna sequence allignment. It does this by inserting gaps into the dna
    sequences, and calculates the best score by minimizing the number of gaps,
    while maximizing the number of dna sequence pair matchings. 
    
    Args:
        first_str (str): First DNA sequence
        second_str (str): Second DNA sequence
        
    Returns:
        int: Length of the alignment score
    """
  # initialize rows and cols for dp array
  rows = len(first_str) + 1
  cols = len(second_str) + 1

  # initialize dp table to save all values
  dp_arr = [[0] * cols for i in range(rows)]

  ## fill in gap points
  for i in range(1,cols):
    dp_arr[0][i] = -1 * i
  for i in range(1,rows):
    dp_arr[i][0] = -1 * i

  # iterate through every row, col and compute the values from the top, left and diagonal
  # top represents inserting a gap in sequence 2 and left represents inserting a gap in sequence 1
  # diagonal represents a character alignment and we compare if it is a match or mismatch
  for row in range(1, rows):
    for col in range(1, cols):
      top = dp_arr[row - 1][col] - 1
      left = dp_arr[row][col - 1] - 1
      diagonal = dp_arr[row - 1][col - 1] + 1 if first_str[row - 1] == second_str[col - 1] else dp_arr[row - 1][col - 1] - 1
      dp_arr[row][col] = max(top,left,diagonal)


  return dp_arr[rows - 1][cols - 1]



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
    query_file = input("Enter query file name: ")
    database_file = input("Enter sequences file name: ")
    
    print(f"\nAttempting to read files in current directory")
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
    
    print("\nEnter 1 for LCSubsequence, 2 for LCSubstring, 3 for Needleman-Wunsch, and 4 for Edit Distance")
    algorithm = input("Choose algorithm: ")

    if algorithm == "1":
        # Find the most similar sequence using LCS
        best_score = -1
        best_header = None
        
        print("\nComputing LCSubsequence scores...")
        for header, sequence in database_sequences.items():
            score = longest_common_subsequence(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = header.split()[0] if len(header.split()) > 0 else header
            print(f"LCSubseqence Score for {short_header}: {score}")
            
            if score > best_score:
                best_score = score
                best_header = header
        
        print("\nMost similar sequence:")
        print(best_header)
        print(f"LCSequence Score: {best_score}")
    
    elif algorithm == "2":
        # Find the most similar sequence using longest common substring
        best_score2 = -1
        best_header2 = None
        
        print("\nComputing LCSubstring scores...")
        for header, sequence in database_sequences.items():
            score = longest_common_substring(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = header.split()[0] if len(header.split()) > 0 else header
            print(f"LCSubstring Score for {short_header}: {score}")
            
            if score > best_score2:
                best_score2 = score
                best_header2 = header
        
        print("\nMost similar sequence:")
        print(best_header2)
        print(f"LCSubstring Score: {best_score2}")
    
    elif algorithm == "3":
        # Find the most similar sequence using needleman wunsch algorithm
        best_score3 = -1
        best_header3 = None
        
        print("\nComputing Needleman-Wunsch scores...")
        for header, sequence in database_sequences.items():
            score = needleman_wunsch(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = header.split()[0] if len(header.split()) > 0 else header
            print(f"N-W Score for {short_header}: {score}")
            
            if score > best_score3:
                best_score3 = score
                best_header3 = header
        
        print("\nMost similar sequence:")
        print(best_header3)
        print(f"Needleman Wunsch Score: {best_score3}")
    
    elif algorithm == "4":
        # Find the most similar sequence using edit distance algorithm
        best_score4 = 10000
        best_header4 = None
        
        print("\nComputing Edit Distance scores...")
        for header, sequence in database_sequences.items():
            score = editDistance(query_sequence, sequence)
            
            # Print a short version of the header for easier reading
            short_header = header.split()[0] if len(header.split()) > 0 else header
            print(f"Edit Dist Score for {short_header}: {score}")
            
            if score < best_score4:
                best_score4 = score
                best_header4 = header
    
        print("\nMost similar sequence:")
        print(best_header4)
        print(f"Edit Distance Score: {best_score4}")