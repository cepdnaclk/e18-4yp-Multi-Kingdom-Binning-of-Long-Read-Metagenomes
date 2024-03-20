def process_hmmout_file(file_path, threshold):
    """
    Process the .hmmout file and get unique reads, their most matching marker gene, and score
    :param file_path: Path to the .hmmout file
    :param threshold: Threshold value for score
    :return: Dictionary containing unique reads, their most matching marker gene, and score
    """
    # Dictionary to store data for each target read
    score_data = {}

    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                
                data = line.split()
                target_read = data[0].split("_")
                
                # Extracting target ID, calculating score, and retrieving marker gene
                target_id = target_read[0]
                temp_score = (int(data[18]) - int(data[17]) + 1 )/int(data[5])*100
                temp_marker_gene = data[3]
                
                # Update data if target_id exists, otherwise add new entry if score is above threshold
                if target_id in score_data:
                    if score_data[target_id]['score'] < temp_score:
                        score_data[target_id]['score'] = temp_score
                        score_data[target_id]['marker_gene'] = temp_marker_gene
                else:
                    if temp_score >= threshold:
                        score_data[target_id] = {'marker_gene': temp_marker_gene, 'score': temp_score}
                        
    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    # Write the processed data to a text file
    output_file_path = 'processed_marker_scores.txt'
    with open(output_file_path, 'w') as output_file:
        for key,value in score_data.items():
            output_file.write(f"{key}\t{value['marker_gene'].ljust(20)}\t{value['score']}\n")

    print("Processed data has been written to 'processed_marker_scores.txt'")
        
    return score_data

