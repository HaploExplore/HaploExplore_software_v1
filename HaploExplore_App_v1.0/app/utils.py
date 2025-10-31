import logging
import os
import time

def ask_confirmation(prompt : str) -> bool:
    """
    Ask a yes/no question to the user : answer is mandatory.
    """
    while True:
        confirmation = input(prompt).strip().lower()
        if confirmation in ['yes', 'y']:
            return True
        if confirmation in ['no', 'n']:
            return False
        print("Invalid input. Please enter 'yes' or 'no'.")

def is_integer(value) -> bool :
    """
    Checks that a value can be converted into an int :
    """
    try:
        int(value)
        return True
    except ValueError:
        return False

def minimum_distance(positions : list[float|int]) -> float|int:
    """
    Not used yet : determines the smallest distance between two SNPs from their list of positions.
    Possible usage : optimize plotting.
    """
    positions.sort()
    differences = [positions[i+1] - positions[i] for i in range(len(positions) - 1)]
    min_distance = min(differences)
    return min_distance, positions

def find_position_within_sorted_list(target : int,  sorted_list : list[int]) -> int : 
    """
    Find the index of a specific integer within a sorted list.
    """
    start, end = 0, len(sorted_list) - 1

    while start <= end :
        mid = (start + end) // 2

        if sorted_list[mid] == target :
            return mid
        if sorted_list[mid] < target :
            start = mid + 1
        else :
            end = mid - 1
        
    return -1

def log_time() -> str:
    """Return the current timestamp as a string"""
    return time.strftime('%Y-%m-%d %H:%M:%S')

def generate_log_filename(base_name: str = "haploblocks_analysis") -> str:
    """Generate a unique log file name with a timestamp in an 'output' directory."""
    timestamp = log_time().replace(" ", "_").replace(":", "-")
    script_dir = os.path.dirname(os.path.abspath(__file__))  
    output_dir = os.path.join(script_dir, "output")  

    os.makedirs(output_dir, exist_ok=True)  
    return os.path.join(output_dir, f"{base_name}_{timestamp}.log")

def setup_logging(log_filepath: str):
    """Set up logging to a log file and console, ensuring logs are written immediately."""
    print(f"Log file path: {log_filepath}")
    os.makedirs(os.path.dirname(log_filepath), exist_ok=True)  

    # Remove all existing handlers to avoid duplicates
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filepath, mode='w', encoding='utf-8'),  # Log to file
            logging.StreamHandler()  # Log to console
        ]
    )

    # Ensure logs are flushed
    logging.getLogger().handlers[0].flush()

    print(f"Logging initialized. Output will be saved to {log_filepath}")
