import os
import sys

def process_data(file_path):
    plate_data = {}
    
    # Read the file and process each line
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 4:
                print(f"Skipping malformed line: {line}")
                continue

            file_name, read_r1, read_r2, synchronized = parts
            read_r1 = int(read_r1)
            read_r2 = int(read_r2)
            synchronized = synchronized == "True"

            if synchronized:
                # Parse plate ID and barcode set
                _, rest = file_name.split('/')
                plate_id, barcode_r1, _ = rest.split('_b', 2)
                barcode_id_str = barcode_r1.split('_')[0]  # Get only the number part
                barcode_id = int(barcode_id_str)

                if plate_id not in plate_data:
                    plate_data[plate_id] = {'set1': 0, 'set2': 0, 'files1': [], 'files2': []}

                # Determine which set the barcode belongs to and accumulate counts
                if 1 <= barcode_id <= 96:
                    plate_data[plate_id]['set1'] += read_r1
                    plate_data[plate_id]['files1'].append(file_name)
                elif 97 <= barcode_id <= 192:
                    plate_data[plate_id]['set2'] += read_r1
                    plate_data[plate_id]['files2'].append(file_name)

    # Decide which set to keep and delete others
    for plate_id, data in plate_data.items():
        if data['set1'] >= data['set2']:
            set_to_keep = '1-96'
            to_delete = data['files2']
        else:
            set_to_keep = '97-192'
            to_delete = data['files1']
        
        for file in to_delete:
            #% print(f"File to delete: {file}")
            try:
                os.remove(file)
                print(f"Deleted: {file}")
            except FileNotFoundError:
                print(f"Could not find file to delete: {file}")

        # Print the decision to retain which barcode set
        print(f"{plate_id}: Retained Barcode Set {set_to_keep}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python select_barcode_set.py <total.read.counts.txt>")
        sys.exit(1)

    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' does not exist.")
        sys.exit(1)

    process_data(input_file)
    
