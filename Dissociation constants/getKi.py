from zeep import Client
from zeep.helpers import serialize_object
import hashlib
import time
import csv
import os

# Define WSDL URL
wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

# Placeholder email and password 
username = "j.doe@example.edu"
password = hashlib.sha256("password".encode("utf-8")).hexdigest()

# Initialize the SOAP client
client = Client(wsdl)

# Retry function
def fetch_with_retry(func, *args, **kwargs):
    for attempt in range(3):  # Retry up to 3 times
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(2)  # Wait before retrying
    raise Exception("All retry attempts failed.")

# Fetch EC numbers
def get_ec_numbers():
    ec_response = fetch_with_retry(client.service.getEcNumbersFromKiValue, 
                                   username, 
                                   password)
    return list(ec_response)  # Convert the response to a list

# Parse Ki response into structured data
def parse_ki_response(ki_response):
    data = []
    if ki_response:
        for entry in ki_response:
            print(f"Parsing Ki entry: {entry}")
            # Access fields as dictionary keys
            data.append({
                "ecNumber": entry.get("ecNumber"),
                "kiValue": entry.get("kiValue"),
                "kiValueMaximum": entry.get("kiValueMaximum"),
                "inhibitor": entry.get("inhibitor"),
                "commentary": entry.get("commentary"),
                "organism": entry.get("organism"),
                "ligandStructureId": entry.get("ligandStructureId"),
                "literature": entry.get("literature")
            })
    return data

# Save data to a CSV file
def save_to_csv(data, filename=None):
    if filename is None:
        # Get path to the current script directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(script_dir, "data")
        os.makedirs(output_dir, exist_ok=True)  # Create the data directory if it doesn't exist
        filename = os.path.join(output_dir, "all_ki_values.csv")
    if data:
        # Define the specific headings for the CSV file
        keys = ["ecNumber", "kiValue", "kiValueMaximum", "inhibitor", "commentary", "organism", "ligandStructureId", "literature"]
        with open(filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()  # Write headers
            writer.writerows(data)  # Write rows
        print(f"Data saved successfully to {filename}")
    else:
        print("No Ki values found to save.")

# Fetch Ki values for a given EC number
def fetch_ki_values(ec_number):
    print(f"Fetching Ki values for EC Number: {ec_number}")
    parameters = (
        username,
        password,
        f"ecNumber*{ec_number}",
        "kiValue*",
        "kiValueMaximum*",
        "inhibitor*",
        "commentary*",
        "organism*Homo sapiens",
        "ligandStructureId*",
        "literature*"
    )
    ki_response = fetch_with_retry(client.service.getKiValue, *parameters)
    if ki_response:
        ki_response = serialize_object(ki_response)  # Serialize the response into a Python object
        print(f"Received Ki response for EC Number: {ec_number}")
        return parse_ki_response(ki_response)
    else:
        print(f"No data found for EC Number: {ec_number}")
        return []

# Main script
def main():
    print("Fetching EC numbers...")
    ec_numbers = get_ec_numbers()
    print(f"Fetched {len(ec_numbers)} EC numbers.")

    ec_numbers = ec_numbers

    all_ki_values = []
    for ec_number in ec_numbers:
        ki_values = fetch_ki_values(ec_number)
        all_ki_values.extend(ki_values)  # Add to the main list

    print("Saving data to ki_values.csv...")
    save_to_csv(all_ki_values)

if __name__ == "__main__":
    main()
