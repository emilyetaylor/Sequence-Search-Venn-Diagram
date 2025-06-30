import requests
import pandas as pd
import time

# submit_id_mapping(from_db, to_db, df) - this function will submit an list of accessions, to uniprot, to get each accession's
# unique id for another database. It will return the number of ids submitted for mapping, and will return a job id string
# INPUT from_db - the database where the provided ids are coming from
# INPUT to_db - the database where the provided ids need to be mapped to (given new ids)
# INPUT df - the dataframe containing the original accessions to be mapped
# RETURN response.json()["jobID"] - job id of the mapping done by uniprot
def submit_id_mapping(from_db, to_db, df):
    accession_list = df["accession"].dropna().unique().tolist()
    url = "https://rest.uniprot.org/idmapping/run"
    headers = {
        "Content-Type" : "application/x-www-form-urlencoded"
    }
    payload = {
        "from" : from_db,
        "to" : to_db,
        "ids" : ",".join(accession_list)
    }

    print("=== SUBMITTING ID MAPPING REQUEST ===")
    print("From:", from_db)
    print("To:", to_db)
    print("First 5 IDs:", accession_list[:5])
    print("Total IDs:", len(accession_list))
    response = requests.post(url, data = payload, headers=headers)
    response.raise_for_status()
    return response.json()["jobId"]


def get_entry(job_id):
    while True:
        status = check_status(job_id)
        if status["jobStatus"] == "FINISHED":
            break
        print("Waiting for mapping to finish")
        time.sleep(3)


def check_status(job_id):
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    response = requests.get(status_url)
    response.raise_for_status()
    return response.json()


def get_results(job_id):
    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    response = requests.get(result_url)
    response.raise_for_status()
    return response.json()

def parse_mapped_results(json_data):
    results = json_data["results"]
    parsed = []
    for entry in results:
        from_id = entry["from"]
        to_id = entry["to"]
        parsed.append({"uniprot": from_id, "mapped_id": to_id})
    return pd.DataFrame(parsed)
