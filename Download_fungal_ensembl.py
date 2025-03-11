import time
import csv
from Bio import Entrez

# Set your email address (required by NCBI)
Entrez.email = "cathalmeehan3@gmail.com"

# List of modern AM fungal genera
genus_list = [
    "Acaulospora", "Ambispora", "Archaeospora", "Cetraspora",
    "Claroideoglomus", "Corymbiglomus", "Dentiscutata", "Diversispora",
    "Dominikia", "Funneliformis", "Gigaspora", "Glomus",
    "Kamienskia", "Lucioglomus", "Pacispora", "Paraglomus",
    "Racocetra", "Redeckera", "Rhizophagus", "Scutellospora",
    "Sclerocystis", "Septoglomus", "Simiglomus"
]

# Containers for detailed records
assembly_rows = []
transcriptome_rows = []

for genus in genus_list:
    print(f"Processing genus: {genus}")
    # -------------------------------
    # 1. Search NCBI Assembly database for assemblies
    # -------------------------------
    assembly_query = f'"{genus}"[Organism]'
    try:
        handle = Entrez.esearch(db="assembly", term=assembly_query, retmax=100)
        assembly_record = Entrez.read(handle)
        handle.close()
        assembly_ids = assembly_record.get("IdList", [])
        print(f"  Found {len(assembly_ids)} assembly record(s).")
    except Exception as e:
        print(f"  Error searching assembly for {genus}: {e}")
        assembly_ids = []

    # Retrieve metadata for each assembly
    for aid in assembly_ids:
        try:
            handle = Entrez.esummary(db="assembly", id=aid)
            summary = Entrez.read(handle)
            handle.close()
            # The summary structure is nested; we extract the first document summary.
            detail = summary["DocumentSummarySet"]["DocumentSummary"][0]
            assembly_rows.append({
                "Genus": genus,
                "AssemblyID": aid,
                "AssemblyAccession": detail.get("AssemblyAccession", ""),
                "Organism": detail.get("Organism", ""),
                "SpeciesName": detail.get("SpeciesName", ""),
                "AssemblyName": detail.get("AssemblyName", ""),
                "AssemblyStatus": detail.get("AssemblyStatus", ""),
                "SubmitDate": detail.get("SubmitDate", ""),
                "Version": detail.get("Version", "")
            })
        except Exception as e:
            print(f"  Error retrieving summary for assembly {aid} for {genus}: {e}")
        time.sleep(0.3)  # small delay between requests
    time.sleep(1)

    # -------------------------------
    # 2. Search NCBI Nucleotide database for transcriptome shotgun assemblies (TSA)
    # -------------------------------
    tsa_query = f'"{genus}"[Organism] AND "transcriptome shotgun assembly"[Title]'
    try:
        handle = Entrez.esearch(db="nuccore", term=tsa_query, retmax=100)
        transcriptome_record = Entrez.read(handle)
        handle.close()
        transcriptome_ids = transcriptome_record.get("IdList", [])
        print(f"  Found {len(transcriptome_ids)} transcriptome record(s).")
    except Exception as e:
        print(f"  Error searching transcriptome for {genus}: {e}")
        transcriptome_ids = []

    # Retrieve metadata for each transcriptome record
    for tid in transcriptome_ids:
        try:
            handle = Entrez.esummary(db="nuccore", id=tid)
            summary = Entrez.read(handle)
            handle.close()
            # esummary for nuccore typically returns a list of document summaries
            detail = summary[0]
            transcriptome_rows.append({
                "Genus": genus,
                "TranscriptomeID": tid,
                "Accession": detail.get("AccessionVersion", ""),
                "Title": detail.get("Title", ""),
                "Organism": detail.get("Organism", ""),
                "ReleaseDate": detail.get("ReleaseDate", "")
            })
        except Exception as e:
            print(f"  Error retrieving summary for transcriptome {tid} for {genus}: {e}")
        time.sleep(0.3)
    time.sleep(1)

# -------------------------------
# Write the detailed assembly metadata to a CSV file
# -------------------------------
assembly_output_file = "ncbi_amf_assembly_details.csv"
with open(assembly_output_file, mode="w", newline="") as csvfile:
    fieldnames = ["Genus", "AssemblyID", "AssemblyAccession", "Organism",
                  "SpeciesName", "AssemblyName", "AssemblyStatus", "SubmitDate", "Version"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in assembly_rows:
        writer.writerow(row)

# -------------------------------
# Write the detailed transcriptome metadata to a CSV file
# -------------------------------
transcriptome_output_file = "ncbi_amf_transcriptome_details.csv"
with open(transcriptome_output_file, mode="w", newline="") as csvfile:
    fieldnames = ["Genus", "TranscriptomeID", "Accession", "Title", "Organism", "ReleaseDate"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in transcriptome_rows:
        writer.writerow(row)

print(f"\nAssembly details written to {assembly_output_file}")
print(f"Transcriptome details written to {transcriptome_output_file}")
