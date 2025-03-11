#!/usr/bin/env python
import os
import time
import csv
import requests
from Bio import Entrez, SeqIO

# ================================
# Configuration
# ================================
Entrez.email = "cathalmeehan3@gmail.com"  # Your NCBI email
metadata_csv = "ncbi_amf_assembly_details.csv"  # CSV file with assembly metadata
download_dir = "/Volumes/Elements/Assembly"  # Directory to save downloaded files
chunks_dir = os.path.join(download_dir, "Chunks")  # Directory for processed sequence chunks
combined_chunks_file = "fungal_chunks.txt"  # Final combined training file
chunk_size = 512  # Chunk size in bp
sleep_interval = 0.5  # Pause between downloads (seconds)

# Create directories if they don't exist
for d in [download_dir, chunks_dir]:
    if not os.path.exists(d):
        os.makedirs(d)


# ================================
# Functions
# ================================
def download_file_with_suffix(assembly_accession, assembly_name, file_suffix, out_filepath):
    """
    Downloads a file for a given assembly using a specified suffix.
    The function first attempts to retrieve the FTP path via esummary (with validate=False);
    if that fails, it constructs the FTP URL manually. Then it converts the URL from ftp:// to https://.
    Returns True if download succeeds, False otherwise.
    """
    ftp_path = ""
    try:
        # Attempt to retrieve the assembly summary (skip DTD validation)
        handle = Entrez.esummary(db="assembly", id=assembly_accession)
        summary = Entrez.read(handle, validate=False)
        handle.close()
        docs = summary.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        if docs:
            detail = docs[0]
            ftp_path = detail.get("FtpPath_GenBank", "") or detail.get("FtpPath_RefSeq", "")
    except Exception as e:
        print(f"Error retrieving esummary for {assembly_accession}: {e}")

    if ftp_path:
        folder_name = ftp_path.rstrip("/").split("/")[-1]
    else:
        # Fallback: manually construct the FTP path using the accession and assembly name.
        try:
            parts = assembly_accession.split("_")
            prefix = parts[0]  # e.g. "GCA"
            digits = parts[1].split(".")[0]  # e.g. "910592055"
            d1, d2, d3 = digits[0:3], digits[3:6], digits[6:9]
        except Exception as e:
            print(f"Error parsing accession {assembly_accession}: {e}")
            return False
        folder_name = f"{assembly_accession}_{assembly_name}"
        ftp_path = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{d1}/{d2}/{d3}/{folder_name}"

    filename = f"{folder_name}{file_suffix}"
    # Convert the FTP URL to HTTPS so that requests can handle it.
    if ftp_path.lower().startswith("ftp://"):
        https_path = "https://" + ftp_path[6:]
    else:
        https_path = ftp_path
    download_url = https_path + "/" + filename
    print(f"Downloading {assembly_accession} file:\n  {filename}\n  {download_url}")

    try:
        response = requests.get(download_url, stream=True, timeout=60)
        if response.status_code == 200:
            with open(out_filepath, "wb") as out_file:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        out_file.write(chunk)
            print(f"Downloaded {assembly_accession} file {filename} to {out_filepath}")
            return True
        else:
            print(f"Error downloading {assembly_accession}: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"Error processing {assembly_accession}: {e}")
        return False


def process_fasta_to_chunks(fasta_filepath, chunk_size, chunks_output_dir):
    """
    Processes a gzipped FASTA file by splitting each sequence into fixed-length chunks.
    Each chunk is written as a line in an output text file.
    If the chunk file already exists, it is reused.
    Returns the path to the generated chunks file, or None on error.
    """
    base_name = os.path.basename(fasta_filepath).replace("_genomic.fna.gz", "")
    chunks_file = os.path.join(chunks_output_dir, f"{base_name}_chunks.txt")

    if os.path.exists(chunks_file):
        print(f"Chunks file for {base_name} already exists, skipping processing.")
        return chunks_file

    chunk_count = 0
    import gzip
    try:
        with gzip.open(fasta_filepath, "rt") as handle, open(chunks_file, "w") as out_f:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                for i in range(0, len(seq) - chunk_size + 1, chunk_size):
                    chunk = seq[i:i + chunk_size]
                    out_f.write(chunk + "\n")
                    chunk_count += 1
        print(f"Processed {fasta_filepath}: {chunk_count} chunks saved to {chunks_file}")
    except Exception as e:
        print(f"Error processing {fasta_filepath}: {e}")
        return None

    return chunks_file


# ================================
# Pipeline: Download and Process Assemblies and Annotations
# ================================
all_chunks_files = []

with open(metadata_csv, mode="r", newline="") as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        assembly_acc = row.get("AssemblyAccession", "").strip()
        assembly_name = row.get("AssemblyName", "").strip()
        if assembly_acc and assembly_name:
            # Define output paths for genomic FASTA and annotation (GFF) files.
            fasta_out_path = os.path.join(download_dir, f"{assembly_acc}.fna.gz")
            gff_out_path = os.path.join(download_dir, f"{assembly_acc}.gff.gz")

            # Download FASTA if not already present.
            if os.path.exists(fasta_out_path):
                print(f"{assembly_acc} FASTA already downloaded, skipping download.")
            else:
                success = download_file_with_suffix(assembly_acc, assembly_name, "_genomic.fna.gz", fasta_out_path)
                if not success:
                    continue  # Skip processing if download fails
                time.sleep(sleep_interval)

            # Download Annotation (GFF) if not already present.
            if os.path.exists(gff_out_path):
                print(f"{assembly_acc} GFF annotation already downloaded, skipping download.")
            else:
                success = download_file_with_suffix(assembly_acc, assembly_name, "_genomic.gff.gz", gff_out_path)
                if not success:
                    print(f"Annotation download failed for {assembly_acc}.")
                time.sleep(sleep_interval)

            # Process FASTA into chunks if not already done.
            chunks_file = process_fasta_to_chunks(fasta_out_path, chunk_size, chunks_dir)
            if chunks_file:
                all_chunks_files.append(chunks_file)
        else:
            print("Missing AssemblyAccession or AssemblyName in row, skipping.")

# ================================
# Combine All Chunk Files into a Single Training File
# ================================
with open(combined_chunks_file, "w") as outfile:
    for cf in all_chunks_files:
        with open(cf, "r") as infile:
            outfile.write(infile.read())
print(f"\nAll chunks combined into {combined_chunks_file}")
