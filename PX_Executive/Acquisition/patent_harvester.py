"""
PX_Executive/Acquisition/patent_harvester.py

The Sovereign Patent Acquisition Engine.
Strategy:
1. Ingest Local Bulk Data (Preferred, 100% Reliable).
2. Fetch Remote Data (Fallback, High Latency, Retries).

Output:
- Normalizes data into PX_Warehouse/Operations/Inputs/patents/
- Decouples the runtime Pipeline from the fickle Internet.
"""
import os
import json
import time
import random
import requests
from datetime import datetime
from pathlib import Path

# --- CONFIGURATION ---
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RAW_DATA_DIR = REPO_ROOT / "PX_Data" / "patents" / "raw"
WAREHOUSE_INBOX = REPO_ROOT / "PX_Warehouse" / "Operations" / "Inputs" / "patents"

# Ensure directories exist
RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)
WAREHOUSE_INBOX.mkdir(parents=True, exist_ok=True)

class PatentHarvester:
    def __init__(self):
        self.session = requests.Session()
        # Pretend to be a browser to avoid some basic bot-blocking
        self.session.headers.update({
            "User-Agent": "PredatorX/Enterprise (Research; contact@predatorx.ai)"
        })

    def ingest_bulk_file(self, file_path):
        """
        Mode 1: BULK INGESTION (100% Reliable)
        Reads a local JSONL file (like seed.jsonl) and stages it for the pipeline.
        """
        print(f"--- [HARVEST] Ingesting Bulk File: {file_path} ---")
        path = Path(file_path)
        if not path.exists():
            print(f"CRITICAL: File not found: {path}")
            return 0

        count = 0
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                try:
                    record = json.loads(line)
                    self._materialize_to_warehouse(record)
                    count += 1
                except json.JSONDecodeError:
                    print(f"WARN: Skipping malformed line in {file_path}")
        
        print(f"--- [HARVEST] Success. Staged {count} patents to Warehouse. ---")
        return count

    def fetch_from_api_robust(self, patent_numbers):
        """
        Mode 2: ROBUST API FETCH (For fresh data updates)
        Uses exponential backoff to handle the flaky API.
        """
        print(f"--- [HARVEST] Attempting API Fetch for {len(patent_numbers)} patents ---")
        success_count = 0
        
        for p_id in patent_numbers:
            # Retry Loop
            max_retries = 5
            for attempt in range(max_retries):
                try:
                    # REPLACE THIS URL with your specific API endpoint (e.g. PatentsView)
                    # This is a generic robust fetch pattern.
                    url = f"https://api.patentsview.org/patents/query?q={{\"patent_number\":\"{p_id}\"}}"
                    
                    response = self.session.get(url, timeout=10)
                    
                    if response.status_code == 200:
                        data = response.json()
                        # Extract the actual record from the API response structure
                        if "patents" in data and data["patents"]:
                            record = data["patents"][0]
                            # Normalize fields to match our schema
                            normalized = {
                                "patent_number": record.get("patent_number"),
                                "title": record.get("patent_title"),
                                "abstract": record.get("patent_abstract"),
                                "grant_date": record.get("patent_date"),
                                "assignees": record.get("assignees", []),
                                "ingested_at": datetime.utcnow().isoformat()
                            }
                            self._materialize_to_warehouse(normalized)
                            success_count += 1
                            print(f"  [OK] Fetched {p_id}")
                            break # Success, exit retry loop
                        else:
                             print(f"  [404] Patent {p_id} not found.")
                             break
                    
                    elif response.status_code == 429:
                        # Rate Limited - Sleep and Retry
                        wait = (2 ** attempt) + random.uniform(0, 1)
                        print(f"  [429] Rate limited. Sleeping {wait:.2f}s...")
                        time.sleep(wait)
                    else:
                        print(f"  [ERR] HTTP {response.status_code} for {p_id}")
                        break

                except Exception as e:
                    print(f"  [EXC] Attempt {attempt+1}/{max_retries} failed: {e}")
                    time.sleep(2) # Basic cool-down
            
            # Global rate limit pause between items
            time.sleep(0.5)

        print(f"--- [HARVEST] API Run Complete. Fetched {success_count}/{len(patent_numbers)} ---")

    def _materialize_to_warehouse(self, record):
        """Writes the record to the Warehouse Inbox. The Pipeline reads from HERE."""
        p_id = record.get("patent_number", "UNKNOWN")
        # Sanitize filename
        safe_id = "".join([c for c in p_id if c.isalnum() or c in ('-','_')])
        
        outfile = WAREHOUSE_INBOX / f"PATENT_{safe_id}.json"
        with open(outfile, "w", encoding="utf-8") as f:
            json.dump(record, f, indent=2)

if __name__ == "__main__":
    harvester = PatentHarvester()
    
    # 1. Try to ingest the local seed file first (The User's specific request)
    seed_file = RAW_DATA_DIR / "seed.jsonl"
    if seed_file.exists():
        harvester.ingest_bulk_file(seed_file)
    else:
        print(f"No seed file found at {seed_file}. Run your PowerShell snippet first!")