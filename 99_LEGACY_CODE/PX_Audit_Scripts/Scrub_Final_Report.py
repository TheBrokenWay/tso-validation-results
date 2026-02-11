import json
import os
from collections import Counter

def generate_yield_report():
    path = "E:/foundation/PX_Warehouse/WorldLines"
    routes = []
    bands = []
    
    files = [f for f in os.listdir(path) if f.endswith(".worldline")]
    
    for file in files:
        with open(os.path.join(path, file), "r") as f:
            data = json.load(f)
            header = data.get("header", {})
            routes.append(header.get("route", "UNKNOWN"))
            # Pulling the scrub band (HIGH/MEDIUM/LOW)
            bands.append(data.get("physics_snapshot", {}).get("scrub_band", "LOW"))

    route_counts = Counter(routes)
    band_counts = Counter(bands)
    
    print("\n=== PREDATOR X: SCRUB YIELD REPORT ===")
    print(f"Total Processed: {len(files)}")
    print("-" * 30)
    print(f"BLUE (CHAF - PROMOTED):   {route_counts['BLUE']}")
    print(f"RED  (WEED - DECOY):      {route_counts['RED']}")
    print(f"TERM (WEED - BLOCKED):    {route_counts['TERMINATED']}")
    print("-" * 30)
    print(f"HIGH QUALITY (>0.90):     {band_counts['HIGH']}")
    print(f"MED QUALITY  (>0.85):     {band_counts['MEDIUM']}")
    print(f"LOW QUALITY  (<0.85):     {band_counts['LOW']}")

if __name__ == "__main__":
    generate_yield_report()
