import json
import os
import csv
import glob

# CONFIGURATION
TARGET_DIR = r"E:\foundation\PX_Warehouse\SMART_Antiviral_Dossiers"
REPORT_FILE = r"E:\foundation\PX_Warehouse\SMART_DIAMOND_LEADERBOARD.csv"

def rank_diamonds():
    print(f">>> INITIATING DIAMOND AUDIT...")
    
    files = glob.glob(os.path.join(TARGET_DIR, "*_DOSSIER.json"))
    print(f"    Found {len(files)} Dossiers.")
    
    leaderboard = []
    
    for filepath in files:
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
                
            # Extract Key Metrics
            cand_id = data['candidate_id']
            potency = data['profile']['base_potency_nm']
            mech = data['profile']['mechanism']
            breadth = data['profile']['breadth_score']
            tox = data['profile']['toxicity_index']
            
            # Safety Flags (Count how many "CAUTION" flags exist)
            safety = data['population_safety']
            risk_flags = sum(1 for v in safety.values() if v != "ACCEPTABLE")
            
            leaderboard.append({
                "Candidate_ID": cand_id,
                "Potency_nM": potency,
                "Mechanism": mech,
                "Breadth": breadth,
                "Toxicity_Index": tox,
                "Safety_Flags": risk_flags, # Lower is better
                "Pregnancy_Risk": safety.get('pregnancy_lactation', 'UNKNOWN'),
                "Pediatric_Risk": safety.get('pediatrics', 'UNKNOWN'),
                "Filepath": os.path.basename(filepath)
            })
            
        except Exception as e:
            print(f"    Error reading {filepath}: {e}")

    # RANKING LOGIC: 
    # 1. Fewest Safety Flags (Safest)
    # 2. Lowest Toxicity Index (Cleanest)
    # 3. Lowest Potency nM (Strongest)
    leaderboard.sort(key=lambda x: (x['Safety_Flags'], x['Toxicity_Index'], x['Potency_nM']))
    
    # Write to CSV
    if leaderboard:
        keys = leaderboard[0].keys()
        with open(REPORT_FILE, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(leaderboard)
            
        print(f"\n>>> AUDIT COMPLETE.")
        print(f"    Top Candidate: {leaderboard[0]['Candidate_ID']}")
        print(f"    Safety Flags: {leaderboard[0]['Safety_Flags']}")
        print(f"    Potency: {leaderboard[0]['Potency_nM']} nM")
        print(f"    Report Generated: {REPORT_FILE}")
        
        print("\n    [TOP 3 RECOMMENDATIONS]")
        for i, cand in enumerate(leaderboard[:3]):
            print(f"    #{i+1}: {cand['Candidate_ID']} (Flags: {cand['Safety_Flags']}, Tox: {cand['Toxicity_Index']})")

if __name__ == "__main__":
    rank_diamonds()
