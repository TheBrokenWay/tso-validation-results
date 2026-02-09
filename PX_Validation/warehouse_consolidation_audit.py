"""
warehouse_consolidation_audit.py
Comprehensive Warehouse Audit for Full Consolidation

Scans all PX_Warehouse folders, identifies all assets, flags issues,
and generates detailed audit logs.
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import hashlib

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

class WarehouseAuditor:
    """Comprehensive warehouse auditor"""
    
    def __init__(self, warehouse_root="PX_Warehouse", log_dir="PX_LOGS"):
        self.warehouse_root = Path(warehouse_root)
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # Audit results
        self.findings = {
            "json_dossiers": [],
            "python_scripts": [],
            "markdown_docs": [],
            "worldlines": [],
            "batch_orders": [],
            "duplicates": [],
            "orphaned": [],
            "invalid_json": [],
            "folder_structure": {},
        }
        
        self.file_hashes = defaultdict(list)
        self.json_ids = defaultdict(list)
        
    def compute_hash(self, filepath):
        """Compute SHA256 hash of file"""
        try:
            with open(filepath, 'rb') as f:
                return hashlib.sha256(f.read()).hexdigest()
        except Exception as e:
            return f"ERROR: {e}"
    
    def scan_folder(self, folder, category="UNKNOWN"):
        """Recursively scan folder and categorize files"""
        if not folder.exists():
            return
        
        print(f"📂 Scanning: {folder.relative_to(self.warehouse_root)}")
        
        file_count = {"json": 0, "py": 0, "md": 0, "other": 0}
        
        for item in folder.rglob("*"):
            if item.is_file() and not item.name.startswith('.'):
                rel_path = item.relative_to(self.warehouse_root)
                file_info = {
                    "path": str(rel_path),
                    "size": item.stat().st_size,
                    "modified": datetime.fromtimestamp(item.stat().st_mtime).isoformat(),
                    "category": category,
                    "parent_folder": folder.name,
                }
                
                # Categorize by extension
                if item.suffix == ".json":
                    file_count["json"] += 1
                    # Try to parse JSON
                    try:
                        with open(item, 'r') as f:
                            data = json.load(f)
                        
                        # Extract ID if exists
                        json_id = data.get("id") or data.get("compound_id") or data.get("dossier_id")
                        if json_id:
                            self.json_ids[json_id].append(str(rel_path))
                            file_info["json_id"] = json_id
                        
                        file_info["json_valid"] = True
                        file_info["json_keys"] = list(data.keys())[:10]  # First 10 keys
                        
                        # Categorize JSON type
                        if "TRIAL_SIMULATION_DOSSIER" in item.name:
                            file_info["json_type"] = "TRIAL_SIMULATION"
                        elif "SMART" in item.name:
                            file_info["json_type"] = "SMART_ANTIVIRAL"
                        elif "BATCH" in item.name:
                            file_info["json_type"] = "BATCH_ORDER"
                        elif "DOSSIER" in item.name:
                            file_info["json_type"] = "COMMERCIAL_DOSSIER"
                        elif "WL-PRV" in item.name:
                            file_info["json_type"] = "WORLDLINE_CANDIDATE"
                        else:
                            file_info["json_type"] = "UNKNOWN"
                        
                        self.findings["json_dossiers"].append(file_info)
                        
                    except json.JSONDecodeError as e:
                        file_count["json"] += 1
                        file_info["json_valid"] = False
                        file_info["error"] = str(e)
                        self.findings["invalid_json"].append(file_info)
                
                elif item.suffix == ".py":
                    file_count["py"] += 1
                    file_info["type"] = "PYTHON_SCRIPT"
                    self.findings["python_scripts"].append(file_info)
                
                elif item.suffix == ".md":
                    file_count["md"] += 1
                    file_info["type"] = "MARKDOWN_DOC"
                    self.findings["markdown_docs"].append(file_info)
                
                elif item.suffix == ".worldline":
                    file_count["other"] += 1
                    file_info["type"] = "WORLDLINE"
                    self.findings["worldlines"].append(file_info)
                
                else:
                    file_count["other"] += 1
                
                # Check for duplicates by hash
                file_hash = self.compute_hash(item)
                if not file_hash.startswith("ERROR"):
                    self.file_hashes[file_hash].append(str(rel_path))
        
        return file_count
    
    def identify_duplicates(self):
        """Identify duplicate files by hash"""
        for file_hash, paths in self.file_hashes.items():
            if len(paths) > 1:
                self.findings["duplicates"].append({
                    "hash": file_hash,
                    "count": len(paths),
                    "paths": paths,
                })
        
        for json_id, paths in self.json_ids.items():
            if len(paths) > 1:
                self.findings["duplicates"].append({
                    "type": "JSON_ID_DUPLICATE",
                    "id": json_id,
                    "count": len(paths),
                    "paths": paths,
                })
    
    def identify_orphans(self):
        """Identify orphaned files in TrialSimulations root"""
        trial_sims_root = self.warehouse_root / "TrialSimulations"
        if not trial_sims_root.exists():
            return
        
        for item in trial_sims_root.iterdir():
            if item.is_file() and item.suffix == ".json":
                self.findings["orphaned"].append({
                    "path": str(item.relative_to(self.warehouse_root)),
                    "reason": "Loose file in TrialSimulations root (should be in LiveRuns/)",
                })
    
    def run_audit(self):
        """Run complete warehouse audit"""
        print("="*70)
        print("🔍 WAREHOUSE CONSOLIDATION AUDIT")
        print("="*70)
        print(f"Started: {datetime.now().isoformat()}")
        print()
        
        # Scan each top-level folder
        folders = [
            ("99_WAREHOUSE_ARCHIVE", "ARCHIVE"),
            ("00_COMMERCIAL_DOSSIERS", "COMMERCIAL"),
            ("00_LIVE_RESEARCH_OUTPUT", "RESEARCH"),
            ("Commercial_Dossiers", "COMMERCIAL"),
            ("SMART_Antiviral_Dossiers", "SMART_ANTIVIRAL"),
            ("TrialSimulations", "TRIAL_SIMULATION"),
            ("Orders", "BATCH_ORDER"),
            ("WorldLines", "WORLDLINE"),
        ]
        
        for folder_name, category in folders:
            folder_path = self.warehouse_root / folder_name
            if folder_path.exists():
                counts = self.scan_folder(folder_path, category)
                self.findings["folder_structure"][folder_name] = counts
        
        # Identify issues
        print("\n🔍 Identifying duplicates...")
        self.identify_duplicates()
        
        print("🔍 Identifying orphaned files...")
        self.identify_orphans()
        
        # Generate report
        self.generate_report()
        
        print("\n✅ Audit complete!")
    
    def generate_report(self):
        """Generate comprehensive audit report"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save detailed JSON log
        json_log_path = self.log_dir / f"consolidation_audit_{timestamp}.json"
        with open(json_log_path, 'w') as f:
            json.dump(self.findings, f, indent=2)
        
        # Generate human-readable summary
        summary_path = self.log_dir / f"consolidation_audit_{timestamp}.txt"
        with open(summary_path, 'w') as f:
            f.write("="*70 + "\n")
            f.write("WAREHOUSE CONSOLIDATION AUDIT REPORT\n")
            f.write("="*70 + "\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write("\n")
            
            # Summary counts
            f.write("SUMMARY COUNTS:\n")
            f.write("-"*70 + "\n")
            f.write(f"Total JSON Dossiers:       {len(self.findings['json_dossiers'])}\n")
            f.write(f"Total Python Scripts:      {len(self.findings['python_scripts'])}\n")
            f.write(f"Total Markdown Docs:       {len(self.findings['markdown_docs'])}\n")
            f.write(f"Total WorldLines:          {len(self.findings['worldlines'])}\n")
            f.write(f"Invalid JSON Files:        {len(self.findings['invalid_json'])}\n")
            f.write(f"Duplicate Files:           {len(self.findings['duplicates'])}\n")
            f.write(f"Orphaned Files:            {len(self.findings['orphaned'])}\n")
            f.write("\n")
            
            # Folder breakdown
            f.write("FOLDER STRUCTURE:\n")
            f.write("-"*70 + "\n")
            for folder, counts in self.findings["folder_structure"].items():
                f.write(f"\n{folder}:\n")
                f.write(f"  JSON:     {counts['json']}\n")
                f.write(f"  Python:   {counts['py']}\n")
                f.write(f"  Markdown: {counts['md']}\n")
                f.write(f"  Other:    {counts['other']}\n")
            
            # JSON type breakdown
            f.write("\n")
            f.write("JSON DOSSIER TYPES:\n")
            f.write("-"*70 + "\n")
            json_types = defaultdict(int)
            for dossier in self.findings["json_dossiers"]:
                json_type = dossier.get("json_type", "UNKNOWN")
                json_types[json_type] += 1
            for json_type, count in sorted(json_types.items()):
                f.write(f"  {json_type:30s}: {count}\n")
            
            # Duplicates
            if self.findings["duplicates"]:
                f.write("\n")
                f.write("DUPLICATES FOUND:\n")
                f.write("-"*70 + "\n")
                for dup in self.findings["duplicates"]:
                    if "hash" in dup:
                        f.write(f"Hash: {dup['hash'][:16]}... ({dup['count']} copies)\n")
                    else:
                        f.write(f"ID: {dup['id']} ({dup['count']} copies)\n")
                    for path in dup['paths']:
                        f.write(f"  - {path}\n")
                    f.write("\n")
            
            # Orphaned files
            if self.findings["orphaned"]:
                f.write("\n")
                f.write("ORPHANED FILES:\n")
                f.write("-"*70 + "\n")
                for orphan in self.findings["orphaned"]:
                    f.write(f"  {orphan['path']}\n")
                    f.write(f"    Reason: {orphan['reason']}\n")
            
            # Invalid JSON
            if self.findings["invalid_json"]:
                f.write("\n")
                f.write("INVALID JSON FILES:\n")
                f.write("-"*70 + "\n")
                for invalid in self.findings["invalid_json"]:
                    f.write(f"  {invalid['path']}\n")
                    f.write(f"    Error: {invalid['error']}\n")
        
        print(f"\n📊 Audit logs saved:")
        print(f"   JSON: {json_log_path}")
        print(f"   Text: {summary_path}")
        
        # Print summary to console
        print("\n" + "="*70)
        print("AUDIT SUMMARY")
        print("="*70)
        print(f"Total JSON Dossiers:       {len(self.findings['json_dossiers'])}")
        print(f"Total Python Scripts:      {len(self.findings['python_scripts'])}")
        print(f"Total Markdown Docs:       {len(self.findings['markdown_docs'])}")
        print(f"Invalid JSON Files:        {len(self.findings['invalid_json'])}")
        print(f"Duplicate Files:           {len(self.findings['duplicates'])}")
        print(f"Orphaned Files:            {len(self.findings['orphaned'])}")
        print("="*70)


if __name__ == "__main__":
    auditor = WarehouseAuditor()
    auditor.run_audit()
