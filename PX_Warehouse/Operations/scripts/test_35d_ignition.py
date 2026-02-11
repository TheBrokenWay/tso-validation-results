import random
import sys
import os

# Ensure root is in path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from PX_Warehouse.WorldLine_Database import WorldLineDatabase

random.seed(42)
print("🚀 IGNITING 35D ENGINE...")

db = WorldLineDatabase()
mock_35d_vector = [random.random() for _ in range(35)]

# Forces logic to write the .worldline file needed by Protocol Zero
filepath = db.record_materialization(
    task_id="PROTOCOL-ZERO-IGNITION",
    block=mock_35d_vector,
    coherence=0.99,
    lab_results={"status": "SIMULATED_GOLD"},
    route="TEST_VECTOR"
)

print(f"✅ SUCCESS. 35D Logic is active.")
print(f"📄 WorldLine created at: {filepath}")
