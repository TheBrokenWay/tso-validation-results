import time
import json
import os
from datetime import datetime

class SystemTelemetry:
    """
    Monitors pipeline throughput, error rates, execution drift, and hash mismatches.
    """
    def __init__(self, telemetry_path=r"E:\foundation\PX_Warehouse\Operations\system_telemetry.json"):
        self.telemetry_path = telemetry_path
        self.data = self._load_telemetry()

    def log_execution(self, pipeline_name, duration, items_processed, errors=0):
        """
        Logs a pipeline execution run.
        """
        if "pipelines" not in self.data:
            self.data["pipelines"] = {}
        
        if pipeline_name not in self.data["pipelines"]:
            self.data["pipelines"][pipeline_name] = []

        entry = {
            "timestamp": datetime.now().isoformat(),
            "duration_sec": duration,
            "throughput": items_processed / duration if duration > 0 else 0,
            "items_processed": items_processed,
            "errors": errors
        }
        self.data["pipelines"][pipeline_name].append(entry)
        self._save_telemetry()

    def log_drift(self, metric_name, expected_value, actual_value):
        """
        Logs execution drift for sensitive metrics.
        """
        if "drift" not in self.data:
            self.data["drift"] = {}
        
        if metric_name not in self.data["drift"]:
            self.data["drift"][metric_name] = []

        entry = {
            "timestamp": datetime.now().isoformat(),
            "expected": expected_value,
            "actual": actual_value,
            "variance": actual_value - expected_value
        }
        self.data["drift"][metric_name].append(entry)
        self._save_telemetry()

    def _load_telemetry(self):
        if os.path.exists(self.telemetry_path):
            with open(self.telemetry_path, 'r') as f:
                return json.load(f)
        return {"pipelines": {}, "drift": {}, "mismatches": []}

    def _save_telemetry(self):
        os.makedirs(os.path.dirname(self.telemetry_path), exist_ok=True)
        with open(self.telemetry_path, 'w') as f:
            json.dump(self.data, f, indent=4)

if __name__ == "__main__":
    telemetry = SystemTelemetry()
    # Example: Log a stratification run
    telemetry.log_execution("STRATIFICATION_V3", 45.2, 265, errors=2)
    print("Telemetry logged.")
