import hashlib
import json

VM_SPEC = {
    "name": "PREDATOR_X_SCIENTIFIC_VM",
    "version": "1.2.0-GAIP",
    "max_parallelism": 8,
    "memory_bandwidth_gbps": 64,
    "latency_model": "constant",
    "float_precision": "float64"
}

def get_vm_spec():
    return VM_SPEC

def get_vm_fingerprint():
    return hashlib.sha256(json.dumps(VM_SPEC, sort_keys=True).encode()).hexdigest()
