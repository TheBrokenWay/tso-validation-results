import time

MURAL = {"nodes": {}, "edges": []}

def update_node(role, version, fingerprint, health="OK"):
    MURAL["nodes"][role] = {
        "version": version,
        "fingerprint": fingerprint,
        "last_seen": time.time(),
        "health": health
    }

def update_edge(source, target, latency, success=True):
    MURAL["edges"].append({
        "source": source,
        "target": target,
        "last_activity": time.time(),
        "latency": latency,
        "success": success
    })

def get_mural():
    return MURAL
