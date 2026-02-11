"""CSA (Constitutional) adapter for PX_System.benchmark. Delegates to PX_Executive.GAIP_Gateway."""
import asyncio
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))


class CSAPentarchy:
    def __init__(self):
        from PX_Executive.GAIP_Gateway import GAIPGateway
        self._gateway = GAIPGateway(mode="REGULATORY")

    async def evaluate_action(self, prop):
        # GAIPGateway.evaluate is sync; run in executor for async benchmark
        loop = asyncio.get_event_loop()
        csa = {"status": "COHERENT", "human_review_required": False}
        aas = {"status": "SUCCESS"}
        v_state = {"amplitude": 1.0}
        result = await loop.run_in_executor(
            None,
            lambda: self._gateway.evaluate(csa, aas, v_state, prop),
        )
        return {"status": "COHERENT" if result.get("authorized") else "REJECTED", "rationale": result.get("rationale", "")}
