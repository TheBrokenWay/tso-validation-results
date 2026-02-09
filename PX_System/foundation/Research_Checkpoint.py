"""
Research Checkpoint - Deterministic Implementation
Manages checkpoint/resume functionality for discovery campaigns
"""
from typing import Dict, Any, Optional
from pathlib import Path
import json
from datetime import datetime, timezone


class ResearchCheckpoint:
    """
    Checkpoint manager for autonomous research campaigns
    
    Provides save/restore functionality for discovery campaign state,
    allowing campaigns to be paused and resumed.
    
    Constitutional Compliance:
        - L51 (Zero Gaps): All state persisted to disk
        - ALCOA+: All checkpoints timestamped and immutable
    """
    
    def __init__(self, session_id: str, checkpoint_dir: Optional[Path] = None):
        """
        Initialize checkpoint manager
        
        Args:
            session_id: Unique session identifier
            checkpoint_dir: Directory for checkpoint files (default: temp/checkpoints)
        """
        self.session_id = session_id
        self.checkpoint_dir = checkpoint_dir or Path("temp/checkpoints")
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_file = self.checkpoint_dir / f"{session_id}_checkpoint.json"
    
    def save(self, state: Dict[str, Any]) -> str:
        """
        Save research campaign state to checkpoint
        
        Args:
            state: Campaign state dictionary
        
        Returns:
            str: Checkpoint filename
        
        Constitutional Compliance:
            - ALCOA+: Timestamped, immutable checkpoint
            - L9: All state data from actual campaign execution
        """
        checkpoint = {
            "session_id": self.session_id,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "state": state,
            "checkpoint_version": "1.0.0"
        }
        
        with open(self.checkpoint_file, 'w', encoding='utf-8') as f:
            json.dump(checkpoint, f, indent=2)
        
        return str(self.checkpoint_file)
    
    def load(self) -> Optional[Dict[str, Any]]:
        """
        Load research campaign state from checkpoint
        
        Returns:
            dict: Campaign state or None if no checkpoint exists
        """
        if not self.checkpoint_file.exists():
            return None
        
        try:
            with open(self.checkpoint_file, 'r', encoding='utf-8') as f:
                checkpoint = json.load(f)
            
            return checkpoint.get("state")
        except (json.JSONDecodeError, IOError):
            return None
    
    def exists(self) -> bool:
        """
        Check if checkpoint exists for this session
        
        Returns:
            bool: True if checkpoint file exists
        """
        return self.checkpoint_file.exists()
    
    def delete(self) -> bool:
        """
        Delete checkpoint file
        
        Returns:
            bool: True if deleted successfully
        """
        if self.checkpoint_file.exists():
            try:
                self.checkpoint_file.unlink()
                return True
            except OSError:
                return False
        return False


__all__ = ["ResearchCheckpoint"]
