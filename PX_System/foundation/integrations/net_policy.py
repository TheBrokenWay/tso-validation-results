"""
Network Policy Utilities
Manages network access policies and rate limiting for external APIs
"""


class NetworkPolicyError(Exception):
    """Exception for network policy violations"""
    pass


def check_external_access_allowed(service: str) -> bool:
    """
    Check if external service access is allowed
    
    Args:
        service: Service name (e.g., "chembl", "pubchem", "uspto")
    
    Returns:
        bool: True if access is allowed
    """
    # Default: allow all external access
    # This can be enhanced with configuration-based policies
    return True


def get_rate_limit(service: str) -> dict:
    """
    Get rate limit configuration for service
    
    Args:
        service: Service name
    
    Returns:
        dict: Rate limit configuration
    """
    rate_limits = {
        "chembl": {"requests_per_second": 10, "requests_per_hour": 1000},
        "pubchem": {"requests_per_second": 5, "requests_per_hour": 500},
        "uspto": {"requests_per_second": 1, "requests_per_hour": 100}
    }
    
    return rate_limits.get(service, {"requests_per_second": 10, "requests_per_hour": 1000})


# Global configuration
_ONLINE_SYNC_ENABLED = True  # Enable external API access by default


def online_sync_enabled() -> bool:
    """
    Return current network policy for external API access.
    """
    return bool(_ONLINE_SYNC_ENABLED)


def set_online_sync_enabled(value: bool) -> None:
    """
    Update network policy for external API access.
    """
    global _ONLINE_SYNC_ENABLED
    _ONLINE_SYNC_ENABLED = bool(value)


__all__ = [
    "NetworkPolicyError", 
    "check_external_access_allowed", 
    "get_rate_limit",
    "online_sync_enabled",
    "set_online_sync_enabled"
]
