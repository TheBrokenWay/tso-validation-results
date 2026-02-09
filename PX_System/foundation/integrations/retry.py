"""
Retry utilities for transient upstream errors
Used by engines when accessing external APIs (ChEMBL, PubChem, etc.)
"""


class TransientUpstreamError(Exception):
    """Exception for temporary upstream service failures"""
    def __init__(self, message: str, retry_after: int = 5):
        super().__init__(message)
        self.retry_after = retry_after


def retry_on_transient(func, max_attempts=3, delay=1.0):
    """
    Retry decorator for functions that may encounter transient errors
    
    Args:
        func: Function to retry
        max_attempts: Maximum retry attempts
        delay: Delay between retries (seconds)
    
    Returns:
        Function result or raises exception
    """
    import time
    
    for attempt in range(max_attempts):
        try:
            return func()
        except TransientUpstreamError as e:
            if attempt < max_attempts - 1:
                time.sleep(e.retry_after or delay)
            else:
                raise
        except Exception:
            raise


def retry_transient_http(max_attempts=3, delay=1.0, backoff=2.0):
    """
    Decorator for HTTP requests that may encounter transient errors
    
    Retries on common transient HTTP errors:
    - Connection errors
    - Timeout errors
    - 5xx server errors
    - 429 rate limiting
    
    Args:
        max_attempts: Maximum retry attempts (default: 3)
        delay: Initial delay between retries in seconds (default: 1.0)
        backoff: Backoff multiplier for exponential delay (default: 2.0)
    
    Returns:
        Decorator function
    
    Example:
        @retry_transient_http(max_attempts=5, delay=2.0)
        def fetch_chembl_data(compound_id):
            response = requests.get(f"https://chembl.org/api/{compound_id}")
            return response.json()
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            import time
            current_delay = delay
            
            for attempt in range(max_attempts):
                try:
                    return func(*args, **kwargs)
                except TransientUpstreamError as e:
                    if attempt < max_attempts - 1:
                        time.sleep(e.retry_after or current_delay)
                        current_delay *= backoff
                    else:
                        raise
                except Exception as e:
                    # Check if it's a retryable HTTP error
                    error_str = str(e).lower()
                    retryable = any(x in error_str for x in [
                        'connection', 'timeout', '5', '429', 'rate limit',
                        'service unavailable', 'gateway'
                    ])
                    
                    if retryable and attempt < max_attempts - 1:
                        time.sleep(current_delay)
                        current_delay *= backoff
                    else:
                        raise
            
            return None  # Should never reach here
        
        return wrapper
    return decorator


__all__ = ["TransientUpstreamError", "retry_on_transient", "retry_transient_http"]
