"""
SMILES Security Validation
Prevents injection attacks and validates SMILES string safety.

Constitutional Requirement: All molecular inputs must be validated
before processing to prevent adversarial attacks.
"""
import re


class SmilesSecurityError(Exception):
    """Exception for SMILES security violations"""
    pass


# Allowed SMILES characters (SMILES grammar)
SMILES_ALLOWED_CHARS = set("@+-=#$:()[]\\/%cCnNoOsSpPFIBr0123456789")

# Maximum reasonable SMILES length
MAX_SMILES_LENGTH = 2000


def validate_smiles_string(smiles: str) -> None:
    """
    Validate SMILES string for security
    
    Args:
        smiles: SMILES string to validate
    
    Raises:
        SmilesSecurityError: If SMILES contains invalid characters or patterns
    """
    if not isinstance(smiles, str):
        raise SmilesSecurityError("SMILES must be a string")
    
    if len(smiles) == 0:
        raise SmilesSecurityError("SMILES cannot be empty")
    
    if len(smiles) > MAX_SMILES_LENGTH:
        raise SmilesSecurityError(f"SMILES too long: {len(smiles)} > {MAX_SMILES_LENGTH}")
    
    # Check for shell injection patterns
    dangerous_patterns = [
        r'\$\(',  # $(command)
        r'`',     # backticks
        r';\s*\w+',  # semicolon command separator
        r'\|\s*\w+',  # pipe to command
        r'>\s*/\w+',  # redirect to file
        r'<\s*/\w+',  # redirect from file
    ]
    
    for pattern in dangerous_patterns:
        if re.search(pattern, smiles):
            raise SmilesSecurityError(f"Dangerous pattern detected: {pattern}")
    
    # Check for disallowed characters
    for char in smiles:
        if char not in SMILES_ALLOWED_CHARS:
            raise SmilesSecurityError(f"Invalid SMILES character: {char}")


def sanitize_smiles(smiles: str) -> str:
    """
    Sanitize SMILES string by removing dangerous characters
    
    Args:
        smiles: Raw SMILES string
    
    Returns:
        Sanitized SMILES string
    """
    # Remove any non-SMILES characters
    return ''.join(c for c in smiles if c in SMILES_ALLOWED_CHARS)


__all__ = ["SmilesSecurityError", "validate_smiles_string", "sanitize_smiles"]
