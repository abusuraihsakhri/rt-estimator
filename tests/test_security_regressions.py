"""Regression tests for integrity checks and CSV boolean parsing."""
import pytest
from agents.base import AuditTrail, SecurityException
from cli import parse_csv_bool

def test_signed_audit_chain_and_defensive_copy():
    audit = AuditTrail(secret_key="controlled-test-key")
    audit.log("tester", "worker", "complete", {"status": "ok"})
    audit.log("tester", "worker", "next", {"status": "ok"})
    assert audit.verify_integrity()
    external = audit.get_trail()
    external[0]["actor"] = "tampered"
    assert audit.verify_integrity()
    audit.logs[0]["actor"] = "tampered"
    assert not audit.verify_integrity()

def test_last_block_signature_tampering_detected():
    audit = AuditTrail(secret_key="controlled-test-key")
    audit.log("tester", "worker", "complete", {"status": "ok"})
    audit.logs[0]["current_hash"] = "0" * 64
    assert not audit.verify_integrity()

def test_actor_phi_does_not_enter_audit():
    audit = AuditTrail(secret_key="controlled-test-key")
    with pytest.raises(SecurityException):
        audit.log("patient-123-45-6789", "worker", "complete", {"status": "ok"})
    assert audit.get_trail() == []

@pytest.mark.parametrize("value,expected", [("False", False), ("0", False), ("no", False), ("TRUE", True), ("1", True), ("yes", True)])
def test_batch_boolean(value, expected):
    assert parse_csv_bool(value) is expected

def test_invalid_batch_boolean():
    with pytest.raises(ValueError):
        parse_csv_bool("perhaps")
