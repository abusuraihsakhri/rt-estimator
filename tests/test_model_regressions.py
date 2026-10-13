"""Validate API audit inputs and explicitly unsupported provider behavior."""
import pytest
from pydantic import ValidationError
from agents.models import SystemTaskPayload
from agents.llm_factory import LLMFactory, MockLLM


def test_no_nan_or_infinite_api_metrics():
    payload = {"task_id": "T1", "target_identifier": "KEY-1", "primary_metric": 1.5}
    for bad in (float("nan"), float("inf"), float("-inf")):
        with pytest.raises(ValidationError):
            SystemTaskPayload(**{**payload, "primary_metric": bad})
        with pytest.raises(ValidationError):
            SystemTaskPayload(**{**payload, "secondary_metric": bad})


def test_model_factory_never_silently_implies_external_inference():
    assert isinstance(LLMFactory.create("mock"), MockLLM)
    for provider in ("openai", "ollama", "claude"):
        with pytest.raises(ValueError):
            LLMFactory.create(provider)
    assert "No model inference was performed" in MockLLM().invoke("test")
