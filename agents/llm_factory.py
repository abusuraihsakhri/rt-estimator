"""A deterministic mock provider for offline testing; no remote model inference."""
from .base import PHIGuard


class MockLLM:
    def __init__(self, system_name: str = "Rt Estimator"):
        self.system_name = system_name

    def invoke(self, prompt: str) -> str:
        PHIGuard.assert_no_phi(prompt)
        return f"[{self.system_name} mock]: No model inference was performed. Query received: '{prompt[:60]}...'."


class LLMFactory:
    """Return a local deterministic mock; reject providers not implemented here."""

    @staticmethod
    def create(provider: str = "mock", system_name: str = "Rt Estimator"):
        prov = str(provider).lower()
        if prov in ("mock", "deterministic", "test"):
            return MockLLM(system_name)
        raise ValueError(f"model provider {provider!r} is not implemented")
