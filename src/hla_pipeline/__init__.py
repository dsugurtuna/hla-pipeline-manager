"""HLA Pipeline Manager: plan, verify and deploy SNP2HLA imputation batches."""

__version__ = "2.0.0"

from .deployer import DeploymentReport, ResultDeployer
from .executor import BatchExecutor, ExecutionPlan, ExecutionResult, PipelineConfig
from .reporter import ClinicalReport, ClinicalReporter
from .verifier import ImputationVerifier, VerificationReport

__all__ = [
    "BatchExecutor",
    "ExecutionResult",
    "ExecutionPlan",
    "PipelineConfig",
    "ImputationVerifier",
    "VerificationReport",
    "ResultDeployer",
    "DeploymentReport",
    "ClinicalReporter",
    "ClinicalReport",
]
