from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.apa_b import _validate_fit_audit, _validate_na_audit


FIELDS = [
    "status", "screening_tests", "screening_na", "screening_na_fraction",
    "confirmation_tests", "confirmation_na", "confirmation_na_fraction",
    "excluded_genes_fewer_than_two_testable_sites", "stageR_input_genes",
    "stageR_input_sites", "stageR_adjusted_na", "stageR_adjusted_na_fraction",
    "na_policy",
]


def write_audit(path: Path, **overrides: str) -> None:
    row = {
        "status": "WARN_UNTESTABLE_PVALUES",
        "screening_tests": "100",
        "screening_na": "2",
        "screening_na_fraction": "0.02",
        "confirmation_tests": "500",
        "confirmation_na": "5",
        "confirmation_na_fraction": "0.01",
        "excluded_genes_fewer_than_two_testable_sites": "3",
        "stageR_input_genes": "95",
        "stageR_input_sites": "490",
        "stageR_adjusted_na": "300",
        "stageR_adjusted_na_fraction": "0.6",
        "na_policy": "untestable hypotheses remain NA and cannot be significant",
    }
    row.update(overrides)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


class ApaBStatisticsTests(unittest.TestCase):
    def test_valid_na_audit_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "contrast.na_audit.tsv"
            write_audit(path)
            _validate_na_audit(path)

    def test_audit_rejects_zero_testable_hypotheses(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "contrast.na_audit.tsv"
            write_audit(path, stageR_input_genes="0")
            with self.assertRaisesRegex(RuntimeError, "no testable stageR hypotheses"):
                _validate_na_audit(path)

    def test_audit_rejects_policy_that_could_promote_na(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "contrast.na_audit.tsv"
            write_audit(path, na_policy="replace NA with zero")
            with self.assertRaisesRegex(RuntimeError, "unsupported NA policy"):
                _validate_na_audit(path)

    def test_audit_rejects_inconsistent_na_fraction(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "contrast.na_audit.tsv"
            write_audit(path, screening_na_fraction="0.5")
            with self.assertRaisesRegex(RuntimeError, "inconsistent fractions"):
                _validate_na_audit(path)

    def test_r_module_uses_supported_allow_na_policy(self) -> None:
        script = Path(__file__).parents[1] / "scripts" / "R" / "drimseq_stager_all_pairs.R"
        text = script.read_text(encoding="utf-8")
        self.assertIn('stageWiseAdjustment(staged, method="dtu", alpha=alpha, allowNA=TRUE)', text)
        self.assertIn('na_policy="untestable hypotheses remain NA and cannot be significant"', text)
        self.assertIn('if ("--self-test" %in% args)', text)

    def test_multifactor_numeric_retry_is_scoped_and_reproducible(self) -> None:
        script = Path(__file__).parents[1] / "scripts" / "R" / "drimseq_stager_all_pairs.R"
        text = script.read_text(encoding="utf-8")
        self.assertIn("fit_drimseq_with_numeric_retry", text)
        self.assertIn("one_way <- !isTRUE(multifactor)", text)
        self.assertIn("set.seed(seed)", text)
        self.assertIn("add_uniform=add_uniform", text)
        self.assertIn('policy <- "deterministic_add_uniform_retry"', text)
        self.assertIn('if (!isTRUE(multifactor) || !recognized_numeric_failure)', text)
        self.assertIn('paste0(con$contrast_id, ".fit_audit.tsv")', text)

    def test_valid_numeric_retry_fit_audit_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "contrast.fit_audit.tsv"
            path.write_text(
                "contrast_id\tstatus\tfit_policy\tmultifactor\tone_way\trandom_seed\t"
                "add_uniform_used\tprimary_error\n"
                "treated_vs_control\tWARN_NUMERIC_RETRY\t"
                "deterministic_add_uniform_retry\tTRUE\tFALSE\t123\tTRUE\t"
                "non-finite value supplied by optim\n",
                encoding="utf-8",
            )
            _validate_fit_audit(path)

    def test_retry_fit_audit_cannot_claim_unpaired_design(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "contrast.fit_audit.tsv"
            path.write_text(
                "contrast_id\tstatus\tfit_policy\tmultifactor\tone_way\trandom_seed\t"
                "add_uniform_used\tprimary_error\n"
                "treated_vs_control\tWARN_NUMERIC_RETRY\t"
                "deterministic_add_uniform_retry\tFALSE\tTRUE\t123\tTRUE\t"
                "non-finite value supplied by optim\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(RuntimeError, "numerical-retry audit is incomplete"):
                _validate_fit_audit(path)


if __name__ == "__main__":
    unittest.main()
