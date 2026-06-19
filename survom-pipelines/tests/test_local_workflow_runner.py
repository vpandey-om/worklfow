import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import local_workflow_runner as runner


class LocalWorkflowRunnerStateTests(unittest.TestCase):
    def test_parse_param_overrides(self):
        params = runner.parse_param_overrides(
            [
                "quality_threshold=25",
                "minimum_read_length=30",
                "trim_poly_g=false",
                "salmon_index=/tmp/salmon_index",
                "custom_value=null",
            ]
        )
        self.assertEqual(params["quality_threshold"], 25)
        self.assertEqual(params["minimum_read_length"], 30)
        self.assertEqual(params["trim_poly_g"], False)
        self.assertEqual(params["salmon_index"], "/tmp/salmon_index")
        self.assertIsNone(params["custom_value"])

    def test_state_save_load_and_completed_output_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            original_runs_root = runner.RUNS_ROOT
            runner.RUNS_ROOT = tmp_path / "runs"
            try:
                run_id = "demo_run"
                output_dir = tmp_path / "outputs"
                output_dir.mkdir()
                (output_dir / "result.txt").write_text("ok")

                state = runner.initial_state(run_id, "trim_only", output_dir)
                state["status"] = "completed"
                runner.save_state(state)

                loaded = runner.load_state(run_id)
                self.assertEqual(loaded["run_id"], run_id)
                self.assertEqual(loaded["workflow"], "trim_only")
                self.assertEqual(loaded["status"], "completed")
                self.assertTrue(runner.outputs_valid(loaded))
            finally:
                runner.RUNS_ROOT = original_runs_root

    def test_update_state_from_execution_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            original_runs_root = runner.RUNS_ROOT
            runner.RUNS_ROOT = tmp_path / "runs"
            try:
                state = runner.initial_state("resume_me", "qc_trim_strandedness", tmp_path / "outputs")
                runner.save_state(state)
                execution_record = tmp_path / "execution_record.json"
                execution_record.write_text(
                    json.dumps(
                        {
                            "step_outputs": [
                                {"step_id": "raw_qc", "process_name": "FASTQC", "status": "completed", "outputs": []},
                                {
                                    "step_id": "trim",
                                    "process_name": "FASTP",
                                    "status": "failed",
                                    "outputs": [],
                                    "missing_outputs": [],
                                },
                            ]
                        }
                    )
                )

                runner.update_state_from_execution(state, execution_record, "failed")

                saved = runner.load_state("resume_me")
                self.assertEqual(saved["status"], "failed")
                self.assertEqual(saved["completed_steps"], ["raw_qc"])
                self.assertEqual(saved["failed_step"], "trim")
                self.assertEqual(saved["current_step"], "trim")
            finally:
                runner.RUNS_ROOT = original_runs_root


if __name__ == "__main__":
    unittest.main()
