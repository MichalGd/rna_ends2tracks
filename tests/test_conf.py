import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.conf import ConfError, project_from_conf, read_conf
from rnaends2tracks.config import ConfigError, build_conf_plan, workflow_requirements

HEADER = "sample_id,description,genome,biological_replicate_id,technical_replicate_id,lane_id,fastq_r1,fastq_r2,condition,batch,subject,library_protocol,library_layout,read_length,kit_catalog,umi_present\n"


class ConfTests(unittest.TestCase):
    def test_core_stage_requirements_follow_enabled_tracks(self):
        project = {
            "modules": {"gene_expression": False, "apa_a": False, "tracks": True},
            "apa_b": {"enabled": False},
            "tracks": {
                "families": {"all_reads": True, "exact_ends": False, "filtered_ends": False,
                             "rejected_ends": False, "active_pas": False},
                "normalizations": {"raw": True, "cpm": True, "deseq2": False, "robust_cpm": False},
            },
        }
        self.assertEqual(workflow_requirements(project), {
            "exact_ends": False, "active_pas": False, "apa_comparison": False,
        })
        project["tracks"]["families"]["filtered_ends"] = True
        self.assertTrue(workflow_requirements(project)["exact_ends"])
        self.assertFalse(workflow_requirements(project)["active_pas"])
        project["tracks"]["normalizations"]["deseq2"] = True
        self.assertTrue(workflow_requirements(project)["active_pas"])

    def test_parser_rejects_shell_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.conf"
            path.write_text("PROJECT_ID=x\nSAMPLESHEET=$(touch bad)\nOUTPUT_DIR=out\n", encoding="utf-8")
            with self.assertRaises(ConfError): read_conf(path)

    def test_ucsc_url_rejects_markdown_and_accepts_plain_http(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "config.conf"
            base = "PROJECT_ID=x\nSAMPLESHEET=samples.csv\nOUTPUT_DIR=results\n"
            config.write_text(
                base + "UCSC_BIGDATA_URL_PREFIX=[http://example.test](http://example.test)\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ConfError, "plain URL"):
                project_from_conf(config)
            config.write_text(
                base + "UCSC_BIGDATA_URL_PREFIX=http://example.test/project\n",
                encoding="utf-8",
            )
            project, _ = project_from_conf(config)
            self.assertEqual(
                project["tracks"]["ucsc_bigdata_url_prefix"],
                "http://example.test/project",
            )

    def test_enrichment_controls_and_apa_b_manifest_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); config = root / "config.conf"
            base = "PROJECT_ID=x\nSAMPLESHEET=samples.csv\nOUTPUT_DIR=results\n"
            config.write_text(base + "ENRICHMENT_MIN_GENESET_SIZE=600\nENRICHMENT_MAX_GENESET_SIZE=500\n",
                              encoding="utf-8")
            with self.assertRaisesRegex(ConfError, "must not exceed"):
                project_from_conf(config)
            config.write_text(base + "RUN_APA_B=true\nAPA_B_PILOT_ACCEPTED=true\nAPA_B_COMMAND_TEMPLATE=adapter\n",
                              encoding="utf-8")
            with self.assertRaisesRegex(ConfError, "APA_B_VALIDATION_MANIFEST"):
                project_from_conf(config)

    def test_apa_b_parallel_controls_are_bounded_by_engine_threads(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); config = root / "config.conf"
            base = "PROJECT_ID=x\nSAMPLESHEET=samples.csv\nOUTPUT_DIR=results\n"
            config.write_text(base + "APA_B_THREADS=8\nAPA_B_CLUSTER_PARALLEL_JOBS=9\n", encoding="utf-8")
            with self.assertRaisesRegex(ConfError, "APA_B_CLUSTER_PARALLEL_JOBS"):
                project_from_conf(config)
            config.write_text(
                base + "APA_B_THREADS=12\nAPA_B_ENDPOINT_SOURCE=exact_ends\n"
                "APA_B_ENDPOINT_PARALLEL_JOBS=6\nAPA_B_CLUSTER_PARALLEL_JOBS=8\nAPA_B_DEEPIP_THREADS=12\n",
                encoding="utf-8",
            )
            project, _ = project_from_conf(config)
            self.assertEqual(project["apa_b"]["endpoint_source"], "exact_ends")
            self.assertEqual(project["resources"]["apa_b"]["cluster_parallel_jobs"], 8)

    def test_downstream_and_method_specific_parallel_controls(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); config = root / "config.conf"
            base = "PROJECT_ID=x\nSAMPLESHEET=samples.csv\nOUTPUT_DIR=results\n"
            config.write_text(
                base + "PARALLEL_DOWNSTREAM_MODULES=true\nDOWNSTREAM_MODULE_PARALLEL_JOBS=2\n"
                "APA_CONTRAST_PARALLEL_JOBS=4\nAPA_A_CONTRAST_PARALLEL_JOBS=3\n",
                encoding="utf-8",
            )
            project, _ = project_from_conf(config)
            self.assertEqual(project["resources"]["downstream"]["parallel_modules"], 2)
            self.assertEqual(project["resources"]["apa_a"]["contrast_parallel_jobs"], 3)
            self.assertEqual(project["resources"]["apa_b"]["contrast_parallel_jobs"], 4)

    def test_mixed_genomes_create_only_within_genome_contrasts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sheet = root / "samples.csv"
            rows = []
            for genome, prefix in (("GRCh38", "H"), ("GRCm39", "M")):
                for condition in ("control", "treat"):
                    for replicate in (1, 2):
                        sample = f"{prefix}_{condition}_{replicate}"
                        rows.append(f"{sample},example,{genome},{sample},T01,L001,/data/{sample}.fastq.gz,,{condition},B1,{sample},quantseq_rev_v2_se,SE,101,REV_V2,false")
            sheet.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")
            config = root / "config.conf"
            config.write_text("\n".join([
                "PROJECT_ID=mixed", f'SAMPLESHEET="{sheet}"', f'OUTPUT_DIR="{root / "results"}"',
                "HG38_STAR_INDEX=/ref/h/star", "HG38_FASTA=/ref/h/fa", "HG38_GTF=/ref/h/gtf",
                "HG38_CHROM_SIZES=/ref/h/sizes", "HG38_PAS_ATLAS=/ref/h/pas",
                "MM39_STAR_INDEX=/ref/m/star", "MM39_FASTA=/ref/m/fa", "MM39_GTF=/ref/m/gtf",
                "MM39_CHROM_SIZES=/ref/m/sizes", "MM39_PAS_ATLAS=/ref/m/pas",
            ]) + "\n", encoding="utf-8")
            plan = build_conf_plan(config, check_inputs=False)
            self.assertTrue(plan.project["tracks"]["early_c0"])
            self.assertEqual(len(plan.contrasts), 2)
            self.assertEqual({row["genome"] for row in plan.contrasts}, {"GRCh38", "GRCm39"})
            self.assertTrue(all(row["contrast_id"].startswith(row["genome"] + ".") for row in plan.contrasts))

    def test_custom_minimum_replicates_is_respected(self):
        # Regression guard: this was previously hard-coded to two in contrast generation.
        with (
            self.assertRaisesRegex(ConfigError, "required biological replication"),
            tempfile.TemporaryDirectory() as temporary,
        ):
                root = Path(temporary); sheet = root / "samples.csv"
                rows = []
                for condition in ("control", "treat"):
                    for replicate in (1, 2):
                        sample = f"S_{condition}_{replicate}"
                        rows.append(f"{sample},x,GRCm39,{sample},T01,L001,/data/{sample}.fastq.gz,,{condition},B1,{sample},quantseq_rev_v2_se,SE,101,REV_V2,false")
                sheet.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")
                conf = root / "config.conf"
                conf.write_text("\n".join(["PROJECT_ID=x", f'SAMPLESHEET="{sheet}"', f'OUTPUT_DIR="{root / "out"}"',
                    "MIN_REPLICATES_PER_CONDITION=3", "MM39_STAR_INDEX=/r/s", "MM39_FASTA=/r/f",
                    "MM39_GTF=/r/g", "MM39_CHROM_SIZES=/r/c", "MM39_PAS_ATLAS=/r/p"]) + "\n", encoding="utf-8")
                build_conf_plan(conf, check_inputs=False)

    def test_condition_order_cannot_silently_drop_a_condition(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); sheet = root / "samples.csv"
            rows = []
            for condition in ("control", "treat", "recovery"):
                for replicate in (1, 2):
                    sample = f"S_{condition}_{replicate}"
                    rows.append(
                        f"{sample},x,GRCm39,{sample},T01,L001,/data/{sample}.fastq.gz,,"
                        f"{condition},B1,{sample},quantseq_rev_v2_se,SE,101,REV_V2,false"
                    )
            sheet.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")
            conf = root / "config.conf"
            conf.write_text("\n".join([
                "PROJECT_ID=x", f'SAMPLESHEET="{sheet}"', f'OUTPUT_DIR="{root / "out"}"',
                "CONDITION_ORDER=control,treat", "MM39_STAR_INDEX=/r/s", "MM39_FASTA=/r/f",
                "MM39_GTF=/r/g", "MM39_CHROM_SIZES=/r/c", "MM39_PAS_ATLAS=/r/p",
            ]) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ConfigError, "omits observed GRCm39 conditions: recovery"):
                build_conf_plan(conf, check_inputs=False)


if __name__ == "__main__":
    unittest.main()
