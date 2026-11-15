import csv
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.prepare_expression import prepare_project
from rnaends2tracks.conf import project_from_conf


class PrepareExpressionTests(unittest.TestCase):
    def test_separate_project_preserves_metadata_and_original_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); old=root/"old"; old.mkdir()
            header="sample_id,description,genome,biological_replicate_id,technical_replicate_id,lane_id,fastq_r1,fastq_r2,condition,batch,subject,library_protocol,library_layout,read_length,kit_catalog,umi_present\n"
            (old/"samples.csv").write_text(header+"".join(f"{condition}{rep},fixture,GRCm39,{condition}{rep},T01,L001,deleted_{condition}{rep}.fq.gz,,{condition},B1,{condition}{rep},quantseq_rev_v2_se,SE,150,REV_V2,false\n" for condition in ("A","B") for rep in (1,2)))
            source=old/"config.conf"
            content='PROJECT_ID=old\nSAMPLESHEET=samples.csv\nOUTPUT_DIR=results\nMM39_GTF=genes.gtf\nMM39_CHROM_SIZES=sizes\nPAIRING_MODE=none\n'
            source.write_text(content)
            metadata=old/"results"/"00_metadata"; metadata.mkdir(parents=True)
            (metadata/"contrasts.tsv").write_text("genome\tcontrast_id\tnumerator\tdenominator\tresolved_design\nGRCm39\tB_vs_A\tB\tA\t~ condition\n")
            target=prepare_project(source,root/"new","new",jobs=1)
            project,sheet=project_from_conf(target)
            self.assertEqual(project["input_mode"],"bam")
            self.assertEqual(Path(project["alignment_source_dir"]),(old/"results"/"02_alignment").resolve())
            self.assertEqual(source.read_text(),content)
            self.assertEqual(project["gene_expression"]["min_count"],5)
            self.assertFalse(project["modules"]["tracks"])
            with Path(sheet).open() as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(len(rows),4)
            self.assertEqual(Path(rows[0]["fastq_r1"]),(old/"deleted_A1.fq.gz").resolve())
            self.assertTrue((target.parent/"source_metadata"/"source_manifest.json").is_file())
            with self.assertRaisesRegex(ValueError,"empty"):
                prepare_project(source,root/"new","new",jobs=1)
